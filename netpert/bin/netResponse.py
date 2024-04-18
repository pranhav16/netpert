#!/usr/bin/env python

import time

import sys
import os.path
import argparse
import logging
import pandas as pd
import numpy as np
import warnings

from math import log, e, isinf, isnan, isclose
from scipy import linalg,integrate,stats
from collections import defaultdict

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

logging.basicConfig(format='%(name)s %(levelname)s: %(message)s')
logger = logging.getLogger('NetResponse')
logger.setLevel(logging.INFO)

current_path = os.path.abspath(__file__)

# Navigate to the parent directory (netpert/bin to netpert)
netpert_path = os.path.dirname(os.path.dirname(current_path))

# Combine with 'databases' to get the path to the databases directory
databases_path = os.path.join(netpert_path, 'databases')
projects_path = os.path.join(netpert_path, 'projects')

TOPLEVEL = '.'
# DBDIR = os.path.join(TOPLEVEL, 'databases')
DBDIR = databases_path
# PROJDIR = os.path.join(TOPLEVEL, 'projects')
PROJDIR = projects_path

class CMDParser():
    #Command line parser

    def __init__(self):
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()

        #NetResponse analysis
        analysis_parser = subparser.add_parser('analysis', help='Run NetResponse.')
        analysis_parser.add_argument('responses',help='Response file name')
        analysis_parser.add_argument('biologicalknowledge',help='Biological Knowledge File Name')
        analysis_parser.add_argument('driver',help='Network driver gene symbol.')
        analysis_parser.add_argument('-s',choices=['human','mouse'],default='human',help='species. Default human')
        analysis_parser.add_argument('-v',choices=['standard','reduced'],default='reduced',help='Network size parameter')
        analysis_parser.set_defaults(func = runAnalysis)

        args = parser.parse_args()
        self.__dict__.update(args.__dict__)
        return
   
def runAnalysis(args):
    """
    Runs NetResponse analysis
    """
    print('*************************************************')
    logger.info('Running NetResponse analysis...')
    print('*************************************************')
    start_time = time.time()
    getNetworkInteractions(args)
    end_time = time.time()
    logger.info('getNetworkInteractions took %.1f seconds',end_time-start_time)
    start_time = time.time()
    driver,intermediates,responses = getNetworkGenes(args,verbose=True)
    end_time = time.time()
    logger.info('getNetworkGenes took %.1f seconds',end_time-start_time)
    if args.v == 'reduced':
        intermediates = getSmallerNetworkInteractions(args,driver,intermediates,responses)
    start_time = time.time()
    getPrioritization(args,driver,intermediates,responses,endpoints=False,r=0.5,nt=8)
    end_time = time.time()
    logger.info('getPrioritization took %.1f seconds',end_time-start_time)
    return None

def getNetworkInteractions(args):
    """
    Compiles interactions into one file and stores it in database dir, if needed.
    Removes unwanted PPI between a protein and itself.
    """
    logger.info('Compiling %s network interactions from data files...',args.s)
    
    #get filepath
    if args.s == 'human':
        interaction_filepath = os.path.join(DBDIR, 'human_interactions.txt')
    elif args.s == 'mouse':
        interaction_filepath = os.path.join(DBDIR, 'mouse_interactions.txt')
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    #check if interaction file exists
    #if os.path.isfile(interaction_filepath):
    #    logger.info('Found compiled %s interactions file %s',args.s,interaction_filepath)
    #else:
        #read PPI and TF
    if args.s == 'human':
      ppi = readIIDHuman()
      tf = readTF('human')
    elif args.s == 'mouse':
      ppi = readIIDMouse()
      tf = readTF('mouse')
    else:
      logger.error('%s is not a supported species.', args.s)
      sys.exit()
        #combine data and write
    ppi_frame = pd.DataFrame(ppi,columns=['source','target'])
    ppi_frame['type'] = 'PPI'
    tf_frame = pd.DataFrame(tf,columns=['source','target'])
    tf_frame['type'] = 'TF'
    frames = [ppi_frame,tf_frame]
    wholeinteract = pd.concat(frames)
    wholeinteract.to_csv(interaction_filepath,sep='\t',index=False,header=False)
    logger.info('Wrote compiled %s interactions file to %s',args.s,interaction_filepath)
    
    if args.v == 'standard':
      edges = readEdges(args,verbose=True)
       
    return None
def getSmallerNetworkInteractions(args, driver, intermediates, responses):
    """
    Compiles interactions into one file and stores it in database dir, if needed.
    Removes unwanted PPI between a protein and itself.
    """
    logger.info('Compiling %s network interactions from data files...',args.s)
    
    #get filepath
    if args.s == 'human':
        interaction_filepath = os.path.join(DBDIR, 'human_interactions.txt')
    elif args.s == 'mouse':
        interaction_filepath = os.path.join(DBDIR, 'mouse_interactions.txt')
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    #check if interaction file exists
    #if os.path.isfile(interaction_filepath):
       # logger.info('Found compiled %s interactions file %s',args.s,interaction_filepath)
            #read PPI and TF
    if args.s == 'human':
        ppi = readIIDHuman()
        tf = readTF('human')
    elif args.s == 'mouse':
        ppi = readIIDMouse()
        tf = readTF('mouse')
    else:
        logger.error('%s is not a supported species.', args.s)
        sys.exit()
    #combine data and write
    ppi_frame = pd.DataFrame(ppi,columns=['source','target'])
    ppi_frame['type'] = 'PPI'
    tf_frame = pd.DataFrame(tf,columns=['source','target'])
    tf_frame['type'] = 'TF'
    frames = [ppi_frame,tf_frame]
    wholeinteract = pd.concat(frames)
    count = 0
    driverToInter = wholeinteract[(wholeinteract['source']== driver) & (wholeinteract['target'].isin(intermediates))]
    InterToResp = wholeinteract[(wholeinteract['source'].isin(intermediates)) & ( wholeinteract['target'].isin(responses))]
    #wholeinteract.to_csv(interaction_filepath,sep='\t',index=False,header=False)
    #driverToInter.to_csv(interaction_filepath,sep='\t',index=False,header=False)
    #InterToResp.to_csv('./networks/Twist1/interToResp.txt',sep='\t',index=False,header=False)
        
    smallerIntermediates = list()
    for inter in driverToInter['target']:
        if inter not in smallerIntermediates:
            smallerIntermediates.append(inter)
    for inter in InterToResp['source']:
        if inter not in smallerIntermediates:
            smallerIntermediates.append(inter)
    
    smallNetwork = list()
    smallNetwork.append(driver)
    smallNetwork.extend(smallerIntermediates)
    smallNetwork.extend(responses)
    
    
   # DriverToInterToResp = pd.concat([driverToInter,InterToResp])
   # InterToInter = wholeinteract[(wholeinteract['source'].isin(smallerIntermediates)) & (wholeinteract['target'].isin(smallerIntermediates))]
   # smallerWholeInteractions = pd.concat([DriverToInterToResp,InterToInter])
    smallerWholeInteractions = wholeinteract[(wholeinteract['source'].isin(smallNetwork)) & (wholeinteract['target'].isin(smallNetwork))]
    smallerWholeInteractions.to_csv(interaction_filepath,sep='\t',index=False,header=False)

    logger.info('Wrote compiled %s interactions file to %s',args.s,interaction_filepath)

    edges = readEdges(args,verbose=True)
    
    fp = open(os.path.join(PROJDIR,'/networkIntermediates.txt'), 'w')
    for intermediate in smallerIntermediates:
        fp.write(intermediate + '\n')
    fp.close()
    logger.info('Found %s intermediate genes in the network. Wrote to %s.' % (str(len(smallerIntermediates)),os.path.join(PROJDIR,'/networkIntermediates.txt')))
    return(smallerIntermediates)




def getNetworkGenes(args,verbose=False):
    """
    Reads user-defined driver and respones. Checks they are in the network.
    Identifies intermediates from interaction file.
    Returns driver, list of intermediates, and list of network respones.
    """
    if verbose:
        logger.info('Compiling network genes...')

    #Read network interactions
    if args.s == 'human':
        gene_filepath = os.path.join(DBDIR, 'human_genes.txt')
        interaction_filepath = os.path.join(DBDIR, 'human_interactions.txt')
    elif args.s == 'mouse':
        gene_filepath = os.path.join(DBDIR, 'mouse_genes.txt')
        interaction_filepath = os.path.join(DBDIR, 'mouse_interactions.txt')
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    if not os.path.isfile(gene_filepath):

        #read interaction file
        wholeinteract = readInteractions(interaction_filepath)

        #get whole gene list
        wholegenelist = wholeinteract[0].tolist() + wholeinteract[1].tolist()
        wholegenelist = list(set(wholegenelist))

        #write to file
        wholegenedf = pd.DataFrame(wholegenelist)
        wholegenedf.to_csv(gene_filepath,sep='\t',index=False,header=False)
        logger.info('Wrote compiled %s gene list file to %s', args.s, gene_filepath)
        
    #check if gene list is imported
    try:
        wholegenelist
    except NameError:
        wholegenelist = readGenes(gene_filepath)
        if verbose:
            logger.info('Found compiled %s gene list file %s',args.s,gene_filepath)

    if verbose:
        logger.info('Network directory: %s',PROJDIR)
    
    #read driver and check it is in network
    driver = checkDriver(args,wholegenelist,verbose)
    
    #get responses in network
    responses = getNetworkResponses(args,driver,wholegenelist,verbose)

    #get intermediates in network
    intermediatesFilepath = os.path.join(PROJDIR,'/networkIntermediates.txt') 
#    if os.path.isfile(intermediatesFilepath):
#        if verbose:
#            logger.info('Found network intermediates file %s. Delete this file if you want the intermediates to be compiled again',intermediatesFilepath)
 #   else:
        #get intermediates
    intermediates = wholegenelist

        #remove driver from intermediates
    intermediates.remove(driver)

        #remove responses from intermediates
    for r in responses:
        intermediates.remove(r)
    
        #output intermediates
    fp = open(intermediatesFilepath, 'w')
    for intermediate in intermediates:
        fp.write(intermediate + '\n')
    fp.close()
    if verbose:
        logger.info('Found %s intermediate genes in the network. Wrote to %s.' % (str(len(intermediates)),intermediatesFilepath))

    #check if intermediates is imported
    try:
        intermediates
    except NameError:
        intermediates = readGenes(intermediatesFilepath)

    #output network size
    if verbose:
        network_size = 1 + len(intermediates) + len(responses)
        logger.info('Network consists of %d genes: the driver %s, %d intermediates, and %d responses.',network_size,driver,len(intermediates),len(responses))

    return(driver,intermediates,responses)

def getPrioritization(args,driver,intermediates,responses,endpoints=False,r=0.5,nt=16):
    """
    Calculates NetResponse weights.
    Ranks network genes.
    Writes results
    """

    #get whole gene list
    wholegenelist = list()
    wholegenelist.append(driver)
    wholegenelist.extend(intermediates)
    wholegenelist.extend(responses)
    
    #get network edges 
    edge_dict = readEdges(args)
    edge_ab, edge_ba = edgeTuplesToDicts(edge_dict)

    #get genename (vertex) to index dict and vice versa
    v2i, i2v = nameToIndex(wholegenelist)

    #get genename (vertex) to indegree dict and outdegree dict
    v2indeg, v2outdeg = nameToDegree(wholegenelist,edge_ba,edge_ab)

    #get genename (vertex) to type (driver, intermediate, response)
    v2type = nameToType(wholegenelist,driver,intermediates,responses)

    #get laplacian matrix
    (adj_mat, lap_mat) = getLaplacian(edge_dict, v2i)

    #get time step size dt
    dt = getTimeStep(driver,v2outdeg,r=r,nt=nt)

    #get G(t) list
    g_list = getGlist(lap_mat,dt,nt)

    #get weights
    v2weight = getWeights(args,g_list,v2i,driver,responses,nt,dt,endpoints=endpoints)

    #get rankings
    (vertices_ranked,v2rank) = getRankings(v2weight,reverse_order=True)

    #get P-value
    v2pVal = readBiologicalKnowledge_pVal(args)

    #get log2FC
    v2log2FC = readBiologicalKnowledge_logFC(args)

    #map human gene to repurposing hub drug
    hgene2repurp,repurp2hgene = getRepurpDrugs()

    #map mouse gene to repurposing hub drug
    v2repurp = getMouseGene2RepurpDrugs(wholegenelist,hgene2repurp)

    #get gene pathway type
    v2pathType = getPathType(args,driver,intermediates,responses)

    #write to file
    writeWeights(args,v2weight,v2rank,v2pathType,v2repurp,repurp2hgene,vertices_ranked, v2pVal,v2log2FC, drug_target_max=5)

    return None

def getPathType(args,driver,intermediates,responses):
    """
    Returns pathway type
    D = driver
    R = response
    DIR = intermediate with edge from driver and edge to response
    DI = intermediate with edge from driver
    IR = intermediate with edge to response
    I = intermediate no direct edges to/from driver/response
    """

    #get whole gene list
    wholegenelist = list()
    wholegenelist.append(driver)
    wholegenelist.extend(intermediates)
    wholegenelist.extend(responses)

    #get type (D,I,R)
    v2type = nameToType(wholegenelist,driver,intermediates,responses)

    interm_dict = dict()
    resp_dict = dict()
    driv_dict = dict()
    v2pathType = dict()

    for gene in wholegenelist:
        typ = v2type[gene]
        if typ == 'I':
            interm_dict[gene] = True
        elif typ == 'R':
            resp_dict[gene] = True
        elif typ == 'D' or typ == 'DR':
            driv_dict[gene] = True
        else:
            print('error in type')

    #read interactions
    if args.s == 'mouse':
        interaction_filepath = os.path.join(DBDIR, 'mouse_interactions.txt')
    elif args.s == 'human':
        interaction_filepath = os.path.join(DBDIR, 'human_interactions.txt')
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    reftable = pd.read_csv(interaction_filepath,sep='\t',header=None,names=['source','target','type'])

    indegFromDriv = defaultdict(dict)
    indegFromInterm = defaultdict(dict)
    outdeg2interm = defaultdict(dict)
    outdeg2resp = defaultdict(dict)

    for index,row in reftable.iterrows():
        src = row['source']
        des = row['target']
        typ = row['type']

        if typ == 'PPI':
            #3 gene types, 2 types per interaction, resampling permitted, order matters, so 9 combinations
            #interm <-> resp
            if src in interm_dict and des in resp_dict:
                #src
                outdeg2resp[src][(src,des)] = True
                #des
                indegFromInterm[des][(src,des)] = True
                outdeg2interm[des][(des,src)] = True
            #resp <-> interm
            elif src in resp_dict and des in interm_dict:
                #src
                indegFromInterm[src][(des,src)] = True
                outdeg2interm[src][(src,des)] = True
                #des
                outdeg2resp[des][(des,src)] = True
            #interm <-> interm
            elif src in interm_dict and des in interm_dict:
                #src
                indegFromInterm[src][(des,src)] = True
                outdeg2interm[src][(src,des)] = True
                #des
                indegFromInterm[des][(src,des)] = True
                outdeg2interm[des][(des,src)] = True
            #driv <-> interm
            elif src in driv_dict and des in interm_dict:
                #src
                indegFromInterm[src][(des,src)] = True
                outdeg2interm[src][(src,des)] = True
                #des
                indegFromDriv[des][(src,des)] = True
            #interm <-> driv
            elif src in interm_dict and des in driv_dict:
                #src
                indegFromDriv[src][(des,src)] = True
                #des
                indegFromInterm[des][(src,des)] = True
                outdeg2interm[des][(des,src)] = True
            #driv <-> resp
            elif src in driv_dict and des in resp_dict:
                #src
                outdeg2resp[src][(src,des)] = True
                #des
                indegFromDriv[des][(src,des)] = True
            #resp <-> driv
            elif src in resp_dict and des in driv_dict:
                #src
                indegFromDriv[src][(des,src)] = True
                #des
                outdeg2resp[des][(des,src)] = True
            #driv <-> driv
            elif src in driv_dict and des in driv_dict:
                #des
                indegFromDriv[des][(src,des)] = True
            #resp <-> resp
            elif src in resp_dict and des in resp_dict:
                #src
                outdeg2resp[src][(src,des)] = True
                #des
                outdeg2resp[des][(des,src)] = True
        elif typ == 'TF':
            #3 gene types, 2 types per interaction, resampling permitted, order matters, so 9 combinations
            #interm -> resp
            if src in interm_dict and des in resp_dict:
                #src
                outdeg2resp[src][(src,des)] = True
                #des
                indegFromInterm[des][(src,des)] = True
            #resp -> interm
            elif src in resp_dict and des in interm_dict:
                #src
                outdeg2interm[src][(src,des)] = True
            #interm -> interm
            elif src in interm_dict and des in interm_dict:
                #src
                outdeg2interm[src][(src,des)] = True
                #des
                indegFromInterm[des][(src,des)] = True
            #driv -> interm 
            elif src in driv_dict and des in interm_dict:
                #src
                outdeg2interm[src][(src,des)] = True
                #des
                indegFromDriv[des][(src,des)] = True
            #interm -> driv
            elif src in interm_dict and des in driv_dict:
                #des
                indegFromInterm[des][(src,des)] = True
            #driv -> resp
            elif src in driv_dict and des in resp_dict:
                #src
                outdeg2resp[src][(src,des)] = True
                #des
                indegFromDriv[des][(src,des)] = True
            #resp -> driv
                #not needed
            #resp -> resp
            elif src in resp_dict and des in resp_dict:
                #src
                outdeg2resp[src][(src,des)] = True
            #driv -> driv
            elif src in driv_dict and des in driv_dict:
                #des
                indegFromDriv[des][(src,des)] = True
        else:
            print('error in interaction')

    v2indegDriv = dict()
    v2indegInterm = dict()
    v2outdegInterm = dict()
    v2outdegResp = dict()

    for v in wholegenelist:

        indeg_driv = str(len(indegFromDriv[v]))
        indeg_interm = str(len(indegFromInterm[v]))
        outdeg_interm = str(len(outdeg2interm[v]))
        outdeg_resp = str(len(outdeg2resp[v]))

        v2indegDriv[v] = indeg_driv
        v2indegInterm[v] = indeg_interm
        v2outdegInterm[v] = outdeg_interm
        v2outdegResp[v] = outdeg_resp

        if v in interm_dict:
            if indeg_driv == "0" and outdeg_resp == "0":
                typ = 'I'
            elif indeg_driv == "0":
                typ = 'IR'
            elif outdeg_resp == "0":
                typ = 'DI'
            else:
                typ = 'DIR'
        elif v in driv_dict:
            typ = 'D'
        elif v in resp_dict:
            typ = 'R'
        else:
            print('error bad type')
        v2pathType[v] = typ

    return(v2pathType)

def readBiologicalKnowledge(args):
    """
    Reads file called biologicalknowledge.
    File may have extension .xlsx, .csv, or .tsv.
    File must have columns named Gene and Quanity.
    Gene is a column of gene symbols.
    Quantity is a column of numbers used to rank genes.
    Returns dict of gene to quantity.
    """
    v2knowledge = dict()

    #find biological knowledge filepath
    knowledgeFilepath = args.biologicalknowledge
    if os.path.isfile(knowledgeFilepath):
        #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',knowledgeFilepath)
        pass
    else:
        knowledgeFilepath = args.biologicalknowledge
        if os.path.isfile(knowledgeFilepath): 
            #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .tsv extension.',knowledgeFilepath)
            pass
        else:
            knowledgeFilepath = args.biologicalknowledge
            if os.path.isfile(knowledgeFilepath): 
                #logger.info('Found biological knowledge file %s.',knowledgeFilepath)
                pass
            else:
                logger.error('No biological knowledge file found. Must provide a file named biologicalknowledge with .xlsx, .csv, or .tsv extension in directory %s.', PROJDIR)
                sys.exit()

    #check if file is empty
    if not os.stat(knowledgeFilepath).st_size:
        logger.error('File %s is empty.',knowledgeFilepath)
        sys.exit()

    #read data
    if knowledgeFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(knowledgeFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Biological knowledge file %s must contain only 1 sheet.', knowledgeFilepath)
            sys.exit()
        df = pd.read_excel(knowledgeFilepath)
    elif knowledgeFilepath.split(".")[-1] == 'csv':
        df = pd.read_csv(knowledgeFilepath, sep=',', engine='python')
    else:
        df = pd.read_csv(knowledgeFilepath, sep='\t', engine='python')

    #check for columns
    header = list(df.columns)
    header_reqd = ['Gene','Quantity']
    if 'Gene' not in header:
        logger.error('Biologial knowledge file %s must have a column named Gene containing gene symbols.',knowledgeFilepath)
        sys.exit()
    if 'Quantity' not in header:
        logger.error('Biological knowldege file %s must contain a column named Quantity containing the quantity used to rank genes by biological knowledge.',knowledgeFilepath)

    #store data into dict
    for index,row in df.iterrows():
        gene = str(row['Gene'])
        quantity = float(row['Quantity'])
        if isinf(quantity) or isnan(quantity):
            continue
        v2knowledge[gene] = float(quantity)

    logger.info('Read biological knowledge for %d genes from %s',len(v2knowledge),knowledgeFilepath)

    return(v2knowledge)
def readBiologicalKnowledge_pVal(args):
    """
    Reads file called biologicalknowledge.
    File may have extension .xlsx, .csv, or .tsv.
    File must have columns named Gene and Quanity.
    Gene is a column of gene symbols.
    Returns dict of gene to pvalue.
    """
    v2pVal = dict()

    #find biological knowledge filepath
    knowledgeFilepath = args.biologicalknowledge
    if os.path.isfile(knowledgeFilepath):
        #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',knowledgeFilepath)
        pass
    else:
        knowledgeFilepath = args.biologicalknowledge
        if os.path.isfile(knowledgeFilepath): 
            #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .tsv extension.',knowledgeFilepath)
            pass
        else:
            knowledgeFilepath = args.biologicalknowledge
            if os.path.isfile(knowledgeFilepath): 
                #logger.info('Found biological knowledge file %s.',knowledgeFilepath)
                pass
            else:
                logger.error('No biological knowledge file found. Must provide a file named biologicalknowledge with .xlsx, .csv, or .tsv extension in directory %s.', PROJDIR)
                sys.exit()

    #check if file is empty
    if not os.stat(knowledgeFilepath).st_size:
        logger.error('File %s is empty.',knowledgeFilepath)
        sys.exit()

    #read data
    if knowledgeFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(knowledgeFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Biological knowledge file %s must contain only 1 sheet.', knowledgeFilepath)
            sys.exit()
        df = pd.read_excel(knowledgeFilepath)
    elif knowledgeFilepath.split(".")[-1] == 'csv':
        df = pd.read_csv(knowledgeFilepath, sep=',', engine='python')
    else:
        df = pd.read_csv(knowledgeFilepath, sep='\t', engine='python')

    #check for columns
    header = list(df.columns)
    header_reqd = ['Gene','P-value']
    if 'Gene' not in header:
        logger.error('Biologial knowledge file %s must have a column named Gene containing gene symbols.',knowledgeFilepath)
        sys.exit()
    if 'P-value' in header:
        for index,row in df.iterrows():
            gene = str(row['Gene'])
            v2pVal[gene] = float(row['P-value'])
    else:
        logger.info('Biological knowledge file %s does not contain a column named P-Value containing the pvalue of each gene.',knowledgeFilepath)

    #store data into dict


    logger.info('Read p-value for %d genes from %s',len(v2pVal),knowledgeFilepath)

    return(v2pVal)

def readBiologicalKnowledge_logFC(args):
    """
    Reads file called biologicalknowledge.
    File may have extension .xlsx, .csv, or .tsv.
    File must have columns named Gene and Quanity.
    Gene is a column of gene symbols.
    Returns dict of gene to pvalue.
    """
    v2log2FC = dict()

    #find biological knowledge filepath
    knowledgeFilepath = args.biologicalknowledge
    if os.path.isfile(knowledgeFilepath):
        #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',knowledgeFilepath)
        pass
    else:
        knowledgeFilepath = args.biologicalknowledge
        if os.path.isfile(knowledgeFilepath): 
            #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .tsv extension.',knowledgeFilepath)
            pass
        else:
            knowledgeFilepath = args.biologicalknowledge
            if os.path.isfile(knowledgeFilepath): 
                #logger.info('Found biological knowledge file %s.',knowledgeFilepath)
                pass
            else:
                logger.error('No biological knowledge file found. Must provide a file named biologicalknowledge with .xlsx, .csv, or .tsv extension in directory %s.', PROJDIR)
                sys.exit()

    #check if file is empty
    if not os.stat(knowledgeFilepath).st_size:
        logger.error('File %s is empty.',knowledgeFilepath)
        sys.exit()

    #read data
    if knowledgeFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(knowledgeFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Biological knowledge file %s must contain only 1 sheet.', knowledgeFilepath)
            sys.exit()
        df = pd.read_excel(knowledgeFilepath)
    elif knowledgeFilepath.split(".")[-1] == 'csv':
        df = pd.read_csv(knowledgeFilepath, sep=',', engine='python')
    else:
        df = pd.read_csv(knowledgeFilepath, sep='\t', engine='python')

    #check for columns
    header = list(df.columns)
    header_reqd = ['Gene','Log2FC']
    if 'Gene' not in header:
        logger.error('Biologial knowledge file %s must have a column named Gene containing gene symbols.',knowledgeFilepath)
        sys.exit()
    if 'Log2FC' in header:
        for index,row in df.iterrows():
            gene = str(row['Gene'])
            v2log2FC[gene] = float(row['Log2FC'])
    else:
        logger.info('Biological knowledge file %s does not contain a column named Log2FC containing the log2 fold change of each gene.',knowledgeFilepath)

    #store data into dict


    logger.info('Read logFC for %d genes from %s',len(v2log2FC),knowledgeFilepath)

    return(v2log2FC)

def getMouseGene2RepurpDrugs(genes,hgene2repurp):
    """
    Returns a dict mapping mouse gene target to repurposing hub drug.
    """
    v2repurp = defaultdict(dict)

    #get mouse to human conversion
    mouse2human = getMouse2Human()

    for gene in genes:
        if gene in mouse2human:
            gene_h = mouse2human[gene]
        else:
            gene_h = gene
        if gene_h in hgene2repurp:
            v2repurp[gene] = hgene2repurp[gene_h]

    return(v2repurp)

def getRepurpDrugs():
    """
    Returns 2 dict mappings: 
    1. human gene symbol to Broad Repurposing Hub drug.
    2. Broad Repurposing Hub drug to human gene symbol.
    """
    hgene2repurp = defaultdict(dict)
    repurp2hgene = defaultdict(dict)
    drugs = dict()
    filepath = os.path.join(DBDIR, 'repurposing_drugs_20200324.txt')
    df = pd.read_csv(filepath, sep='\t', engine='python',comment='!')

    for index,row in df.iterrows():
        drug = str(row['pert_iname'])
        genes = str(row['target']).split('|')
        for gene in genes:
            if gene != 'nan':
                hgene2repurp[gene][drug] = True
                repurp2hgene[drug][gene] = True
                drugs[drug] = True

    logger.info('Read %d drugs and %d targets from the Broad Drug Repurposing Hub.', len(drugs),len(hgene2repurp))
    return(hgene2repurp,repurp2hgene)

def getMouse2Human():
    """
    Returns a dict mapping of mouse gene symbol to homolog human gene symbol. 
    """

    #Check if mouse to human gene mapping file exists
    filepath = os.path.join(DBDIR, 'mouse_to_human_genes.tsv')
    if not os.path.isfile(filepath):
        #read data
        rawfilepath = os.path.join(DBDIR, 'HOM_MouseHumanSequence.rpt')
        #logger.info('Compiling mouse to human gene mapping from %s...',rawfilepath)
        data = pd.read_csv(rawfilepath, sep='\t')

        num2species = dict()
        m2h = dict()
        mouse2human = dict()

        for index,row in data.iterrows():
            num = row['DB Class Key']
            if num not in num2species.keys():
                num2species[num] = dict()
                num2species[num]['mouse'] = dict()
                num2species[num]['human'] = dict()
            species = row['Common Organism Name']
            symbol = row['Symbol']
            if species == 'mouse, laboratory':
                num2species[num]['mouse'][symbol] = True
            if species == 'human':
                num2species[num]['human'][symbol] = True

        for num in num2species.keys():
            for m in num2species[num]['mouse'].keys():
                if m not in m2h.keys():
                    if num2species[num]['human']:
                        m2h[m] = dict()
                for h in num2species[num]['human'].keys():
                    m2h[m][h] = True
            
        for m in m2h.keys():
            if len(m2h[m]) > 1:
                if m.upper() in m2h[m].keys():
                    mouse2human[m]= m.upper()
                else:
                    choice = str(list(m2h[m].keys())[0])
                    mouse2human[m] = str(choice)
            else:
                mouse2human[m] = str(list(m2h[m].keys())[0])

        hgenes = dict()

        f = open(filepath,'w')
        f.write("mouse\thuman\n")
        for m in mouse2human.keys():
            h = mouse2human[m]
            f.write("%s\t%s\n" % (m,h))
            hgenes[h] = True
        f.close()
        logger.info('Mapped %d mouse genes to %d human genes. Wrote to %s.',len(mouse2human),len(hgenes),filepath)
    else:
        #read the existing file
        data = pd.read_csv(filepath, sep='\t')

        mouse2human = dict()

        for index,row in data.iterrows():
            m = row['mouse']
            h = row['human']
            mouse2human[m] = h

    return(mouse2human)

def checkKeys(dict,key):
    if key in dict.keys():
        return "{:.2e}".format(dict[key])
    else:
        return ""
    
def writeWeights(args,v2weight,v2rank,v2pathType,g2repurp,repurp2hgene,vertices_ranked,v2pVal,v2log2FC,drug_target_max=5):
    """
    Writes NetResponse weights to network directory.
    Adds Repurposing Hub compounds.
    """
    filepath = os.path.join(PROJDIR,'/geneRankings.tsv' )
    fp = open(filepath,'w')
    if len(v2pVal)!=0 and len(v2log2FC)!=0:
        fields = ['gene','rank','type','weight','pvalue','log2FC','repurposinghub_cpds']
        fp.write('\t'.join(fields) + '\n')
        for v in vertices_ranked:
            repurps = dict()
            repurplist = list()
            if v in g2repurp:
                repurps = g2repurp[v]
            for drug in repurps.keys():
                l = len(repurp2hgene[drug])
                if l <= drug_target_max:
                    repurplist.append(drug)
            repurpstr = ','.join(repurplist)

            toks = [v,str(v2rank[v]),v2pathType[v],"{:.2e}".format(v2weight[v]),checkKeys(v2pVal,v),checkKeys(v2log2FC,v),repurpstr]
            fp.write('\t'.join(toks) + '\n')
        fp.close()
        logger.info('Wrote ranked genes to %s.',filepath)
    elif len(v2pVal)!=0 and len(v2log2FC)==0:
        fields = ['gene','rank','type','weight','pvalue','repurposinghub_cpds']
        fp.write('\t'.join(fields) + '\n')
        for v in vertices_ranked:
            repurps = dict()
            repurplist = list()
            if v in g2repurp:
                repurps = g2repurp[v]
            for drug in repurps.keys():
                l = len(repurp2hgene[drug])
                if l <= drug_target_max:
                    repurplist.append(drug)
            repurpstr = ','.join(repurplist)

            toks = [v,str(v2rank[v]),v2pathType[v],"{:.2e}".format(v2weight[v]),checkKeys(v2pVal,v),repurpstr]
            fp.write('\t'.join(toks) + '\n')
        fp.close()
        logger.info('Wrote ranked genes to %s.',filepath)        
    elif len(v2pVal)==0 and len(v2log2FC)!=0:
        fields = ['gene','rank','type','weight','log2FC','repurposinghub_cpds']
        fp.write('\t'.join(fields) + '\n')
        for v in vertices_ranked:
            repurps = dict()
            repurplist = list()
            if v in g2repurp:
                repurps = g2repurp[v]
            for drug in repurps.keys():
                l = len(repurp2hgene[drug])
                if l <= drug_target_max:
                    repurplist.append(drug)
            repurpstr = ','.join(repurplist)

            toks = [v,str(v2rank[v]),v2pathType[v],"{:.2e}".format(v2weight[v]),checkKeys(v2log2FC,v),repurpstr]
            fp.write('\t'.join(toks) + '\n')
        fp.close()
        logger.info('Wrote ranked genes to %s.',filepath)
    else:
        fields = ['gene','rank','type','weight','repurposinghub_cpds']
        fp.write('\t'.join(fields) + '\n')
        for v in vertices_ranked:
            repurps = dict()
            repurplist = list()
            if v in g2repurp:
                repurps = g2repurp[v]
            for drug in repurps.keys():
                l = len(repurp2hgene[drug])
                if l <= drug_target_max:
                    repurplist.append(drug)
            repurpstr = ','.join(repurplist)

            toks = [v,str(v2rank[v]),v2pathType[v],"{:.2e}".format(v2weight[v]),repurpstr]
            fp.write('\t'.join(toks) + '\n')
        fp.close()
        logger.info('Wrote ranked genes to %s.',filepath)                
    return None

def getRankings(v2weight,ties='high',reverse_order=True):
    """
    Returns the following:
    vertices_ranked: list of genes ranked by v2weight. Default order is descending (reverse), and
    v2rank: a dict of gene names to rank. Default order is descending (reverse).
    """
    pairs = sorted([ (v2weight[v], v) for v in v2weight ], reverse = reverse_order)
    vertices_ranked = [ v for (d,v) in pairs ]
    v2rank = dict()
    if ties == 'high':
        if reverse_order:
            rank = 1
            (d0,v0) = pairs[0]
            v2rank[v0] = rank
            previous_d =d0
            counter = 0
            for (d,v) in pairs:
                if isclose(d,previous_d) and v != v0:
                    v2rank[v] = rank
                    counter += 1
                elif d < previous_d:
                    rank = rank + counter + 1
                    counter = 0
                    v2rank[v] = rank
                    previous_d = d
                elif d > previous_d:
                    logger.error('Pairs are not ranked')
                    sys.exit()
        else:
            rank = 1
            (d0,v0) = pairs[0]
            v2rank[v0] = rank
            previous_d =d0
            counter = 0
            for (d,v) in pairs:
                if isclose(d,previous_d) and v != v0:
                    v2rank[v] = rank
                    counter += 1
                elif d > previous_d:
                    rank = rank + counter + 1
                    counter = 0
                    v2rank[v] = rank
                    previous_d = d
                elif d < previous_d:
                    logger.error('Pairs are not ranked')
                    sys.exit()

    return(vertices_ranked,v2rank)


def getWeights(args,g_list,v2i,driver,responses,nt,dt,endpoints=False):
    """
    Returns the weight for each gene in the network in a dict.
    """
    logger.info('Calculating NetResponse weights for network genes...')
    v2w = dict()
    vertices = sorted(v2i.keys())

    for v in vertices:
        v2w[v] = 0.0
    d = driver
    a = v2i[d]
    for r in responses:
        b = v2i[r]
        for v in vertices:
            if not endpoints:
                if v == d or v == r:
                    continue
            x = v2i[v]
            y = [(g_list[nt-i][b,x])*(g_list[i][x,a]) for i in range(nt+1)]
            z = integrate.romb(y,dx=dt)
            v2w[v] = v2w[v] + z
    norm = 0
    for v in vertices:
        norm += v2w[v]
    if norm == 0.:
        norm = 1.
    for v in vertices:
        v2w[v] = v2w[v] / norm

    return(v2w)

def getGlist(lap_mat,dt,nt):
    """
    Returns the G(t) as a list of matrices.
    """
    logger.info('Calculating the Green\'s function...')
    ret = [None] * (nt + 1)
    (nr,nc) = lap_mat.shape
    assert(nr == nc), 'Bad laplacian shape: %s' % str(lap_mat.shape)
    my_eye = np.eye(nr)
    epsilon = - dt * lap_mat
    g1 = linalg.expm(epsilon)
    curval = my_eye
    ret[0] = curval
    for i in range(nt):
        curval = np.matmul(curval, g1)
        ret[i + 1] = curval
    return(ret)

def getTimeStep(driver,v2outdeg,r=0.5,nt=16):
    """
    Returns the time step size.
    """
    # r is the relaxation parameter (the reduction in driver activity)
    # nt is the number of time steps to reach r (must be a non-negative power of 2)
    
    #calculate dt (time step size)
    den = 1 - r
    dt = -log(den,e) / (v2outdeg[driver] * nt)

    #logger.info('dt calculated to be ~%f. Driver density will reach %g in %d steps. Driver outdeg is %d.', round(dt,5),den,nt,v2outdeg[driver])
    return(dt)

def getLaplacian(edge_dict,v2i):
    """
    Returns the adjacency matrix and Laplacian.
    """
    nvert = len(v2i)
    adj_mat = np.zeros( (nvert,nvert) )
    indeg_arr = np.zeros(nvert)
    outdeg_arr = np.zeros(nvert)
    lap_mat = np.zeros( (nvert,nvert) )
    
    
    for (src,dest) in edge_dict.keys():
        (a,b) = (v2i[src],v2i[dest])
        adj_mat[b,a] += 1.0
        indeg_arr[b] += 1.0
        outdeg_arr[a] += 1.0
    for a in range(nvert):
        lap_mat[a,a] = outdeg_arr[a]
    lap_mat = lap_mat - adj_mat
    return(adj_mat, lap_mat)


def nameToType(genelist,driver,intermediates,responses):
    """
    Returns a dict of gene names to type.
    D: driver.
    I: intermediate.
    R: response.
    """
    v2type = dict()
    for v in genelist:
        if v == driver:
            v2type[v] = 'D'
        elif v in intermediates:
            v2type[v] = 'I'
        elif v in responses:
            v2type[v] = 'R'
        else:
            logger.warning('Gene %s is not in the driver, intermediate, or response lists.',v)
    return(v2type)

def nameToDegree(genelist,edge_ba,edge_ab):
    """
    Returns indegree and outdegree dicts.
    Gene name is the key.
    """
    v2indeg = dict()
    v2outdeg = dict()
    for v in genelist:
        v2indeg[v] = 0
        if v in edge_ba:
            v2indeg[v] = len(edge_ba[v].keys())
        v2outdeg[v] = 0
        if v in edge_ab:
            v2outdeg[v] = len(edge_ab[v].keys())
    return(v2indeg,v2outdeg)


def nameToIndex(genelist):
    """
    Reads a list of gene names.
    Returns dicts of names to indices and vice versa
    """
    n2i = dict()
    i2n = dict()
    names = sorted(genelist)
    for (n,i) in zip( names, list(range(len(names))) ):
        n2i[n] = i
        i2n[i] = n
    return(n2i,i2n)

def edgeTuplesToDicts(edge_dict):
    """
    Reads dict of edges.
    Returns two dict-of-dicts.
    a is source, be is destination.
    ab uses sources as first key.
    ba uses dests as first key.
    """
    edge_ab = dict()
    edge_ba = dict()
    for (a,b) in edge_dict.keys():
        if a not in edge_ab:
            edge_ab[a] = dict()
        edge_ab[a][b] = True
        if b not in edge_ba:
            edge_ba[b] = dict()
        edge_ba[b][a] = True

    return(edge_ab, edge_ba)

def readEdges(args,verbose=False):
    """
    Reads interactions.
    Returns a dictionary of directed edges.
    Undirected edges are recorded as two directed edges in the dictionary.
    """
    if args.s == 'human':
        filepath = os.path.join(DBDIR, 'human_interactions.txt')
    elif args.s == 'mouse':
        filepath = os.path.join(DBDIR, 'mouse_interactions.txt')

    #check filepath exists
    if not os.path.isfile(filepath):
        logger.error('Missing database file %s. Run NetResponse analysis.',filepath)
        sys.exit()

    ret = dict()
    ndirected = 0
    nundirected = 0
    fp = open(filepath, 'r')
    for line in fp:
        my_list = [ ]
        toks = line.strip().split()
        if (len(toks) < 2) or (len(toks) >3):
            if verbose:
                logger.warning('Bad line, wrong token count: %s', line)
            continue
        (a, b) = (toks[0], toks[1])
        if len(toks)==2:
            my_list.append( (a,b) )
            ndirected += 1
        elif len(toks)==3 and toks[2]=='TF':
            my_list.append( (a,b) )
            ndirected += 1
        elif len(toks)==3 and toks[2]=='PPI':
            my_list.append( (a,b) )
            my_list.append( (b,a) )
            nundirected += 1
        else:
            if verbose:
                logger.warning('Bad line, unknown edge type: %s', line)
            continue
        for key in my_list:
            ret[key] = True

    nedge = len(ret)
    if verbose:
        logger.info('Network consists of %d edges: %d directed edges and %d undirected edges from %s.', nedge, ndirected, nundirected, filepath)
    return(ret)

def checkDriver(args,wholegenelist,verbose=False):
    """
    Reads the driver from the commandline.
    Checks if driver is in the network.
    Outputs driver as a string.
    """
    #read driver
    ret = args.driver

    #check if driver is in the network
    if ret in wholegenelist:
        if verbose:
            logger.info('Found driver %s in the network.' % ret)
    else:
        logger.error('Driver %s was not found in the network. Pick another driver.' % ret)
        sys.exit()

    return(ret)

def readDriver(args):
    """
    Reads the driver from the commandline.
    Outputs driver as a string.
    """
    #read driver
    ret = args.driver

    return(ret)

def readResponses(args,driver):
    """
    Reads a file provided by the user containing the gene symbols of the response genes.
    Removes driver from the response genes.
    Returns a list of response genes.
    """

    responseFilepath = args.responses
    if os.path.isfile(responseFilepath):
        #logger.info('Found response file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',responseFilepath)
        pass
    else:
        responseFilepath = args.responses
        if os.path.isfile(responseFilepath):
            #logger.info('Found response file %s. Delete this file if you want to use a file with a .tsv extension.',responseFilepath)
            pass
        else:
            responseFilepath = args.responses
            if os.path.isfile(responseFilepath):
                #logger.info('Found response file %s.',responseFilepath)
                pass
            else:
                logger.error('No response file found. Must provide a file named responses with .xlsx, .csv, or .tsv extension in directory %s.',PROJDIR)
                sys.exit()

    #check if file is empty
    if not os.stat(responseFilepath).st_size:
        logger.error('The response file %s is empty. It must contain at least one gene symbol.',responseFilepath)
        sys.exit()

    #read responses
    if responseFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(responseFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Response file %s must contain only 1 sheet.',responseFilepath)
            sys.exit()
        responses = pd.read_excel(responseFilepath,header=None)
    elif responseFilepath.split(".")[-1] == 'csv':
        responses = pd.read_csv(responseFilepath,sep=',',header=None)
    else:
        responses = pd.read_csv(responseFilepath,sep='\t',header=None)

    responses = list(set(responses[0].tolist()))

    #remove driver from responses
    if driver in responses:
        responses.remove(driver)

    return(responses)



def getNetworkResponses(args,driver,wholegenelist,verbose=False):
    """
    Reads a file provided by the user containing the gene symbols of the response genes.
    Removes driver from the response genes.
    Checks if each response gene is in the network.
    Outputs list of network response genes called networkResponse.
    Returns a list of network responses.
    """

    #check if network response file exists
    networkResponsesFilepath = os.path.join(PROJDIR,'/networkResponses.txt' ) 
    if os.path.isfile(networkResponsesFilepath):
        if verbose:
            logger.info('Found network responses file %s. Delete this file if you want the network respones to be compiled again.',networkResponsesFilepath)
        
        #check if file is empty
        if not os.stat(networkResponsesFilepath).st_size:
            logger.error('The network responses file is empty. It must contain at least one gene symbol. Please delete this file and run NetResponse again.')
            sys.exit()

        responses = pd.read_csv(networkResponsesFilepath,sep='\t',header=None)
        responses = list(set(responses[0].tolist()))

        #remove driver from network responses
        if driver in responses:
            responses.remove(driver)
            if verbose:
                logger.warning('Found driver %s in the network response file. It was removed.' % (driver))
    
        #check if responses are in network and output list
        ret = list()
        for response in responses:
            if response in wholegenelist:
                ret.append(response)
            else:
                logger.error('Not all network response genes in %s are in the network. Please delete this file and run NetResponse analysis, again.', networkResponsesFilepath)
                sys.exit()

    else:
        #check if user response file exists
        responseFilepath = args.responses
        if os.path.isfile(responseFilepath):
            #logger.info('Found response file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',responseFilepath)
            pass
        else:
            responseFilepath = args.responses
            if os.path.isfile(responseFilepath):
                #logger.info('Found response file %s. Delete this file if you want to use a file with a .tsv extension.',responseFilepath)
                pass
            else:
                responseFilepath = args.responses
                if os.path.isfile(responseFilepath):
                    #logger.info('Found response file %s.',responseFilepath)
                    pass
                else:
                    logger.error('No response file found. Must provide a file named responses with .xlsx, .csv, or .tsv extension in directory %s.',PROJDIR)
                    sys.exit()

        #check if file is empty
        if not os.stat(responseFilepath).st_size:
            logger.error('Response file %s is empty. Must contain at least one gene symbol.' % responseFilepath)
            sys.exit()
        
        #read responses and convert to list
        if responseFilepath.split(".")[-1] == 'xlsx':
            xl = pd.ExcelFile(responseFilepath)
            res = len(xl.sheet_names)
            if res != 1:
                logger.error('Response file %s must contain only 1 sheet.',responseFilepath)
                sys.exit()
            responses = pd.read_excel(responseFilepath,header=None)
        elif responseFilepath.split(".")[-1] == 'csv':
            responses = pd.read_csv(responseFilepath,sep=',',header=None)
        else:
            responses = pd.read_csv(responseFilepath,sep='\t',header=None)
        
        responses = list(set(responses[0].tolist()))
        if verbose:
            logger.info('Read %s response genes from %s.' % (str(len(responses)),responseFilepath))
    
        #remove driver from responses
        if driver in responses:
            responses.remove(driver)
            if verbose:
                logger.info('Found driver %s in the response list. It was removed. %s response genes remain.' % (driver,str(len(responses))))
        else:
            if verbose:
                logger.info('Did not find driver %s in the response list.' % driver)
    
        #check if responses are in network and output list
        ret = list()
        outfilepath = os.path.join(PROJDIR,'/networkResponses.txt')
        fp = open(outfilepath, 'w')
        for response in responses:
            if response in wholegenelist:
                ret.append(response)
                fp.write(response + '\n')
        fp.close()
        if verbose:
            logger.info('Found %s of %s response genes in the network. Wrote to %s.' % (str(len(ret)),str(len(responses)),outfilepath))
    return(ret)

def readInteractions(filepath):
    """
    Reads a NetResponse Interaction file.
    Returns a pandas dataframe
    """
    ret = pd.read_csv(filepath, sep='\t',header=None)
    return(ret)

def readGenes(filepath):
    """
    Read a NetResponse Gene file.
    Returns a list
    """
    ret = pd.read_csv(filepath,sep='\t',header=None)
    ret = list(set(ret[0].tolist()))
    return(ret)

def readIIDMouse():
    """
    Reads mouse ppi file database/mouse_annotated_PPIs.txt.
    Converts ppi from uniprot id to mgi symbol.
    Returns a list of unique ppi interactions as tuple (a,b).
    a and b are the genes of the two interacting proteins.
    Removes unwanted ppi between a protein and itself.
    """
    filename = os.path.join(DBDIR, 'mouse_annotated_PPIs.txt')

    #get mapping from mouse uniprot id to mgi symbol
    uniprot2mgisymbol = get_uniprot2mgisymbol()

    #read the ppi file
    #read_iid_ppi
    logger.info('Reading mouse PPI from %s. Removing PPI between a protein and itself...', filename)
    reftable = pd.read_csv(filename,sep='\t',low_memory=False)
    ret = dict()
    for index,row in reftable.iterrows():
        (a,b) = ( row['uniprot1'], row['uniprot2'] )
        if a == b:
            continue
        if a in uniprot2mgisymbol.keys() and b in uniprot2mgisymbol.keys():
            for genea in uniprot2mgisymbol[a]:
                for geneb in uniprot2mgisymbol[b]:
                    if genea == geneb:
                        continue
                    if (genea,geneb) not in ret and (geneb,genea) not in ret:
                        ret[(genea,geneb)] = True
    logger.info('Mapped %d mouse PPI interactions from uniprot id to mgi symbol.', len(ret))
    return(list(ret.keys()))

def get_uniprot2mgisymbol():
    """
    Returns a dict mapping from mouse uniprot id to mgi symbol.
    """
    #get uniprot to mgi num mapping
    uniprot2mginum = get_uniprot2mginum()
    #get mgi num to mgi symbol mapping
    mginum2symbol = get_mginum2symbol()
    ret = dict()
    counter = 0
    for prot in uniprot2mginum.keys():
        ret[prot] = dict()
        for num in uniprot2mginum[prot]:
            if num in mginum2symbol:
                ret[prot][mginum2symbol[num]] = True
                counter += 1
    logger.info('Mapped %d mouse uniprot ids to %d mgi symbols.', len(ret), counter)
    
    return(ret)

def get_mginum2symbol():
    """
    Reads database/MRK_List2.rpt'
    Returns a dict mapping from mouse mgi num to symbol.
    """
    mginum2symbol_file = os.path.join(DBDIR, 'MRK_List2.rpt')
    logger.info('Reading mgi marker accession number to mgi marker symbol from %s...', mginum2symbol_file)
    reftable = pd.read_csv(mginum2symbol_file,sep='\t')
    ret = dict()
    for index,row in reftable.iterrows():
        num = row['MGI Accession ID']
        symbol = row['Marker Symbol']
        ret[num] = symbol
    logger.info('Mapped %d mgi marker accession numbers to symbols.',len(ret))
    return(ret)

def get_uniprot2mginum():
    """
    Reads database/MOUSE_10090_idmapping.dat.
    Returns a dict mapping from mouse uniprot id to MGI number.
    """
    uniprot2mginum_file = os.path.join(DBDIR, 'MOUSE_10090_idmapping.dat')
    logger.info('Reading uniprot to mgi marker accession number from %s...', uniprot2mginum_file)
    ret = dict()
    ret2 = dict()
    reftable = pd.read_csv(uniprot2mginum_file, sep='\t')
    for index,row in reftable.iterrows():
        name = row[0]
        if name not in ret:
            ret[name] = dict()
        id_type = row[1]
        if id_type == 'MGI':
            mgi = row[2]
            ret[name][mgi] = True
    for prot in ret.keys():
        if ret[prot]:
            ret2[prot] = ret[prot]
    logger.info('Mapped %d of %d uniprot ids to a MGI number.', len(ret2), len(ret))
    return(ret2)

def readIIDHuman():
    """
    Reads ppi from Human IID file in database dir.
    Returns list of unique ppi interactions as tuple (src,des).
    Removes unwanted PPI between a protein and itself.
    """
    filename = os.path.join(DBDIR, 'human_annotated_PPIs.txt')
    logger.info('Reading human PPI from %s. Removing PPI between a protein and itself...', filename)
    ref = pd.read_csv(filename, sep='\t',low_memory=False)

    ret = dict()

    for index,row in ref.iterrows():
        (a,b) = (row['symbol1'],row['symbol2'])
        if len(a.strip().split()) != 1:
            continue
        if len(b.strip().split()) != 1:
            continue
        if a == b:
            continue
        if (a,b) not in ret and (b,a) not in ret:
            ret[(a,b)] = True

    logger.info('Read %d unique human PPI.', len(ret))
    return(list(ret.keys()))

def readTF(species):
    """
    Reads regulator interactions from files with format: GENE1  GENE2   Mode    PMID.
    Files stored in database dir. 
    File should not have a header.
    """
    if species == 'human': 
        filename = os.path.join(DBDIR, 'trrust_rawdata.human.tsv')
    elif species == 'mouse':
        filename = os.path.join(DBDIR, 'trrust_rawdata.mouse.tsv')
    else:
        logger.error('%s is not a supported species.',species)
        sys.exit()

    ret = dict()
    fp = open(filename,'r')
    for line in fp:
        toks = line.strip().split('\t')
        if len(toks)!=4:
            logger.warning('Bad line: wrong token count %s',line)
            continue
        (genea,geneb) = (toks[0],toks[1])
        ret[(genea,geneb)] = True
    logger.info('Read %d %s TF regulatory interactions.',len(ret),species)
    return(list(ret.keys()))

def main():
    args = CMDParser()
    args.func(args)
    return None

if __name__=='__main__':
    main()
