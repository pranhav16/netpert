"""
netresponse snakefile (c) Joel Bader Lab 2023
"""

import os

"""
Everything is underneath the root directory
"""
toplevel = config['toplevel']
project = config['project']

print(f'toplevel={toplevel}\nproject={project}')

# bin is a command in python to convert to binary, so use bindir
bindir = os.path.join(toplevel, 'bin')
dbdir = os.path.join(toplevel, 'databases')
projects = os.path.join(toplevel, 'projects')
projdir = os.path.join(projects, project)
print(f'bindir={bindir}\ndbdir={dbdir}\nprojects={projects}\nprojdir={projdir}')

netresponse = os.path.join(bindir, 'netResponse.py')
# can activate a conda environment calling snakemake with --use-conda
# and in a rule supplying either a yaml or an existing named environment
netresponseyaml = os.path.join(bindir, 'netResponseEnv.yaml')
netresponseenv = 'netResponseEnv'

# databases
human_ppi = os.path.join(dbdir, 'human_annotated_PPIs.txt')
mouse_ppi = os.path.join(dbdir, 'mouse_annotated_PPIs.txt')
human_tf = os.path.join(dbdir, 'trrust_rawdata.human.tsv')
mouse_tf = os.path.join(dbdir, 'trrust_rawdata.mouse.tsv')
repurposinghub = os.path.join(dbdir, 'repurposing_drugs_20200324.txt')

# project data
driver = 'Twist1'
georgess = os.path.join(projdir, 'Georgess_Prkd1_2019_Data_Final.xlsx')
jcb = os.path.join(projdir, 'supp_jcb.201306088_JCB_201306088_TableS1.xlsx')

# output files
generankings = os.path.join(projdir, 'geneRankings.tsv')
table1 = os.path.join(projdir, 'Table1.tsv')
table2 = os.path.join(projdir, 'Table2.tsv')

"""
Underneath the named project are reads and results
"""

# dummy file to force execution
dummy = os.path.join(toplevel, 'dummy.txt')

rule all:
    input:
        human_ppi,
        mouse_ppi,
        human_tf,
        mouse_tf,
        repurposinghub,
        generankings
#        table1,
#        table2

rule build_dbs:
    input:
    output:
        human_ppi,
        mouse_ppi,
        human_tf,
        mouse_tf,
        repurposinghub
    shell:
        """
        # the -p flag creates directories along the full path as required
        mkdir -p {dbdir}
        # Download and unzip large files
        cd {dbdir}
        curl -O http://iid.ophid.utoronto.ca/static/download/mouse_annotated_PPIs.txt.gz ; gunzip mouse_annotated_PPIs.txt.gz
        curl -O http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz ; gunzip human_annotated_PPIs.txt.gz
        curl -O https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
        curl -O https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv
        curl -O https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
        curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz ; gunzip MOUSE_10090_idmapping.dat.gz
        curl -O https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt
        curl -O https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt
        cd {toplevel}
        """

rule build_twist1:
    input:
    output:
        georgess = georgess,
        jcb = jcb
    shell:
        """
        cd {projdir}
        curl -OL https://github.com/EwaldLab/2019_Prkd1/raw/master/Georgess_Prkd1_2019_Data_Final.xlsx
        curl -O https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3941052/bin/supp_jcb.201306088_JCB_201306088_TableS1.xlsx
        """
        
rule twist1_analysis:
    conda: netresponseenv
    threads: workflow.cores
    input:
        mouse_ppi,
        mouse_tf,
        jcb
    output: generankings
    shell:
        """
        /usr/bin/env python {netresponse} analysis -s mouse {projdir} {driver}
        """

rule twist1_comparison:
    conda: netresponseenv
    input: generankings
    output:
        table1,
        table2
    shell:
        """
        /usr/bin/env python {netresponse} comparison -s mouse {projdir} {driver}
        """
