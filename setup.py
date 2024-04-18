from setuptools import setup, find_packages, Command
from setuptools.command.install import install
import subprocess
import os

class CustomInstallCommand(install):
    def run(self):
        # Get the directory where setup.py is located
        base_dir = os.path.abspath(os.path.dirname(__file__))
        # Path to build.sh relative to setup.py
        build_script = os.path.join(base_dir, 'netpert', 'bin','setup.sh')
        # Execute build.sh
        subprocess.call(['bash', build_script])
        install.run(self)
setup(
    name="netpert",
    author="Pranhav Sundararajan",
    author_email="pranhav16@gmail.com",
    description="Official NetPert Tool",
    url="https://github.com/joelbaderlab/netresponse_dev",
    cmdclass={
        'install': CustomInstallCommand,
    },
    packages=find_packages(),
    package_data={
        'netpert': ['databases/*.rpt', 'databases/*.tsv','databases/*.dat', 'databases/*.xlsx', 'databases/*.txt' ,'databases/*','projects/*','bin/*'], # Include all .txt files in the 'data' directory
    },
    entry_points={
        'console_scripts': [
            'netpert = netpert.bin.netResponse:main',
        ],
    },
    install_requires=['pandas','numpy','scipy','openpyxl'],
    
)
