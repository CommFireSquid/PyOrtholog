# PyOrtholog
Python script for collection of representative genes for a given gene name (i.e. ACAN)

These scripts should be platform independant, but it is primarily tested on an Ubuntu 16.04 LTS +system

The data directory is not included as it is too large, but scripts will be setup to set this up soon. The only important part is the database as the .txt files are residuals of the init script runs.

################################
File Structure
+ PyOrtholog
  + Data
    > dbLocal.sqlite3
    > gene2accession.txt
    > gene_orthologs.txt
    > taxNames.txt
  + FASTA
    + GeneName...
  + Logs
    + OrthoCollectLogs
      > LogFilePerGeneName...
    > defaultError.log
  + Scripts
    + Process
      + __pycache__
        > OrthologObjects.cpython-35.pyc
      > __init__.py
      > geneOrthologCollect.py
      > OrthologObjects.py
    + Init
      > geneAccessionLoad.py
      > orthologTableLoad.py
      > taxFTP_and_Load.py
