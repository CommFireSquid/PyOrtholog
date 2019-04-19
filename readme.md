# PyOrtholog
Python script for collection of representative genes for a given gene name (i.e. ACAN)

These scripts should be platform independent, but it is primarily tested on an Ubuntu 16.04 LTS +system

The data directory is not included as it is too large, but scripts will be setup to create it soon as I manually instantiated it. The only important part is the database as the .txt files are residuals of the init script runs. To initialize the database only use the dbLocal_Load.py script as it is a combination of all other init scripts. Read the usage function within the script for information about using the arguments.

To run the geneOrthologCollect.py script, include the arguments for gene name, email, and API key (you will need to register with NCBI to get one) as the script will not work without an API key due to the query volume.

Scripts included under Misc were utilized as transitional steps between collection and alignment, or to take the alignment and create a data table for analysis in R. The descriptions of these scripts and how to use them is in the final thesis.

Finally, there have been several changes since the creation of the ACAN files, so do not trust the results you see there. A quick glance will show that the method used to generate these files does not include complements in cases where the coding strand is not the strand provided by the genome. This error was corrected for with the implementation of Biopython additions. Essentially, the files created by the scripts in their current state will not match the sample until I get a chance to update them.

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
