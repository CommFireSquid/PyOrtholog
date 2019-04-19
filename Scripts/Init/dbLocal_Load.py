# Script to load all dbLocal.sqlite3 tables into the database as an initial load
#
#	Bob Anderson
#

import os
import sys
import gzip
import ftplib
import getopt
import shutil
import sqlite3

ftpEmail = ''
datFileDirectory = './Data/'
datTempDirectory = datFileDirectory + 'temp/'
database = datFileDirectory + 'dbLocal.sqlite3'
commitSize = 100000

# FTP files that need to be downloaded from NCBI
ftpTargets = [('/gene/DATA/','gene2refseq.gz'),('/pub/taxonomy/','taxdump.tar.gz'),('/gene/DATA/','gene_orthologs.gz')]
desiredFiles = ['gene2refseq','gene_orthologs','./taxdump/names.dmp']

def main():
	global ftpEmail

	try:
		options, arguments = getopt.getopt(sys.argv[1:],"e:")
		for (opt, arg) in options:
			if (opt == "-e"):
				ftpEmail = arg
	except getopt.GetoptError as OPT:
		usage(OPT)
	if ftpEmail == '':
		usage('-e Not Entered')

	'''# Create the temp data directory
	try:
		os.mkdir(datTempDirectory)
	except FileExistsError:
		shutil.rmtree(datTempDirectory)
		os.mkdir(datTempDirectory)

	for target in ftpTargets:
		ftpHarvest(target)
		print("Completed FTP of target: {}".format(target))'''

	if verifyFTPFiles(desiredFiles):
		writeToDB(desiredFiles)
		# Delete temp files
		shutil.rmtree(datTempDirectory)
	else:
		print("Failed to Validate Files")

def ftpHarvest(ftpFile):
	path, target = ftpFile

	fileOut = open(datTempDirectory + target, "wb")

	ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov','anonymous',ftpEmail)
	ftp.cwd(path)

	ftp.retrbinary("RETR " + target, fileOut.write, 8*1024 )

	ftp.close()
	fileOut.close()


	unzippedTarget = target[:-3]

	with gzip.open(datTempDirectory + target) as fileIn:
		with open(datTempDirectory + unzippedTarget, 'wb') as fileOut:
			shutil.copyfileobj(fileIn, fileOut)

	if unzippedTarget[-4:] == '.tar':
		untarredTarget = unzippedTarget[:-4]
		shutil.unpack_archive(datTempDirectory + unzippedTarget, datTempDirectory + untarredTarget)

def verifyFTPFiles(files):
	for file in files:
		if not os.path.isfile(datTempDirectory + file):
			return False
	return True

def writeToDB(files):
	print('Loading {} Files into Database: {}'.format(len(files),'RefSeq, Orthologs, and Taxonomy'))
	# Write to database for all files in file list:
	#	Load Order
	#		1. refseq
	#		2. ncbi_orthologs
	#		3. taxTable
	connection = sqlite3.connect(database)
	cursor = connection.cursor()
	refSeqID = 0	# Refseq
	orthologID = 1	# Orthologs
	taxonomyID = 2	# Taxonomy

	#---------------------------
	# Prepare tables
	# Create (IF NOT EXISTS) all table and truncate existing values
	sqls = [
		'CREATE TABLE IF NOT EXISTS refseq (taxID INTEGER, geneID INTEGER, status VARCHAR(15), '+\
		'RNA_nucleotide_Accession_Version VARCHAR(20), RNA_Nucleotide_GI INTEGER, protein_accession_version VARCHAR(20), '+\
		'protein_gi INTEGER, genomic_nucleotide_accession_version VARCHAR(20), genomic_nucleotide_gi INTEGER, '+\
		'start_position_on_the_genomic_accession INTEGER, end_position_on_the_genomic_accession INTEGER, '+\
		'orientation VARCHAR(1), assembly VARCHAR(20), mature_peptide_accession_version VARCHAR(20), '+\
		'mature_peptide_gi INTEGER, symbol VARCHAR(20))',
		'CREATE TABLE IF NOT EXISTS ncbi_ortholog (orthID integer PRIMARY KEY, taxID integer, geneID integer, '+\
		'relationship varchar(15), othTaxID integer, othGeneID integer)',
		'CREATE TABLE IF NOT EXISTS taxTable (taxID INTEGER PRIMARY KEY, scientificName VARCHAR(50) NOT NULL, commonName VARCHAR(50))',
		'DELETE FROM refseq',
		'DELETE FROM ncbi_ortholog',
		'DELETE FROM taxTable'
		]
	for statement in sqls:
		cursor.execute(statement)
	connection.commit()

	#---------------------------
	# Load SQLs
	loadSQLs = [
		'INSERT INTO refseq (taxID, geneID, status, RNA_nucleotide_accession_version, '+\
		'RNA_nucleotide_accession_version, protein_accession_version, protein_gi, '+\
		'genomic_nucleotide_accession_version, genomic_nucleotide_gi, start_position_on_the_genomic_accession, '+\
		'end_position_on_the_genomic_accession, orientation, assembly, mature_peptide_accession_version, '+\
		'mature_peptide_gi, symbol)',
		'INSERT INTO ncbi_ortholog (orthID, taxID, geneID, relationship, othTaxID, othGeneID)',
		'INSERT INTO taxTable (taxID, scientificName, commonName)'
	]

	# RefSeq
	commitCTR = 0
	with open(datTempDirectory + files[refSeqID], 'r') as fileIn:
		line = fileIn.readline() # First line will be the table headers which I already have, so it will be ignored
		line = fileIn.readline()

		while line:
			valueString = " VALUES ("
			for x in (line.rstrip()).split('\t'):
				if x == '-':
					valueString = valueString + "NULL, "
				elif isInteger(x):
					valueString = valueString + "{}, ".format(x)
				else:
					x = x.replace("'", "''")
					valueString = valueString + "\'{}\', ".format(x)
			valueString = valueString[:-2] + ")"
			cursor.execute(loadSQLs[refSeqID] + valueString)

			commitCTR += 1
			if commitCTR % commitSize == 0:
				connection.commit()

			line = fileIn.readline()
	connection.commit()
	print("Completed reference sequence load with {} records".format(commitCTR))


	# Orthologs
	# Find the next available PK value for the table and default to 0 if no value is found
	maxPKval = 0
	commitCTR = 0
	with open(datTempDirectory + files[orthologID],'r') as fileIn:
		line = fileIn.readline()	# First line will be a header line that contains the names of the columns, not useful because we already know what the data will contain.
		line = fileIn.readline()

		while line:
			tax, gene, desc, otax, ogene = line.split('\t')
			maxPKval += 1
			cursor.execute(loadSQLs[orthologID] + ' VALUES ({},{},{},\'{}\',{},{})'.format(maxPKval,tax,gene,desc,otax,ogene))
			
			commitCTR += 1
			if commitCTR % commitSize == 0:
				connection.commit()

			line = fileIn.readline()
		connection.commit()
	print("Completed ortholog table load with: {} records".format(commitCTR))

	# Taxonomy
	with open(datTempDirectory + files[taxonomyID], 'r') as fileIn:
		valueStrings = []
		currentTID = -1
		scientificName = ""
		commonName = ""
		line = fileIn.readline()
		# While loop to iterate through the document. Suffers from the fact that it only writes when a new taxID is encountered. This requires an extra write command at the end to handle the last case
		while line:
			values = line.split('|') # Split the fields based on the field seperator
			ctr = 0
			for val in values:	# Clean up the data (i.e. remove all the extra whitespace)
				if ctr == 0:
					try:
						values[ctr] = int(val.rstrip().lstrip())
					except TypeError:
						print('Error handling TaxID: {}\n'.format(val))
						exit(0)
				else:
					values[ctr] = val.rstrip().lstrip().replace("'", "''") # Temporary solution to cases where there are single quotes in the species names
				ctr += 1

			if currentTID > 0 and values[0] != currentTID: # Escape Case -> The case where we have grabbed all the information on a given TaxID and it is ready to be written
				if scientificName != "":
					valString = " VALUES ({}, \'{}\', ".format(currentTID, scientificName)
					if commonName == "":
						valString = valString + "NULL )"
					else:
						valString = valString + "\'{}\' )".format(commonName)
					valueStrings.append(valString)

					
				else:	# If there is no scientific then we will not write anythinig
					logFile.write("Bad TaxID: {} --> No Scientific Name Found\n".format(currentTID))

				if len(valueStrings) == commitSize:
					writeTaxValues(valueStrings, cursor, connection, loadSQLs, taxonomyID)
					valueStrings = []

				# Reassign values so we can get the next TID
				currentTID = -1
				scientificName = ""
				commonName = ""

			# Actual data processing, collect the TID if needed ... so on
			if currentTID < 0: # This is the case when we have moved on to a new TaxID and need to record
				currentTID = values[0]
			if values[3] == "scientific name":	# Case where the scientific name is found
				scientificName = values[1]
			elif values[3] == "common name":	# Case where there is a common name that can be grabbed
				commonName = values[1]

			# Next line reads
			line = fileIn.readline()

	# --------------- END While --------------- #
	# Now handle the potential last case (because of the issue described before while loop)
	if currentTID > 0 and scientificName != "":
		valString = " VALUES ({}, \'{}\', ".format(currentTID, scientificName)
		if commonName == "":
			valString = valString + "NULL )"
		else:
			valString = valString + "\'{}\' )".format(commonName)
		valueStrings.append(valString)

	else:	# If there is no scientific then we will not write anythinig
		print("Bad TaxID: {} --> No Scientific Name Found\n".format(currentTID))

	if len(valueStrings) > 0:
		writeTaxValues(valueStrings, cursor, connection, loadSQLs, taxonomyID)
		valueStrings = []

	print("Completed taxonomy load with {} records".format(currentTID))

def writeTaxValues(valueStrings, databaseC, databaseConn, loadSQLs, taxonomyID):
	# Writes values to the taxTable utilizing the partial SQLs made in the method to insert these values. Commits when the list is exausted
	#		Written on: 11 April 2018
	#		Author: Bob Anderson
	#	No Required methods (Required by parseTaxValues() method)
	# ARGUMENTS;
	#	valueStrings	-	List of values for insert strings in the format: " VALUES ({}, \'{}\', \'{}\')" No quotes on fina value if NULL
	#	databaseC 		-	Cursor from database described in parent method --> described in main()
	#	databaseConn	-	Connection object for the database (needed for commit)
	#	logFile			-	Handle for the log file associated with the run of the script
	for valString in valueStrings:
		databaseC.execute(loadSQLs[taxonomyID] + valString)
	databaseConn.commit()

# Utility Functions
def isInteger(val):
	try:
		val += 1
		return True
	except TypeError:
		return False

# Line Argument Error Message Function
def usage(errMessage):
	# Prints usage information to the log file if there is an issue retrieving the command arguments
	#		Writes to the default error as geneName is needed for regular error file
	#	errMessage	-	error message from getopt (could also be any string passed)
	# 21 May 2018
	print("Below is the error resulting from dbLocal_Load.py")
	print("{}\n".format(errMessage))
	print("Usage:\tdbLocal_Load.py -e email@email.com")
	sys.exit(1)

if __name__ == '__main__':
	main()