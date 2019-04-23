import os
import re
import sys
import gzip
import ftplib
import getopt
import shutil
import sqlite3

# Output Documentation
availabilityOut = [['Pathway','Gene','TaxCount']]

# Directory values
alignDir = 'Alignments/'
alignIn = alignDir + 'Uncleaned Alignments/' # Where to get the files that need their species names cleaned
alignOut = alignDir + 'Clean Alignments/' # Where to put the files with cleaned species names
dataDir = 'Data/'
docDir = 'Documentation/Alignment Metrics'
alignAvail = 'AvailableAlignments.csv'
datTempDirectory = 'temp/'

# Database values
dbname = dataDir + 'dbLocal.sqlite3'
commitSize = 100000
connection = sqlite3.connect(dbname)
cursor = connection.cursor()

# FTP values
ftpEmail = 'username@unomaha.edu'
ftpTarget = ('/pub/taxonomy/','taxdump.tar.gz')
desiredFile = 'taxdump/names.dmp'

def main():
	dbLoad = False
	ftpNeeded = False

	try:
		options, arguments = getopt.getopt(sys.argv[1:],'d')
		for (opt,arg) in options:
			if opt == '-d':
				dbLoad = True
			elif opt == '-f':
				ftpNeeded =True
	except getopt.GetoptError as _:
		pass

	print(f'------------------------------------------\nBeginning Alignment Cleaning\n------------------------------------------')
	print(f'Parameters:\n  Database load:{dbLoad}\n  FTP Needed:\t{ftpNeeded}\n------------------------------------------')
	if ftpNeeded:
		print(f'Begin:\tHarvest FTP File {desiredFile}')
		# Create the temp file directory
		try:
			os.mkdir(datTempDirectory)
		except FileExistsError as _:
			shutil.rmtree(datTempDirectory)
			os.mkdir(datTempDirectory)
		# Collect the FTP File
		ftpHarvest(ftpTarget)
		print(f'Completed:\tFTP Harvest of target file {desiredFile}')
	if dbLoad:
		if not verifyFTPFiles(desiredFile):
			print(f'ERROR:\tFTP file {desiredFile} could not be verified in the temp directory')
			exit(1)
		# Create the database with the collected file
		print(f'Completed:\tVerification of FTP file {desiredFile}')
		print(f'\nBegin:\tInitializing DB')
		dbInit(dbname)
		print(f'Completed:\tInitializing DB')
		# Delete temp files
		shutil.rmtree(datTempDirectory)
		print(f'Completed removal of temp file directory at: {datTempDirectory}')
	# Verify the database is where I expect it to be, if so clean the alignment
	if not os.path.isfile(dbname):
		print(f'ERROR:\nCould not find database where expected, please run with option -d to initialize database')
		exit(1)
	else:
		print(f'\nBegin:\tCleaning Tree Alignment Taxonomy Names')
		cleanAlignments()
		print(f'Completed\tCleaning Tree Alignment Taxonomy Names')
	# Write the documentation to show how many taxons were available for each alignmnet
	print(f'Begin:\tWriting Documentation')
	writeAlignmentOut()
	print(f'Completed:\tWriting Documentation')


#--------------------------------------------#
#		Writing Documenttion Functions
#--------------------------------------------#
def writeAlignmentOut():
	with open(docDir + alignAvail,'w') as fileOut:
		for path, gene, count in availabilityOut:
			fileOut.write(f'{path},{gene},{count}\n')

#--------------------------------------------#
#		Cleaning Alignment Functions
#--------------------------------------------#
def cleanAlignments():
	for file in os.listdir(alignIn):
		speciesCount = cleanFile(alignIn,alignOut,file)

	pathway = file[:4]
	gene = file[5:]
	availabilityOut.append([pathway,gene,availabilityOut])

def cleanFile(directory,directoryOut,filename):
	print(f'\tCleaning File: {directory}{filename}')
	speciesCount = 0
	if not os.path.isfile(directory + filename):
		print(f'ERROR:\tCould not clean file {directory}{filename} because it was not in this location')
		exit(1)
	with open(directoryOut + filename, 'w') as fileOut:
		with open(directory + filename,'r') as fileIn:
			line = fileIn.readline()
			while line:
				if line[0] in ['(',':',')']: # Lines that wont have a organism name
					fileOut.write(line)
				else:
					speciesCount += 1
					fileOut.write(cleanSpecies(line))
				line = fileIn.readline()
	return speciesCount

def cleanSpecies(uncleanSpecies):
	# Uses Regex statement to get species name out of fasta description and compares it to entries in the taxonomy database
	cleanSpecies = ''
	# reSpecies assumes that there is a range in front of the species, this range is placed there for cases where the complete
	#		accession was not collected during the collection algorithm (as it was not needed). Adding it to the regex 
	#		statement allowed for more specificity when looking for the sequence.
	# Breakdown
	# .*			Catches the accession number and alignment number (assigned to each sequence during alignment by MAFFT)
	# _\d+_			Catches the version number for the accession (a reference point for cases where there is no range)
	# (?:\d+-\d+_)?	Catches the accession range (start and end index for the protein coding sequence)
	# ([A-Z][a-z]+_[a-z]+(?:_[a-z]+)?)	Capture the species name (may or may not be 3 elements to it)
	# (?:breed|unplaced|isolate|chromosome|strain|ecotype|linkage)	All the possible words that may follow the species name
	reSpecies = '.*_\d+_(?:\d+-\d+_)?([A-Z][a-z]+_[a-z]+(?:_[a-z]+)?)_(?:breed|unplaced|isolate|chromosome|strain|ecotype|linkage)'
	# re2Species is the same as reSpecies, except it accounts for cases where there is a cross between two species
	re2Species = '.*_\d+_(?:\d+-\d+_)?([A-Z][a-z]+_[a-z]+(?:_[a-z]+)?_x_[A-Z][a-z]+_[a-z]+(?:_[a-z]+)?)_(?:breed|unplaced|isolate|chromosome|strain|ecotype|linkage)'

	result = re.match(reSpecies,uncleanSpecies)
	if result:
		foundSpecies = result.group(1).replace('_',' ')
		sql = f'SELECT COUNT(*) FROM taxTable WHERE scientificName = \'{foundSpecies}\''
		count = int(cursor.execute(sql).fetchone()[0])
		if count == 1:
			cleanSpecies = foundSpecies
		elif count > 1:
			print(f'ERROR:\tSpecies {foundSpecies} returned more than one result in the database')
			exit(1)
		else: # Special Case for oddly named organism
			if foundSpecies == 'Caprimulgus carolinensis':
				cleanSpecies = 'Antrostomus carolinensis'
	if cleanSpecies != '':
		return cleanSpecies
	else:
		result = re.match(re2Species,uncleanSpecies) # Check for second regex case
		if result:
			foundSpecies = result.group(1)
			return foundSpecies
		else:
			print(f'ERROR:\tIssue cleaning species {uncleanSpecies.rstrip()}')
			exit(1)
			return uncleanSpecies


#---------------------------------------#
#		Database Loading Functions
#---------------------------------------#
def dbInit(dbname):

	# Initialize and truncate table
	sqls = [
		'CREATE TABLE IF NOT EXISTS taxTable (taxID INTEGER PRIMARY KEY, scientificName VARCHAR(50) NOT NULL, commonName VARCHAR(50))',
		'DELETE FROM taxTable'
	]
	for statement in sqls:
		cursor.execute(statement)
	connection.commit()

	loadSQL = 'INSERT INTO taxTable (taxID, scientificName, commonName)'

	# Taxonomy
	with open(datTempDirectory + desiredFile, 'r') as fileIn:
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
						print('ERROR:\t The TaxID: {} could not be loaded'.format(val))
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
					writeTaxValues(valueStrings, cursor, connection, loadSQL)
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
		writeTaxValues(valueStrings, cursor, connection, loadSQL)
		valueStrings = []

	print("Completed taxonomy load with {} records".format(currentTID))

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

def verifyFTPFiles(file):
	if not os.path.isfile(datTempDirectory + file):
		return False
	return True

def writeTaxValues(valueStrings, databaseC, databaseConn, loadSQL):
	# Writes values to the taxTable utilizing the partial SQLs made in the method to insert these values. Commits when the list is exausted
	#		Written on: 11 April 2018
	#		Author: Bob Anderson
	#	No Required methods (Required by parseTaxValues() method)
	# ARGUMENTS;
	#	valueStrings	-	List of values for insert strings in the format: " VALUES ({}, \'{}\', \'{}\')" No quotes on fina value if NULL
	#	databaseC 		-	Cursor from database described in parent method --> described in main()
	#	databaseConn	-	Connection object for the database (needed for commit)
	#	logFile			-	Handle for the log file associated with the run of the script
	print(f'\tCommiting {len(valueStrings)}')
	for valString in valueStrings:
		databaseC.execute(loadSQL + valString)
	databaseConn.commit()

# Utility Functions
def isInteger(val):
	try:
		val += 1
		return True
	except TypeError:
		return False

if __name__ == '__main__':
	main()