import os
import ftplib
import sqlite3
import datetime

lineDivider = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
tempFolder = 'TempTarFolder/'
dataDirectory = 'Data/'

def main():
	targetArchive = 'taxdump.tar.gz'
	targetFile = 'names.dmp'
	filename = 'taxNames.txt'
	database = './Data/dbLocal.sqlite3'

	logFile = open('./Logs/init_taxLoad.log', 'w')

	# Processes
	now = datetime.datetime.now()
	logFile.write('Begin Log: FTP pull and dbLocal load of taxInformation\nTime: {}\n{}\n'.format(str(now), lineDivider))
	ftpHarvest_TarGz(targetArchive, targetFile, filename, logFile)
	now = datetime.datetime.now()
	logFile.write('{}\nCompleted dearchiving and file removal at: {}\n'.format(lineDivider, str(now)))
	logFile.write('File: {} located in the data directory: {} ready for load processing\n'.format(filename, dataDirectory))
	parseTaxNames(dataDirectory + filename, database, logFile)

	logFile.close()

def ftpHarvest_TarGz(targetArchive, targetFile, filename, logFile):
	# FTP Harvest variant that pulls .tar.gz files. Specifically handles the useless files that come in the archive utilizing a ./temp folder that will be deleted after the target has been taken out
	#		Written on: 11 April 2018
	#		Author: Bob Anderson
	#	No required methods - All handled "in house"
	# ARGUMENTS:
	#	targetArchive	-	Archive (.tar.gz) of interest as it is named on NCBI's FTP site
	#	targetFile		-	File in the targetArchive that is of interest
	#	filename		-	Desired final file name for the (usually .txt) output that will be placed in ../Data/filename (from the script in Initial load scripts folder)
	ftpSite = 'ftp.ncbi.nlm.nih.gov'
	ftpUser = 'anonymous'
	ftpEmailPass = 'robertjanderson@unomaha.edu'
	ftpCWD = '/pub/taxonomy/'

	now = datetime.datetime.now()

	logFile.write('Current Time: {}\nFTP Site: {}\nFTP User: {}\nFTP Pass (email): {}\nFTP Working Directory: {}\n'.format(str(now),ftpSite,ftpUser,ftpEmailPass,ftpCWD))
	logFile.write('{}\n'.format(lineDivider))

	os.system('mkdir {}'.format(tempFolder))
	fileOut = open(tempFolder + targetArchive, "wb")

	ftp = ftplib.FTP(ftpSite,ftpUser,ftpEmailPass)
	ftp.cwd(ftpCWD)

	ftp.retrbinary("RETR " + targetArchive, fileOut.write, 8*1024 )

	ftp.close()
	fileOut.close()

	# Decompress, Dearchive, and Extract the file I need from temp folder and delete the others
	os.system('gunzip {}'.format(tempFolder + targetArchive)) # Unzip the tar archive
	os.system('tar -xf {} -C {}'.format(tempFolder + targetArchive[:-3], tempFolder)) # Dearchive the tar file into the temp directory
	os.system('[ -f {} ] && rm {} '.format(dataDirectory + filename, dataDirectory + filename)) # If there is a file where we want to move the names.dmp delete it
	os.system('mv {} {}'.format(tempFolder + targetFile, dataDirectory + filename))	# Move the names.dmp file
	os.system('[ -d {} ] && rm -r {}'.format(tempFolder,tempFolder))	# if the temp directory exists (it does) delete it recursively


def parseTaxNames(filename, database, logFile):
	# Designed to iterate through names.dmp for loading the taxTable into the database
	#		Written on: 11 April 2018
	#		Author: Bob Anderon
	# 	Requires writeValues() function
	# ARGUMENTS:
	#	filename	-	Name of the fileIn (usually FTP pulled tax file)
	#	database 	-	Path to the database that needs the taxTable created and/or populated
	logFile.write('{}\nBegin Taxonomy File Processing\nTarget database: {}\n{}\n'.format(lineDivider, database,lineDivider))
	fileIn = open(filename, "r")

	databaseConn = sqlite3.connect(database)
	c = databaseConn.cursor()
	c.execute('CREATE TABLE IF NOT EXISTS taxTable (taxID INTEGER PRIMARY KEY, scientificName VARCHAR(50) NOT NULL, commonName VARCHAR(50))')

	valueStrings = []
	commitSize = 10000

	# Values collected for a particular TID
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
					logFile.write('Error handling TaxID: {}\n'.format(val))
					exit(0)
			else:
				values[ctr] = val.rstrip().lstrip().replace("'", "`") # Temporary solution to cases where there are single quotes in the species names
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
				writeTaxValues(valueStrings, c, databaseConn, logFile)
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
		logFile.write("Bad TaxID: {} --> No Scientific Name Found\n".format(currentTID))

	if len(valueStrings) > 0:
		writeTaxValues(valueStrings, c, databaseConn, logFile)
		valueStrings = []


def writeTaxValues(valueStrings, databaseC, databaseConn, logFile):
	# Writes values to the taxTable utilizing the partial SQLs made in the method to insert these values. Commits when the list is exausted
	#		Written on: 11 April 2018
	#		Author: Bob Anderson
	#	No Required methods (Required by parseTaxValues() method)
	# ARGUMENTS;
	#	valueStrings	-	List of values for insert strings in the format: " VALUES ({}, \'{}\', \'{}\')" No quotes on fina value if NULL
	#	databaseC 		-	Cursor from database described in parent method --> described in main()
	#	databaseConn	-	Connection object for the database (needed for commit)
	#	logFile			-	Handle for the log file associated with the run of the script
	now = datetime.datetime.now()
	logFile.write('Writing {} values to database\nAt Time: {}\n'.format(len(valueStrings), str(now)))
	partialSQL = 'INSERT INTO taxTable (taxID, scientificName, commonName) '
	for valString in valueStrings:
		databaseC.execute(partialSQL + valString)
	databaseConn.commit()
	

if __name__ == "__main__":
	main()