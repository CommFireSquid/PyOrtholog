#!/bin/py -3
#############################################
#											#
#	NCBI DB Entrez Access for Gene Name 	#
#											#
#				Bob Anderson				#
#				1 April 2018				#
#											#
#############################################
#
#	Updates
#--------------------------------------------
#
#	- 31 March 2018		-	Constructed Query using taxID rather than organism name for easy relation to dbLocal.sqlite3
#	- 1 April 2018		-	Addition of regex to pull the accession out of the returned Gene Record
#	- 21 May 2018		-	Refactoring and application of multiprocessing
#	- 27 August 2018	-	Error corrections i.e. 
#								There were several syntax errors arrose seemingly out of nowhere such as a 
#								violation of entrez.esearch issues with a missing email that were not present before. 
#								Aditionally, it seems that there were some missing specifications when it came to the 
#								transition of general globals to GeneObject type variables. Due to poor refactoring 
#								by some fool -_-
#						-	Resolved issues with syntax and variable naming and eliminated database reclamation of
#								missing representative genes to resolve false database corruption crashes. Still no
#								complete run
#						+	Log files cannot be written by the child processes because the main process has a write 
#								lock on the log file.
#						+	Database is listed as corrupted when SELECT statments are run by the child processess
#								This may require changes to be made accounting for missing entries AFTER multiprocessing.
#								Luckily, accounts of the geneMembers are stored in the geneObject and written after
#									multiprocessing has completed. HOWEVER this may benefit from multiprocessing by
#									writing from the parent process as the children are gathering information?
#								Nightcap thought: Add the unkowns as geneMembers and process them differently... OR BETTER
#									looking into passing commands and results between child/parent processes.
#	- 28 August 2018	-	Resolved issues to get a complete run and database entries (even for unknown sequences)
#						+	Not getting the fastas for most results
#						+	Had to remove the foreign key reference for geneID to account for null entries on unknown seqs
#	- 29 August 2018	-	Rewrote the fasta creation process and representative gene lookup section using Biopython.
#								This allowed me to remove the sequenceTrimming() function that was inefficient and the
#									pullFasta() function that was unreliable in multiprocessing. I am not getting all files
#								Produced FASTA was checked for Mus musculus ACAN and it was identical, more analysis needed
#
#--------------------------------------------
# Process count = 25 hovers usage at <20% (<40% spikes) so increased to 50
# Process count = 50 Too many request error
#--------------------------------------------
import os
import re
import sys
import getopt
import sqlite3
import datetime
from shutil import rmtree # Used to remove a directory and all its sub contents
from Bio import Entrez, SeqIO
from multiprocessing import Pool
from OrthologObjects import Gene, GeneMember # Cool custom object for storing the results for a given gene and makes passing them all to a function easy later
#--------------------------------------------
# Global Values
#--------------------------------------------
# General Parameters
geneObject = Gene('None')
startTaxID = 9606
processCount = 40 # Number of processes to be spawned for the ortholog section
outputLineLength = 70 # Line widths for the final fasta outputs (ncbi standard is 70)
Entrez.email = ''

# Log file formatting
recordDivider = '*******************************************************'
altDivider = '-------------------------------'

# Global Directories and Files
logDirectory = './Logs/OrthoCollectLogs' # Where log files conventionally go
defaultErrorFile = './Logs/defaultError.log' # Used when conventions fail to create log file
fastaDirectory = './FASTA' # FASTA File Locations
logFile = '' # Global Log File

# Database Values
databaseName = './Data/dbLocal.sqlite3' # Database
connection = None
cursor = None
#--------------------------------------------

def main():
	global geneObject
	global Entrez
	startTime = datetime.datetime.now() # Time execution begins

	# Collect Command Line Arguments
	try:
		options, arguments = getopt.getopt(sys.argv[1:],"g:e:")
		for (opt, arg) in options:
			if (opt == '-g'):
				geneName = arg
			elif (opt == "-e"):
				Entrez.email = arg
	except getopt.GetoptError as OPT:
		usage(OPT)
	if geneName == '' or '*' in geneName:
		with open(defaultErrorFile,'w') as err:
			err.write("Below is the error resulting from fastaCollect.py on {}\n".format(str(startTime)))
			err.write("No Gene ID Argument Given: Please Use -g in the script call")
			exit(1)

	geneObject.geneName = geneName # Instantiate the Gene object with the supplied geneName parameter

	# Initialization of Storage Files
	createLogFile()	# Opening the global logFile "geneName.log"
	createFastaDir() # Create the FASTA file storage location (Delete if there is a directory under the gene name)
	connectDatabase() # Connecting to the database and creating global connection and cursor

	

	# Notify log file that process has begun
	logFile.write("Time of Execution: {} ||| Parameter Supplied: {}\n".format(str(startTime),geneName))
	logFile.write("{}\nBegin Initial Search TaxID = {}\n".format(recordDivider, startTaxID))

	# Begin Processing
	initialSearch() # Process the gene name for Homo sapiens
	now = datetime.datetime.now()
	logFile.write("Competion intial search for geneName = {}\n{}\nCurrent Time: {}\n".format(geneObject.geneName, recordDivider, str(now)))

	processOrthologs() # Generates a list of orthologs for multiprocessing

	# Condense all members and write their values into the database
	now = datetime.datetime.now()
	logFile.write("Completed Processing: Begin Write to DB Section\n{}\n{}\nCondensing Members\n".format(now, recordDivider))
	ctr = 0
	for _ in geneObject.members:
		condenseMember(ctr)
		ctr += 1
	insertIntoDB()

	now = datetime.datetime.now()
	logFile.write("{}\nProgram Completed at: {}\nRuntime: {}".format(recordDivider,now,now-startTime))

# Process Main Functions
def initialSearch():
	# Locates the Homo sapiens gene for the geneName and creates its gene member
	# 21 May 2018
	initialMember = GeneMember(startTaxID)
	initialMember.fastaLocation = "{}/{}_{}.fasta".format(fastaDirectory, geneObject.geneName, startTaxID)
	
	# Run a search on the NCBI gene database to find the gene matching passed name for the species with tax ID passed --> Should always use taxID = 9606 (Homo sapien) and then use this to get orthologs
	handle = Entrez.esearch(db="gene", retmax=10, 
		term="{}[Gene Name]) AND {}[Taxonomy ID]".format(geneObject.geneName, startTaxID), idtype="acc")
	record = Entrez.read(handle)
	logFile.write("{}\nInitial Entrez Gene DB Search (to get the geneID number for the symbol): \n{}\n{}\n".format(altDivider,altDivider,record))

	resultFound, initialMember.geneID, initialMember.accession, initialMember.range_, initialMember.complement = processFasta(record, initialMember.fastaLocation)

	handle.close()

	geneObject.addMember(initialMember) # Create a gene member object for humans

def processOrthologs():
	# Generates a list of orthologs for multiprocessing and executes it
	# 21 May 2018
	geneID = geneObject.members[0].geneID
	sql = "SELECT othTaxID, othGeneID FROM ncbi_ortholog WHERE geneID = {} AND taxID = \'{}\'".format(geneID, startTaxID)
	result = cursor.execute(sql)
	result = result.fetchall() # Fetches all orthologs for the given sequence
	now = datetime.datetime.now()
	logFile.write("Found {} Results In Ortholog Table\n{}\nBegin Multiprocessing\nAt time: {}\n".format(len(result),
		recordDivider,now)) # Write out the total results

	# Generate the process pool
	p = Pool(processes=processCount)
	orthologObjects = p.map(searchOrthologs, result)
	p.close()

	geneObject.addMembers(orthologObjects)

def searchOrthologs(result):
	# Takes a tuple (result) with tID and gID and processes it into a GeneMember object
	tID, gID = result
	geneMember = GeneMember(tID)
	geneMember.geneID = gID
	
	fastaDestination = "{}/{}_{}.fasta".format(fastaDirectory, geneObject.geneName, tID)
	geneMember.fastaLocation = "{}/{}_{}.fasta".format(fastaDirectory, geneObject.geneName, tID)
	
	resultFound, _, geneMember.accession, geneMember.range_, geneMember.complement = processFasta(gID, fastaDestination)

	if resultFound:
		logFile.write("Found annotated accession for gene id: {} tax id: {}".format(gID, tID))
	else:
		logFile.write("Found no annotated accession for gene id: {} tax id: {}".format(gID, tID))

	return geneMember


# FASTA/Entrez Utility Functions
def regexAccession(stringIn):
	accession = []
	# Uses Regex to identify an accession number in the stringIn
	#	stringIn	-	Entrez search result that is parsed to find the accession
	# 21 May 2018
	rAccession = re.compile("Annotation:\s.*\s(\w{2}_[\d\.]+) ?(\(?.*\)?)")
	rAccession = rAccession.search(stringIn)
	if rAccession:
		accession.append(rAccession.group(1))
		if rAccession.group(2) != '':
			accessionInd = (rAccession.group(2)).rstrip(')').lstrip('(') # Remove the parenthesis
			accessionInd = accessionInd.split('..') # Split on .. into a list
			accession.append(int(accessionInd[0]))
			try:	# Determine if the end of the accession end range contains the complement tag generally "#####, complement"
				accession.append(int(accessionInd[1])) # If the complement tag is present then this will fail
				logFile.write("\n{}\nAccession Found By Regex:\n{}\n{} in the range {}:{} Non-complement\n".format(altDivider,altDivider,accession[0],accession[1],accession[2]))
				return accession, False, True
			except ValueError as _: # Resolution of complement tags
				accession.append(int(accessionInd[1].split(',')[0])) # This will take just the accession number ignoring complement tags
				logFile.write("\n{}\nAccession Found By Regex:\n{}\n{} in the range {}:{} complement\n".format(altDivider,altDivider,accession[0],accession[1],accession[2]))
				return accession, True, True # Accounts for the complement tag with the first bool?
		else:
			logFile.write("\n{}\nAccession Found By Regex:\n{}\n{} no range identified\n".format(altDivider,altDivider,accession[0]))
			return accession, False, True
	else:
		return False, False, False

def processFasta(geneValue, fastaDestination):
# Takes in the gene ID (for orthologs) or annotation result for initial search, runs entrez search for the annotated representative gene,
#	then processes it into the trimmed FASTA file (and writes it).
# RELIES on regexAccession()
# geneValue	-	Integer ID number for gene (For orthologs where this is known) or an Entrez annotation for the initial search
#####
# Return true for found files and false for no result found and accession/range/complement
	foundAccession = None
	foundcomplement = None
	foundRange = []

	if isInteger(geneValue):
		geneID = geneValue
		with Entrez.efetch(db="gene", id=geneID, rettype="gb", retmode="text") as handle:
			result = handle.read()

			# Find the accession with regexAccession and process
			accession, complement, returned = regexAccession(result)
	else:
		# Get the Gene ID from the previous search and use it to find the gene record
		with Entrez.efetch(db="gene", id=geneValue.get("IdList",['No Result'])[0], rettype="gb", retmode="text") as handle:
			geneID = geneValue.get("IdList",['No Result'])[0]
			result = handle.read()
			logFile.write("\n{}\nEntrez Gene Record Fetch (To get the accession and range):\n{}\n{}".format(altDivider,altDivider,result))

			# Find the accession with regexAccession and process
			accession, complement, returned = regexAccession(result) # Currently not doing anything with returned

	# Universal fasta process #
	if returned: # Where an annotation for the entrez information was found
		foundAccession = accession[0]
		foundcomplement = complement
		if len(accession) == 3:
			foundRange = [accession[1], accession[2]]

		# End replicate now find a way to process the incoming FASTA
		with Entrez.efetch(db="nucleotide", id="{}".format(foundAccession), 
				rettype="fasta", retmode="text") as fastaHandle:
			for record in SeqIO.parse(fastaHandle, "fasta"):
				if complement:
					sequence = (record.seq).complement()
				else:
					sequence = record.seq
				with open(fastaDestination, 'w') as fastaOut:
					# Preprocess record - shorten to range and modify description
					if len(foundRange) == 2: # With accession range
						rDescript = re.compile("\w{2}_\d+\.*\d* (.*)")
						rDescript = rDescript.search(record.description.rstrip())
						if rDescript: # If regex can seperate the accession and the description
							fastaOut.write(">{}:{}-{} {}".format(record.id, foundRange[0], 
								foundRange[1], rDescript.group(1)))
						else:
							fastaOut.write(">{}".format(record.description))
						lineToWrite = str(sequence[foundRange[0]-1:foundRange[1]])
					else: # Without accession range
						fastaOut.write(">{}".format(record.description))
						lineToWrite = str(sequence)
					# Write sequence to file
					while len(lineToWrite) >= outputLineLength:	
						if len(lineToWrite) == outputLineLength:
							fastaOut.write('\n' + lineToWrite)
							lineToWrite = ""
						else:
							fastaOut.write('\n' + lineToWrite[:outputLineLength])
							lineToWrite = lineToWrite[outputLineLength:]
					if len(lineToWrite) > 0:
						fastaOut.write('\n' + lineToWrite) # Write the possibly remaining line < line length
					fastaOut.write('\n\n') # NCBI FASTAs come with new line at end of file
		return True, geneID, foundAccession, foundRange, complement
	else: # No entrez accession annotation
		return False, geneID, 'Unknown', [0, 0], False

# Database insert functions
def condenseMember(index):
	# Takes in an index from geneObjects.members and condenses its values into a value string for SQL
	#	index 	-	Integer index value for geneObjects.members
	# 21 May 2018
	hold = geneObject.members[index]
	
	# Condense Accession Number
	if len(hold.range_) == 2:
		if hold.complement:
			accVal = "{} ({} - {}, complement)".format(hold.accession, hold.range_[0], hold.range_[1]) #I am unsure why I was still using accession at this point accession[1], accession[2])
		else:
			accVal = "{} ({} - {})".format(hold.accession, hold.range_[0], hold.range_[1])
	else:
		accVal = hold.accession

	geneObject.members[index].valueString = [hold.taxID, "\'{}\'".format(accVal),
		"\'{}\'".format(hold.fastaLocation), hold.geneID] # removed hold.geneName from here because I dont think that I care about that and it is not currently populated in the GeneMember objects

def insertIntoDB():
	# Writes all values stored in geneObject.members into the database and commits
	#		Also init database connection and sequences table
	initSequencesTable() # Creates (if needed) the sequences table
	now = datetime.datetime.now()
	logFile.write('Writing {} values to database\nAt Time: {}\n'.format(len(geneObject.members), str(now)))
	partialSQL = 'INSERT INTO Sequences (seqID, taxID, accession, fastaLocation, geneID, geneName) '
	
	sqlMax = "SELECT MAX(seqID) FROM Sequences"
	seqID = (cursor.execute(sqlMax).fetchone())[0]
	if seqID == None:
		seqID = 0
	for member in geneObject.members:
		seqID += 1
		try:
			cursor.execute(partialSQL + 'VALUES ( {}, {}, {}, {}, {}, \'{}\' )'.format(seqID, member.valueString[0], 
				member.valueString[1], member.valueString[2], member.valueString[3], geneObject.geneName))
		except sqlite3.OperationalError as e:
			with open(defaultErrorFile,'w') as err:
				err.write("Error Occured at: {}".format(datetime.datetime.now()))
				err.write(partialSQL + 'VALUES ( {}, {}, {}, {}, {}, \'{}\' )'.format(seqID, member.valueString[0], 
				member.valueString[1], member.valueString[2], member.valueString[3], geneObject.geneName))
				err.write(e)
				exit(1)
	connection.commit()

# Storage File Init Functions
def connectDatabase():
	# Establishes the connection to the database and defines global values for connection and cursor
	# 21 May 2018
	global connection
	global cursor
	global exitCode
	try:
		connection = sqlite3.connect(databaseName)
		cursor = connection.cursor()
		return connection, cursor
	except Exception as _:
		logFile.write('Error Contacting the database')
		exit(1)

def createFastaDir():
	# Creates the directory structure as described in diagram i.e. FASTA/ACAN. If needed it will remove an existing directory (and contents)
	#		Then reassigns the global to the gene specific directory made
	# 21 May 2018
	global fastaDirectory
	geneFastaDirectory = fastaDirectory + '/{}'.format(geneObject.geneName)
	if not os.path.isdir(fastaDirectory):
		os.makedirs(fastaDirectory)
	if os.path.isdir(geneFastaDirectory):
		rmtree(geneFastaDirectory) # Remove this directory and all subcontents
	os.makedirs(geneFastaDirectory)	
	fastaDirectory = geneFastaDirectory	

def createLogFile():
	# Creates the log file in the logDirectory (creates logDirectory if it does not exist)
	# 21 May 2018
	global logFile
	global exitCode
	logName = geneObject.geneName + '.log'
	if not os.path.isdir(logDirectory):
		os.makedirs(logDirectory)
	try:
		logFile = open('{}/{}'.format(logDirectory,logName),'w')
	except OSError as e:
		with open(defaultErrorFile,'w') as err:
			err.write(e)
			err.write('\nError opening/creating the log file at: {}/{}'.format(logDirectory,logName))
		exit(1)

def initSequencesTable():
	# Create the Sequences Table if it does not exist and insert this record
	try:
		sql = "CREATE TABLE IF NOT EXISTS Sequences (seqID INTEGER PRIMARY KEY, taxID INTEGER NOT NULL, accession VARCHAR(75) NOT NULL, fastaLocation VARCHAR(100), geneID INTEGER, geneName VARCHAR(50) NOT NULL, CONSTRAINT fk_tax FOREIGN KEY (taxID) REFERENCES taxTable(taxID))" #, CONSTRAINT fk_orth FOREIGN KEY (geneID) REFERENCES ncbi_ortholog(geneID))
		cursor.execute(sql)
		sql = "DELETE FROM Sequences WHERE GeneName = \'{}\'".format(geneObject.geneName)
		cursor.execute(sql)
	except sqlite3.OperationalError:
		with open(defaultErrorFile,'w') as err:
			err.write("Error Occured at: {}".format(datetime.datetime.now()))
			err.write("Error in the following SQL:\n{}".format(sql))
			err.write(e)
			exit(1)
	connection.commit()

# Line Argument Error Message Function
def usage(errMessage):
	# Prints usage information to the log file if there is an issue retrieving the command arguments
	#		Writes to the default error as geneName is needed for regular error file
	#	errMessage	-	error message from getopt (could also be any string passed)
	# 21 May 2018
	time = str(datetime.datetime.now())
	with open(defaultErrorFile,'w') as err:
		err.write("Below is the error resulting from fastaCollect.py on {}\n".format(time))
		err.write("{}\n".format(errMessage))
		err.write("Usage:\tfastaCollect.py -g geneID")
	sys.exit(1)

def isInteger(value):
# Determines if the passed value is an integer and returns true if it is
	try:
		int(value)
		return True
	except (ValueError, TypeError) as _:
		return False

if __name__ == "__main__":
	main()