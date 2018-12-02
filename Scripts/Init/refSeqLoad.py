# Initial load script for the refseq table pulled from the ncbi FTP site
#
#	Bob Anderson - Latest: 26 November 2018
#

import os
import ftplib
import getopt
import sqlite3

ftpEmail = ''
target = 'gene2refseq.gz'
datfiledirectory = './Data/'
datfilename = './gene2refseq.txt'
database = './Data/dbLocal.sqlite3'

def main():
	global ftpEmail

	try:
		options, arguments = getopt.getopt(sys.argv[1:],"e:")
		for (opt, arg) in options:
			if (opt == "-e"):
				ftpEmail = arg
	except getopt.GetoptError as OPT:
		usage(OPT)

	ftpHarvest(target, datfilename)
	orthoLogIterWrite(database, datfilename)

def writeToDB(database, filename):
	connection = sqlite3.connect(database)
	cursor = connection.cursor()

	sql = 'CREATE TABLE IF NOT EXISTS refseq (taxID INTEGER, geneID INTEGER PRIMARY KEY, status VARCHAR(15), ' +
		'RNA_nucleotide_Accession_Version VARCHAR(20), RNA_Nucleotide_GI INTEGER, protein_accession_version VARCHAR(20), ' +
		'protein_gi INTEGER, genomic_nucleotide_accession_version VARCHAR(20), genomic_nucleotide_gi INTEGER, ' +
		'start_position_on_the_genomic_accession INTEGER, end_position_on_the_genomic_accession INTEGER, ' +
		'orientation VARCHAR(1), assembly VARCHAR(20), mature_peptide_accession_version VARCHAR(20), ' +
		'mature_peptide_gi INTEGER, symbol VARCHAR(20))'
	cursor.execute(sql)
	cursor.execute('DELETE FROM refseq')
	conn.commit()

	commitCTR = 0
	with open(datfiledirectory + datfilename 'r') as fileIn:
		line = fileIn.readline() # First line will be the table headers which I already have, so it will be ignored
		line = fileIn.readline()

		while line:
			hold = " VALUES ("
			for x in (line.rstrip()).split('\t'):
				if x == '-':
					hole = hold + "NULL, "
				elif isInteger(x):
					hold = hold + "{}, ".format(x)
				else:
					x = x.replace("'", "''")
					hold = hold + "\'{}\', ".format(x)
			hold = hold[:-2] + ")"

			sql = "INSERT INTO refseq (taxID, geneID, status, RNA_nucleotide_accession_version, RNA_nucleotide_accession_version, protein_accession_version, protein_gi, genomic_nucleotide_accession_version, genomic_nucleotide_gi, start_position_on_the_genomic_accession, end_position_on_the_genomic_accession, orientation, assembly, mature_peptide_accession_version, mature_peptide_gi, symbol)"
			cursor.execute(sql + hold)

			commitCTR += 1
			if commitCTR % 10000 == 0:
				connection.commit()

			line = fileIn.readline()

	connection.commit()
	os.remove(datfiledirectory + filename)

def ftpHarvest(target, filename):
	fileOut = open(target, "wb")

	ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov','anonymous',ftpEmail)
	ftp.cwd('/gene/DATA/')

	ftp.retrbinary("RETR " + target, fileOut.write, 8*1024 )

	ftp.close()
	fileOut.close()

	shutil.unpack_archive(target)
	#os.system('gunzip {}'.format(target))
	shutil.move(target[:-3], datfiledirectory + filename)
	#os.system('mv {} {}'.format(target[:-3], "./Data/" + filename))

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
	time = str(datetime.datetime.now())
	with open(defaultErrorFile,'w') as err:
		err.write("Below is the error resulting from fastaCollect.py on {}\n".format(time))
		err.write("{}\n".format(errMessage))
		err.write("Usage:\tfastaCollect.py -e email@email.com")
	sys.exit(1)