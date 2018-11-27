# Initial load script for the gene accession table pulled from the ncbi FTP site on 3/24/2018
#
#	Bob Anderson - Latest: 31 March 2018
#
# -------------------------------------
# UPDATES
# -------------------------------------



import sqlite3
import ftplib
import time

def main()
	ftpFile = 'gene2accession.gz'
	datfilename = 'gene2accession.txt'
	datfiledirectory = './Data/'
	startTime = time.time()

	# FTP Pull
	ftpHarvest(ftpFile, datfilename, datfiledirectory)

	conn = sqlite3.connect('./Data/dbLocal.sqlite3')
	c = conn.cursor()

	c.execute("CREATE TABLE IF NOT EXISTS ncbi_geneAccession (geneAccessionID INTEGER PRIMARY KEY, taxID INTEGER, geneID INTEGER, status VARCHAR(15), RNA_nucleotide_accession_version VARCHAR(20), protein_accession_version VARCHAR(20), protein_gi INTEGER, genomic_nucleotide_accession_version VARCHAR(20), genomic_nucleotide_gi INTEGER, start_position_on_the_genomic_accession INTEGER, end_position_on_the_genomic_accession INTEGER, orientation VARCHAR(1), assembly VARCHAR(20), mature_peptide_accession_version VARCHAR(20), mature_peptide_gi INTEGER, symbol VARCHAR(20)) ")
	c.execute('DELETE FROM ncbi_geneAccession')
	conn.commit()

	# Find the next available PK value for the table and default to 0 if no value is found
	maxPKval = 0
	commitCTR = 0	# Tracks the number of INSERT transactions. Used to commit every X transactions


	with open(datfiledirectory + datfilename 'r') as fileIn:
		line = fileIn.readline() # First line will be the table headers which I already have, so it will be ignored
		line = fileIn.readline()

		while line:
			maxPKval += 1
			
			# Replace - with NULL in input and assemble the value string
			hold = " VALUES ({}, ".format(maxPKval)
			for x in (line.rstrip()).split('\t'):
				if x == '-':
					hold = hold + "NULL, "
				elif isInteger(x):
					hold = hold + "{}, ".format(x)
				else:
					x = x.replace("'","''")
					hold = hold + "\'{}\', ".format(x)
			hold = hold[:-2] + ")"

			sql = "INSERT INTO ncbi_geneAccession (geneAccessionID, taxID, geneID, status, RNA_nucleotide_accession_version, RNA_nucleotide_accession_version, protein_accession_version, protein_gi, genomic_nucleotide_accession_version, genomic_nucleotide_gi, start_position_on_the_genomic_accession, end_position_on_the_genomic_accession, orientation, assembly, mature_peptide_accession_version, mature_peptide_gi, symbol)"
			c.execute( sql + hold )
							
			
			commitCTR += 1
			if commitCTR % 10000 == 0:
				conn.commit()
				timePassed = time.time() - startTime
				hours, rem = divmod(timePassed, 3600)
				minutes, seconds = divmod(rem, 60)
				print("Rows Committed: {} || Runtime: {:0>2}:{:0>2}:{:05.2f}".format(commitCTR,int(hours),int(minutes),seconds))
			
			line = fileIn.readline()

	conn.commit()
	os.remove(datfiledirectory + filename)
	#os.system('rm {}'.format(datfiledirectory + filename))

def isInteger(val):
	try:
		val += 1
		return True
	except TypeError:
		return False

def ftpHarvest(target, filename, fileDirectory):
	fileOut = open(target, "wb")

	ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov','anonymous','robertjanderson@unomaha.edu')
	ftp.cwd('/gene/DATA/')

	ftp.retrbinary("RETR " + target, fileOut.write, 8*1024 )

	ftp.close()
	fileOut.close()

	shutil.unpack_archive(target)
	#os.system('gunzip {}'.format(target))
	shutil.move(target[-3], fileDirectory + filename)
	#os.system('mv {} {}'.format(target[:-3], fileDirectory + filename))

