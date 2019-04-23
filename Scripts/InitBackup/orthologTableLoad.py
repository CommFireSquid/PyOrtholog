# Initial load script for the ortholog table pulled from the ncbi FTP site on 3/24/2018
#
#	Bob Anderson - Latest: 31 March 2018
#
# -------------------------------------
# UPDATES
#
#	1. Moved the commit out of the with and while sections to increase the transaction performance. This was tested
#			after the load was completed with the commit inside the loop. The testing was performed with an external
#			script with a similar design. It seems that the transactions are stored before the program closes.
#		-	Changed to commit with every 10,000 iterations to ensure the memory isnt flooded, remains untested.
#
#--------------------------------------
#
# This load strategy will be used on the much larger dataset for organisms.

import os
import ftplib
import shutil
import sqlite3

def main():
	target = 'gene_orthologs.gz'
	filename = './geneOrthologs.txt'
	database = './Data/dbLocal.sqlite3'

	ftpHarvest(target, filename)
	orthoLogIterWrite(database, filename)

def ftpHarvest(target, filename):
	fileOut = open(target, "wb")

	ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov','anonymous','robertjanderson@unomaha.edu')
	ftp.cwd('/gene/DATA/')

	ftp.retrbinary("RETR " + target, fileOut.write, 8*1024 )

	ftp.close()
	fileOut.close()

	shutil.unpack_archive(target)
	#os.system('gunzip {}'.format(target))
	shutil.move(target[:-3], "./Data/" + filename)
	#os.system('mv {} {}'.format(target[:-3], "./Data/" + filename))

def orthoLogIterWrite(database, filename):
	conn = sqlite3.connect(database)
	c = conn.cursor()

	c.execute("CREATE TABLE IF NOT EXISTS ncbi_ortholog (orthID integer PRIMARY KEY, taxID integer, geneID integer, relationship varchar(15), othTaxID integer, othGeneID integer)")


	# Find the next available PK value for the table and default to 0 if no value is found
	maxPKval = 0
	commitCTR = 0
	with open('../Data/' + filename, 'r') as fileIn:
		line = fileIn.readline()	# First line will be a header line that contains the names of the columns, not useful because we already know what the data will contain.
		line = fileIn.readline()

		while line:
			tax, gene, desc, otax, ogene = line.split('\t')
			maxPKval += 1
			c.execute("INSERT INTO ncbi_ortholog (orthID, taxID, geneID, relationship, othTaxID, othGeneID) VALUES ({},{},{},\'{}\',{},{})".format(maxPKval,tax,gene,desc,otax,ogene))
			
			commitCTR += 1
			if commitCTR % 10000 == 0:
				conn.commit()

			line = fileIn.readline()


	conn.commit() # Final commit to catch any remaining transactions in the ether
	os.remove("./Data/" + filename)
	#os.system('rm {}'.format("./Data/" + filename))

if __name__ == "__main__":
	main()