import os

directoryIn = './FASTA'
directoryOut = './AlignmentFiles'

for directory in os.listdir(directoryIn):
	fileName = '{}_all.fasta'.format(directory)
	with open('{}/{}'.format(directoryOut,fileName),'w') as fileOut:
		for file in os.listdir('{}/{}'.format(directoryIn,directory)):
			with open('{}/{}/{}'.format(directoryIn,directory,file),'r') as fileIn:
				line = fileIn.readline()
				while line:
					fileOut.write(line)

					line = fileIn.readline()