import os
def main():
	directories = os.listdir('./FASTA')
	for gene in directories:
		with open('./Data/{}_all2.fasta'.format(gene), 'w') as fastaOut:
			files = os.listdir('./FASTA/{}'.format(gene))
			for file in files:
				with open('./FASTA/{}/{}'.format(gene, file), 'r') as fastaIn:
					line = fastaIn.readline()
					while line:
						fastaOut.write(line)
						line = fastaIn.readline()

if __name__ == '__main__':
	main()