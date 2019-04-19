import re

def main():
	repo = []
	reRefSeq = '.*\(RefSeq\) (\w+), '
	reRefSeq2 = '.*\(RefSeq\) (\w+)'
	with open('genes.txt','r') as fileIn:
		line = fileIn.readline()
		while line:
			reMatch = re.match(reRefSeq,line)
			if reMatch:
				if reMatch.group(1) not in repo:
					repo.append('{}'.format(reMatch.group(1)))
				else:
						print('Skipped Duplicate: {}'.format(reMatch.group(1)))
			else:
				reMatch = re.match(reRefSeq2, line)
				if reMatch:
					if reMatch.group(1) not in repo:
						repo.append('{}'.format(reMatch.group(1)))
					else:
						print('Skipped Duplicate: {}'.format(reMatch.group(1)))
				else:
					print("Cannot Find: {}".format(line))
			line = fileIn.readline()

	with open('EtherLipidMetabolism.txt','w') as fileOut:
		for gene in repo:
			fileOut.write('{}\n'.format(gene))


if __name__ == '__main__':
	main()