import os

limit = 117

alignDir = 'Alignments/'
inDir = alignDir + 'Data Files/'
datDir = 'Data/'

Freq = {
	'GLY':{},
	'ELM':{},
	'FAB':{},
	'SLM':{},
	'TRY':{},
	'CTR':{}
}

for file in os.listdir(inDir):
	pathway = file.split('_')[0]
	taxCount = 0
	with open(inDir + file) as fileIn:
		line = fileIn.readline()
		while line:
			taxCount += 1
			line = fileIn.readline()
	if taxCount > limit:
		with open(inDir + file) as fileIn:
			line = fileIn.readline()
			while line:
				line = line.rstrip()
				line = line.split(',')
				if line[0] in Freq[pathway].keys():
					Freq[pathway][line[0]].append(float(line[1]))
				else:
					Freq[pathway][line[0]] = [float(line[1]),]



				line = fileIn.readline()

for pathway,fDict in Freq.items():
	outFile = ['Species,Average',]
	outFile2 = []
	outFile2Spec = ''
	datFileOut = f'{pathway}.csv'
	datFileOut2 = f'{pathway}_Variance.csv'
	for species, freqList in fDict.items():
		freqTotal = 0
		freqCount = 0
		outDist += f'{species}'
		outFile2Spec += ','
		for freq in freqList:
			freqTotal += freq
			freqCount += 1
			outDist += f',{freq}'

		outFile2.append(outDist)

		outFile.append(f'{species},{freqTotal/freqCount}')

	with open(datDir + datFileOut,'w') as fileOut:
		for line in outFile:
			fileOut.write(f'{line}\n')

	with open(datDir + datFileOut2, 'w') as fileOut:
		fileOut.write(outFile2Spec)
		for line in outFile:
			fileOut.write(f'{line}')