import re
import os

alignDir = 'Alignments/'
inDir = alignDir + 'Clean Alignments/'
outDir = alignDir + 'Data Files/'
docDir = 'Documentation/'
files = os.listdir(inDir)

NodeToParent = {} # Distance(value) from Node(key) to unnamed parent
NodeLinksParent = {} # Children ID or Name (key) to ParentID (value)
ParentLinksNode = {} # Parent ID (key) to List of direct Children IDs or Names

def main():
	for file in files:
		global NodeToParent
		global NodeLinksParent
		global ParentLinksNode

		NodeToParent = {}
		NodeLinksParent = {}
		ParentLinksNode = {}

		speciesInFile = [] # List of all species as leaves in file
		newickList = []	# File contents in a list 1 item = 1 line
		base = 'Homo sapiens' # Organism used as reference
		
		with open(inDir + file) as fileIn:
			# Read In all Possible Species
			# Read In file to list
			line = fileIn.readline()
			while line:
				newickList.append(line.rstrip())

				#Read Species From line if it is there
				reSpecies = '((?:[A-Z][a-z]+ [a-z]+(?: [a-z]+)?|[A-Z][a-z]+_[a-z]+_x_[A-Z][a-z]+_[a-z]+)):(\d\.\d{4,5})'
				result = re.search(reSpecies,line)
				if result:
					speciesInFile.append(result.group(1).rstrip())
					# Initial Dict Population
					#	Fill the dictionary with all leaves distance values to first parent
					NodeToParent[result.group(1).rstrip()] = float(result.group(2))
				line = fileIn.readline()

		

		# Find Distances for internal nodes to their parents
		newick = ''.join(newickList)
		parentID = 0
		#reSimpleGroup = '(\((\w+):\d\.\d{3,5},(\w+):\d\.\d{3,5}\):(\d\.\d{3,5}))'
		#reSimpleGroup = '(\(((?:[A-Z][a-z]+ [a-z]+(?: [a-z]+)?|[A-Z][a-z]+_[a-z]+_x_[A-Z][a-z]+_[a-z]+)):\d\.\d{3,5},((?:[A-Z][a-z]+ [a-z]+(?: [a-z]+)?|[A-Z][a-z]+_[a-z]+_x_[A-Z][a-z]+_[a-z]+)):\d\.\d{3,5}\):(\d\.\d{3,5}))'
		reSimpleGroup = '(\(((?:\w| )+):\d\.\d{3,5},((?:\w| )+):\d\.\d{3,5}\):(\d\.\d{3,5}))'
		results = re.findall(reSimpleGroup,newick)
		while results:
			for result in results:
				NodeToParent['Group{}'.format(parentID)] = float(result[3]) # New distance for the found parent
				NodeLinksParent[result[1]] = 'Group{}'.format(parentID)
				NodeLinksParent[result[2]] = 'Group{}'.format(parentID)
				ParentLinksNode['Group{}'.format(parentID)] = [result[1],result[2]]
				newick = newick.replace(result[0],'Group{}:0.00000'.format(parentID))
				parentID += 1
			results = re.findall(reSimpleGroup,newick)
		newick = newick[1:-1].split(',')
		if len(newick) == 2:
			child1 = (newick[0].split(':'))[0]
			child2 = (newick[1].split(':'))[0]
			NodeLinksParent[child1] = 'Root'
			NodeLinksParent[child2] = 'Root'
			ParentLinksNode['Root'] = [child1, child2]
		elif len(newick) == 3: # Weird case with some data ending with three remaining children
			# This must be a case where there is an odd number of distances
			#	Did not have time to properly investigate the implications of each item before the split
			child1 = (newick[0].split(':'))[0]
			child2 = (newick[1].split(':'))[0]
			child3 = (newick[2].split(':'))[0]
			NodeLinksParent[child1] = 'Group{}'.format(parentID)
			NodeLinksParent[child2] = 'Group{}'.format(parentID)
			ParentLinksNode['Group{}'.format(parentID)] = [child1, child2]
			NodeToParent['Group{}'.format(parentID)] = float((newick[1].split(':'))[1])
			NodeLinksParent['Group{}'.format(parentID)] = 'Root'
			NodeLinksParent[child3] = 'Root'
			ParentLinksNode['Root'] = ['Group{}'.format(parentID), child3]


		# Extend ParentLinksNode down to the roots
		for parent, children in ParentLinksNode.items():
			childrenReplaced = findLeaves(children)
			ParentLinksNode[parent] = childrenReplaced

		with open(docDir + 'NodeLinksParent.csv','w') as fileOut:
			for x,y in NodeLinksParent.items():
				fileOut.write('{},{}\n'.format(x,y))

		# COMPLETED BREAKING DOWN NEWICK
		#------------------------------#
		# Calculate distance between everything else and humans
		with open(outDir + file, 'w') as fileOut:
			for species in speciesInFile:
				if species != base:
					parent = lowestCommParent(species,base,species,base)
					distance = distanceBetween(species,parent) + distanceBetween(base,parent)
					fileOut.write('{},{}\n'.format(species,distance))
		print('Completed File: {}'.format(file))



def findLeaves(children):
	# Resolves all children to leaves <- Like a supervillan
	newChildren = []
	for child in children:
		reGroupChild = '(Group\d{1,3})'
		mChild = re.search(reGroupChild, child)
		if mChild:
			hold = findLeaves(ParentLinksNode[mChild.group(1)])
			for x in hold:
				newChildren.append(x)
		else:
			newChildren.append(child)
	return newChildren

def lowestCommParent(species1,species2,xParent,yParent):
	# Finds the lowest common parent for two given species
	#	Species will need to be given twice: 
	#		one set to retain as comparison
	#		one set to progress up the tree
	xParent = NodeLinksParent[xParent]
	yParent = NodeLinksParent[yParent]
	if species2 in ParentLinksNode[xParent]:
		return xParent
	elif species1 in ParentLinksNode[yParent]:
		return yParent
	else:
		return(lowestCommParent(species1,species2,xParent,yParent))

def distanceBetween(species,reference):
	sParent = NodeLinksParent[species]
	sDistance = NodeToParent[species]
	if sParent == reference:
		return sDistance
	else:
		return sDistance + distanceBetween(sParent,reference)

if __name__ == '__main__':
	main()