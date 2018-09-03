#!/bin/py -3
#############################################
#											#
#				  Objects 					#
#											#
#				Bob Anderson				#
#				 5 May 2018					#
#											#
#############################################
#
#	Updates
#--------------------------------------------
#
#
#--------------------------------------------
class GeneMember:
	# Holds all values and performs operations for a species within a gene
	def __init__(self,taxID):
		self.taxID = taxID
		self.geneID = 0
		self.accession = ''
		self.range_ = []
		self.complement = False
		self.fastaLocation = ''
		self.valueString = ''

class Gene:
	# Holds all the values associated with a gene and performs basic operations
	def __init__(self,geneName):
		self.geneName = geneName
		self.members = []
		self.orthologsList = []
	def addMember(self,member):
		# Adds a member to a list of gene members
		self.members.append(member)
	def addMembers(self,memberList):
		# Adds many members to a list of gene members
		for member in memberList:
			self.members.append(member)