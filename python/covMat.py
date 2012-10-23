import numpy as np

class covMat():

	def __init__(self, aname):

		self.matrix = []
		self.paramNames=[]
		self.size =0
		self.name=aname
		self.comment=''

	def paramNameString(self):
		return " ".join(self.paramNames)

	def loadFromFile(self, filename):
		textFileHandle = open(filename)
		textFileLines=textFileHandle.readlines()
		textFileHandle.close()
		first=textFileLines[0].strip()
		if first.startswith('#'):
			paramNames=first[1:].split()
			self.size=len(paramNames)
		else:
			raise Exception('.covmat must now have parameter names header')
		matrix = [[0 for col in range(self.size)] for row in range(self.size)]
		used = []
		for i in range(self.size):
			splitLine = textFileLines[i+1].split()
			for j in range(len(splitLine)):
			 	matrix[i][j] = float(splitLine[j])
			if matrix[i].count(0.) != self.size:
				used.append(i)

		self.size = len(used)
		self.matrix = np.empty((self.size,self.size))
		self.paramNames=[]
		for i in range(self.size):
			self.paramNames.append(paramNames[used[i]])
			for j in range(self.size):		
				self.matrix[i,j] = matrix[used[i]][used[j]]
	

