"""
blacklist.py
blacklist generation for ENCODE

Original Author: Alan Boyle / Python Conversion: Brydon P. G. Wall
Blacklist Copyright (c) 2018 Alan Boyle
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions.
"""
import math
import os
import argparse
import pysam
import statistics








class SequenceData():
	def __init__(self, bamFile: str = None, bamIndexFile: str = None) -> None:
		self.binsInput: list[int] = []
		self.binsSpikes: list[int] = []
		self.binsMultimapping: list[int] = []
		self.binsTemp: list[int] = []
		self.totalReads: int = 0
		self.zeroMulti: int = 0


		self.bamFile = bamFile
		self.bamIndexFile = bamIndexFile















	

	def getInputBins(self, mappability: list[int], binSize: int, binOverlap: int, refName: str) -> None:
	
		self.binsInput = []
		self.binsSpikes = []
		self.binsMultimapping = []
		self.binsTemp = []

		

		tempCounts: list[int] = [0 for i in range(len(mappability))]
		multiCounts: list[int] = [0 for i in range(len(mappability))]
		self.totalReads = 0
		testCntr: int = 0
		zmultiCntr: int = 0

		reader = pysam.AlignmentFile(filename=self.bamFile, index_filename=self.bamIndexFile)

		# This restricts output to a specific chromosome
		if reader.get_tid(refName) == -1:
			pass #print("No Reads to load!\n") #DEBUG
		else:



			# Keep vector of locations and counts of reads
			# filter on mappability based on each read length
			for al in reader.fetch(region=refName):
				if mappability[al.reference_start] > 0 and mappability[al.reference_start] <= al.reference_length: # Filter mappability
					tempCounts[al.reference_start] = tempCounts[al.reference_start] + 1; # Count reads at each position that are uniquely mappable
					if tempCounts[al.reference_start] == 1: # For lambda calculation
						testCntr += 1

				else: # A Multimapping read
					multiCounts[al.reference_start] = multiCounts[al.reference_start] + 1 # Count reads at each position that are not uniquely mappable
					zmultiCntr += 1
				self.totalReads += 1
		
		print(self.totalReads)

		reader.close()
		
		if zmultiCntr < 10000:
			self.zeroMulti += 1
			# raise AssertionError("No multimapping")












		readCntr: int  # collapsed reads per bin

		multiCntr: int # Multimapping reads per bin
		for i in range(0, len(tempCounts) - binSize, binOverlap):
			readCntr = 0

			multiCntr = 0
			for j in range(binSize):
				if tempCounts[i+j] > 0:
					readCntr += tempCounts[i+j]




				if multiCounts[i+j] > 0:
					multiCntr += multiCounts[i+j]


			
			self.binsInput.append(readCntr)

			self.binsMultimapping.append(multiCntr)
			self.binsTemp.append(0)











# Takes as input the mappability at each base and returns
# binned counts of mappability
def getMappabilityBins(retVal: list[int], mappability: list[int], binSize: int, binOverlap: int) -> list[int]:

	v: list[int] = []
	uniqueCntr: int
	uniqueLength: int = 36  # This is arbitraty and defines how long a read needs
							# to be to be considered unique
							# Should be set to something actually calculated in the uint8 files

	for i in range(0, (len(mappability) - binSize) + 1, binOverlap): 
		uniqueCntr = 0
		for j in range(binSize):
			if mappability[i+j] > 0 and mappability[i+j] <= uniqueLength:
				uniqueCntr += 1


		v.append(uniqueCntr)


	return v


def getdir(dirname: str, files: list[str], filetype: str) -> None:

    try:
        for entry in os.scandir(dirname):
            if entry.is_file() and entry.name.endswith(filetype):
                files.append(entry.name)
    except OSError as e:
		# could not open directory
        print("Unable to read input files!")
        exit(1)












def quantile(v, q):
    v.sort()
    n = len(v)
    h = (n - 1) * q
    i = int(h)  # Integer part of h
    f = h - i    # Fractional part of h

    if i + 1 >= n:
        return v[i]

    return v[i] + f * (v[i + 1] - v[i])




# Function to find rank
def rankify(A: list[float]) -> list[float]:
	n: int = len(A)
	R: list[float] = [0.0 for i in range(n)]
	T: list[tuple[float, int]] = []
	r: int = 1

	# Create array of tuples storing value and index
	for j in range(n):
		T.append((A[j], j))


	# Sort tuples by data value
	T.sort(key=lambda t: t[0])



	i: int = 0

	while i < n:
		j = i

		# Get elements of same rank
		while j < n - 1 and T[j][0] == T[j+1][0]:
			j += 1


		m: int = j - i + 1

		for j in range(m):
			# For each equal element use .5
			index = T[i + j][1]
			R[index] = r + (m - 1) * 0.5


		# Increment rank and index
		r += m
		i += m


	return R


def quantileNormalize(data: list[list[float]]) -> list[list[float]]:
	cellCount: int = len(data)
	binCount: int = len(data[0])

	# First calculate rank means
	rankedMean: list[float] = [0.0 for i in range(binCount)]
	for cellID in range(cellCount):
		x: list[float] = [0.0 for i in range(binCount)]
		for binID in range(binCount):
			x[binID] = data[cellID][binID]


		x.sort()

		for binID in range(binCount):
			rankedMean[binID] += x[binID]


	for binID in range(binCount):
		rankedMean[binID] /= float(cellCount)


	# calculate half value for ties
	rankedMeanTie: list[float] = [0.0 for i in range(binCount - 1)]
	for binID in range(binCount - 1):
		rankedMeanTie[binID] = float(rankedMean[binID] + rankedMean[binID + 1]) / 2


	# Iterate through each cell line
	for s in range(cellCount):
		bins: list[float] = [0.0 for i in range(binCount)]
		for p in range(binCount):
			bins[p] = data[s][p]

		bins = rankify(bins)

		binsQuantileNormalized: list[float] = [0.0 for i in range(binCount)]
		for p in range(binCount):
			if bins[p] % 1 != 0:
				binsQuantileNormalized[p] = rankedMeanTie[math.floor(bins[p]) - 1]
			else:
				binsQuantileNormalized[p] = rankedMean[int(bins[p] - 1)]


			data[s][p] = binsQuantileNormalized[p]
	return data






def getAbnormalRegions(inputData: list[SequenceData], binsMap: list[int], type_: int, normalize: bool) -> list[float]:
	"""Types:
	    1 - Read
	  	    - normalize = Reads / mapability
	    2 - Multimapping
	        - normalize = multimapping reads / total reads
	    3 - Spike (not used)"""


	result: list[float] = [0.0 for i in range(len(binsMap))]
	norm_temp: float = 0.0

	data: list[list[float]] = [[0.0] * len(binsMap) for i in range(len(inputData))]

	for j in range(len(inputData)): # process by column
		
		# Generalize this a little
		if type_ == 1: # reads
			inputData[j].binsTemp = inputData[j].binsInput
		elif type_ == 2: # multimapping
			inputData[j].binsTemp = inputData[j].binsMultimapping




		for i in range(len(binsMap)):  # track all rows in this column
			norm_temp = 0.0
			if normalize:
				if type == 2:  # multimapping
					if inputData[j].binsInput[i] != 0:
						norm_temp = inputData[j].binsTemp[i] / inputData[j].binsInput[i]
				elif binsMap[i] != 0:
					norm_temp = inputData[j].binsTemp[i] / binsMap[i]
			elif inputData[j].totalReads != 0:
				norm_temp = inputData[j].binsTemp[i] / inputData[j].totalReads * 1000000

			
			
			data[j][i] = norm_temp

	




	data = quantileNormalize(data)

	# Now we collapse rows
	means: list[float] = [0.0 for i in range(len(inputData))]
	for i in range(len(binsMap)): # over each row
		for j in range(len(inputData)): # over each column
			means[j] = data[j][i]
		
		result[i] = quantile(means, 0.5) # statistics.median(means) # This is median signal


	return result


def main(argc: int, argv: list[str]):

	mappability: list[int] = []
	binsMap: list[int] = []
	inputFileList: list[str] = []
	inputData: list[SequenceData] = []
	mappabilityFile: str = ""
	bamFile: str = ""
	bamIndexFile: str = ""
	refName: str = ""

	readNormList: list[float] = []
	multiList: list[float] = []

	# Parameters
	binSize: int = 1000
	binOverlap: int = 100

	if argc < 2:
		print("Blacklist is used to generate the ENCODE blacklists for various species.")
		print("Usage is ./Blacklist <chr>")
		print("The program requires an input/ folder containing indexed bam files.")
		print("The program requires a mappability/ folder containing Umap mappability files.")
		exit(0)
	else:
		refName = argv[1]


	# Load mappability file -- we need this multiple times
	mappabilityFile = "mappability/" + refName + ".uint8.unique"
	try:
		with open(mappabilityFile, "rb") as infile:
			mappability = [byte for byte in infile.read()]
	except IOError:
		print("Unable to read mappability files!")
		exit(1)






	binsMap = getMappabilityBins(binsMap, mappability, binSize, binOverlap)

	# Here we need to accumulate putative sites based on spikes and input levels
	getdir("input/", inputFileList, ".bam")
	for i in range(len(inputFileList)):
		bamFile = "input/" + inputFileList[i]
		bamIndexFile = "input/" + inputFileList[i] + ".bai"
		inputData.append(SequenceData(bamFile, bamIndexFile))

		inputData[-1].getInputBins(mappability, binSize, binOverlap, refName)


	# Fix error of one larger bins
	if(len(binsMap) > len(inputData[0].binsTemp)):
		binsMap.pop()


	readNormList = getAbnormalRegions(inputData, binsMap, 1, True)
	multiList = getAbnormalRegions(inputData, binsMap, 2, False)






	# Output the raw count information
	tempFirst: int
	tempLast: int
	inRegion: int = 0
	tempHit: int = 0
	hitCode: int = 0
	hitCounter: int = 0
	miss: int = 0
	regionClass: str = ""

	# Generate the threshold levels for weak and strong hits
	readWeakThresh: float = quantile(readNormList, 0.99)
	readStrongThresh: float = quantile(readNormList, 0.999)
	multiWeakThresh: float = quantile(multiList, 0.99)
	multiStrongThresh: float = quantile(multiList, 0.999)
	minThresh: float = min(readNormList) # We can bridge regions that have 0 signal as well.
			      # Optional to use flag these as blacklist.



	for i in range(len(binsMap)):
		if readNormList[i] >= readWeakThresh or multiList[i] >= multiWeakThresh or readNormList[i] <= minThresh:
			# tracking for overlapping regions
			miss = 0
			tempLast = i

			# If this is a new region, record it
			if inRegion == 0:
				tempFirst = i
				inRegion = 1


			# Check to see if this bin passes a threshold
			if readNormList[i] >= readStrongThresh:
				tempHit = 1
				hitCode = hitCode | 1
				hitCounter += 1
			

			if multiList[i] > multiStrongThresh:
				tempHit = 1
				hitCode = hitCode | 2
				hitCounter += 1
			
		else:
			if miss < (binSize // binOverlap) + 200: # bridge over adjacent bins plus 100 * 200 = 20kb
				miss += 1                           # recommend 5k for the smaller genomes
			else: # nothing in this distance
				inRegion = 0

				# If we hit a threshold we output the whole region
				if tempHit == 1:
					if hitCode == 2:
						regionClass = "Low Mappability"
					else:
						regionClass = "High Signal Region"

					print(f"{refName}\t{tempFirst*binOverlap}\t{tempLast*binOverlap + binSize}\t{regionClass}\n")
					tempHit = 0
					hitCode = 0
					hitCounter = 0
				
			
		
	

	# If we were in a region when when hit the end of the chromosome, output region
	# This may go past chromosome end!!
	if tempHit == 1:
		if hitCode == 2:
			regionClass = "Low Mappability"
		else:
			regionClass = "High Signal Region"

		print(f"{refName}\t{tempFirst*binOverlap}\t{tempLast*binOverlap + binSize}\t{regionClass}\n")	

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('args', nargs=argparse.REMAINDER)
	args = parser.parse_args()
	argv = ['blacklist.py']
	for item in args.args:
		argv.append(item)
	argc = len(argv)
	main(argc, argv)
