#!/usr/bin/env python

#Author: Jose Lugo-Martinez
#File: generateCGBayesNetsDataset.py
#Date: December 12, 2017
#Advisor Profs. Ziv Bar-Joseph and Giri Narasimhan
#Description: This function creates sample data sets for constructing a bayesian network and dynamic bayesian notworks with and without time-warping.

#Last Modified: February 18, 2021

#Sample call:
#python generateCGBayesNetsDataset.py ./infant_gut/infant_gut_microbiota_annotated.tsv ./infant_gut/selectedTaxa.tsv ./infant_gut/infant_gut_microbiota_optimal_alignment.tsv infant_gut_microbiota_noalignment infant_gut_microbiota_alignment

import sys, copy, math
import numpy as np
from scipy import interpolate

TAXA_OFFSET = 10 #Index offset in dataset where taxon data starts
MINIMUN_NUMBER_MEASURED_TIMEPOINTS = 9 #Minimum number of required measured timepoints per subject/sample
SAMPLING_RATE = 3 #1, 3, 5, 7, 14
UPPER_BOUND = 1.0 #Maximum value in relation to relative abundance for bounding continuous representation range
LOWER_BOUND = 0.0 #Minimum value in relation to relative abundance for bounding continuous representation range

def warpFunction(a, b, s, warpType='linear'):
	if warpType == 'exponential':
		return np.exp((s - b) / a)
	else:
		return (s - b) / a

def warpFunctionInverse(a, b, t, warpType='linear'):
	if warpType == 'exponential':
		return a * np.log(t) + b
	else:
		return (a * t) + b

def getTaxaList(taxaListFilename):
	try:
		#Open input file
		infile = open(taxaListFilename, "r")
	except(IOError), e:
		print "<<ERROR>> Unable to open the file", taxaListFilename, "\nThis program will be quiting now.", e
		sys.exit()
	
	selectedTaxa = []
	for line in infile:
		line = line.strip()
		currTaxa = line.split('\t')[0]
		if not (currTaxa in selectedTaxa):
			selectedTaxa.append(currTaxa)
	#Close file
	infile.close()

	return selectedTaxa

def getSamples(dataFilename, selectedTaxa, alignmentFilename, outfilenameDBNnoAlignment, outfilenameDBNalignment):
	try:
		#Open input file
		infile = open(alignmentFilename, "r")
	except(IOError), e:
		print "<<ERROR>> Unable to open the file", alignmentFilename, "\nThis program will be quiting now.", e
		sys.exit()

	headers = infile.readline().strip().split('\t')
	
	alignmentParameters = {}
	for line in infile:
		line = line.strip()
		tokens = line.split('\t')
		referenceSampleID = tokens[0]
		currentSampleID = tokens[1]
		if len(tokens) < 3:
			continue
		globalAlignmentError = tokens[2]
		a = float(tokens[3])
		b = float(tokens[4])
		alpha = float(tokens[5])
		beta = float(tokens[6])
		if not (currentSampleID in alignmentParameters):
			alignmentParameters[currentSampleID] = (a, b, alpha, beta)
		else:
			print "Duplicate sample", currentSampleID
	#Close file
	infile.close()

#	print len(alignmentParameters)

	try:
		#Open input file
		infile = open(dataFilename, "r")
	except(IOError), e:
		print "<<ERROR>> Unable to open the file", dataFilename, "\nThis program will be quiting now.", e
		sys.exit()

	headers = infile.readline().strip().split('\t')
	taxaNames = copy.copy(headers[TAXA_OFFSET:])

	headersLine = 'SubjectID' + '\t' + 'Day of life sample obtained' + '\t' + 'Gestational age at birth' + '\t' + 'Postconceptional age sample obtained' + '\t' + 'Gender' + '\t' + 'Mode of birth' + '\t' + 'Room type' + '\t' + 'Human milk used' + '\t' + 'Days of antibiotics'
	for taxaName in selectedTaxa:
		headersLine += '\t' + taxaName
	headersLine += '\n'

	subjectIDs = []
	samplesPerSubject = {}
	samplesPerSubjectInfo = {}
	gestationalAgePerSubject = {}
	genderPerSubject = {}
	modeOfBirthPerSubject = {}
	periodOfStudyPerSubject = {}
	roomTypePerSubject = {}
	humanMilkUsedPerSubject = {}
	daysOfAntibioticsPerSubject = {}
	currDaysOfAntibiotics = []
	currSubjectSample = {}
	previousSubjectID = ''
	#Iterate over file
	for line in infile:
		tokens = line.split('\t')
		subjectID = tokens[0]
		if not (subjectID in subjectIDs):
			subjectIDs.append(subjectID)
		if previousSubjectID != '' and previousSubjectID != subjectID:
			daysOfAntibioticsPerSubject[previousSubjectID] = copy.copy(currDaysOfAntibiotics)
			samplesPerSubject[previousSubjectID] = copy.copy(currSubjectSample)
			sampleSizePerSubject = len(currSubjectSample[taxaNames[0]])
			samplesPerSubjectInfo[previousSubjectID] = (dayOfLifeFirstSample, dayOfLifeLastSample, sampleSizePerSubject)
			currDaysOfAntibiotics = []
			currSubjectSample = {}
		dayOfLifeSampleObtained = float(tokens[1])
		gestationalAgeAtBirth = float(tokens[2])
		if not (subjectID in gestationalAgePerSubject):
			gestationalAgePerSubject[subjectID] = gestationalAgeAtBirth
		postconceptionalAgeSampleObtained = float(tokens[3])
		gender = tokens[4]
		if not (subjectID in genderPerSubject):
			genderPerSubject[subjectID] = gender
		modeOfBirth = tokens[5]
		if not (subjectID in modeOfBirthPerSubject):
			modeOfBirthPerSubject[subjectID] = modeOfBirth
		periodOfStudy = tokens[6]
		if not (subjectID in periodOfStudyPerSubject):
			periodOfStudyPerSubject[subjectID] = periodOfStudy
		roomType = tokens[7] #0-Single room; 1-
		if not (subjectID in roomTypePerSubject):
			roomTypePerSubject[subjectID] = roomType
		humanMilkUsed = tokens[8] #0-none, 1-<10%, 2-between 10% and 50%, 3->50%
		if not (subjectID in humanMilkUsedPerSubject):
			humanMilkUsedPerSubject[subjectID] = humanMilkUsed
		daysOfAntibiotics = float(tokens[9])
		currDaysOfAntibiotics.append((dayOfLifeSampleObtained, daysOfAntibiotics))
		currentAbundancePerTaxa = copy.copy(tokens[TAXA_OFFSET:])
		for taxaIndex in xrange(len (taxaNames)):
			taxaName = taxaNames[taxaIndex]
			abundance = float(currentAbundancePerTaxa[taxaIndex].strip())
			if not (taxaName in currSubjectSample):
				currSubjectSample[taxaName] = [(dayOfLifeSampleObtained, abundance)]
				dayOfLifeFirstSample = dayOfLifeSampleObtained
			else:
				currSubjectSample[taxaName].append((dayOfLifeSampleObtained, abundance))
				dayOfLifeLastSample = dayOfLifeSampleObtained
				
		previousSubjectID = subjectID
	#Close file
	infile.close()

	daysOfAntibioticsPerSubject[previousSubjectID] = copy.copy(currDaysOfAntibiotics)
	samplesPerSubject[previousSubjectID] = copy.copy(currSubjectSample)
	sampleSizePerSubject = len(currSubjectSample[taxaNames[0]])
	samplesPerSubjectInfo[previousSubjectID] = (dayOfLifeFirstSample, dayOfLifeLastSample, sampleSizePerSubject)

	daysOfAntibioticsSplinesPerSubject = {}
	taxonSplinesPerSubject = {}
	selectedTaxaByMean = {}
	for subjectID in subjectIDs:
		abundanceByTaxa = samplesPerSubject[subjectID]
		splinesPerTaxa = {}
		dayFirstSample, dayLastSample, numSamples = samplesPerSubjectInfo[subjectID]
		if numSamples < MINIMUN_NUMBER_MEASURED_TIMEPOINTS:
			del samplesPerSubject[subjectID]
			continue
		timepoints = []
		daysOfAntibioticsPerSample = []
		#Get spline for days of antibiotics across timepoints. NOTE: This is not a continuous function!
		for timepoint, daysOfAntibiotics in daysOfAntibioticsPerSubject[subjectID]:
			timepoints.append(float(timepoint))
			daysOfAntibioticsPerSample.append(float(daysOfAntibiotics))
		#Use B-spline to extrapolate values. NOTE: Parameter s must be adjusted appropriately to avoid over-fitting. 
		tck = interpolate.splrep(timepoints, daysOfAntibioticsPerSample, s=0.05, xb=dayFirstSample, xe=dayLastSample)
		daysOfAntibioticsSplinesPerSubject[subjectID] = copy.copy(tck)
		#Get splines for each taxa across timepoints
		for taxaName, abundanceLevelPerTimepoint in abundanceByTaxa.iteritems():
			timepoints = []
			relativeAbundances = []
			for timepoint, abundance in abundanceLevelPerTimepoint:
				timepoints.append(timepoint)
				relativeAbundances.append(abundance)
			mean = getMean(relativeAbundances)
			variance = getVariance(relativeAbundances)
			#Use B-spline to extrapolate values. NOTE: Parameter s must be adjusted appropriately to avoid over-fitting.
			tck = interpolate.splrep(timepoints, relativeAbundances, k=3, s=0.001, xb=dayFirstSample, xe=dayLastSample)
			splinesPerTaxa[taxaName] = copy.copy(tck)
		taxonSplinesPerSubject[subjectID] = copy.copy(splinesPerTaxa)

	outfileSuffix = '_sr' + str(SAMPLING_RATE) + 'd' + '.tsv'
	outfilenameDBNalignment = outfilenameDBNalignment + outfileSuffix
	outfilenameDBNnoAlignment = outfilenameDBNnoAlignment + outfileSuffix
	outfileAligned = open(outfilenameDBNalignment, 'w')
	outfileUnaligned = open(outfilenameDBNnoAlignment, 'w')
	outfileAligned.writelines(headersLine)
	outfileUnaligned.writelines(headersLine)
	referenceSample = samplesPerSubject[referenceSampleID]
	referenceSampleDaysOfAntibioticsSplineParameters = daysOfAntibioticsSplinesPerSubject[referenceSampleID]
	referenceSampleTaxonSplineParameters = taxonSplinesPerSubject[referenceSampleID]
	dayFirstSample, dayLastSample, numSamples = samplesPerSubjectInfo[referenceSampleID]
	abundanceLevelPerTimepoint = referenceSample[taxaNames[0]]
	timepointsReferenceSample = []
	for timepoint, abundance in abundanceLevelPerTimepoint:
		timepointsReferenceSample.append(timepoint)
	#Add day of first sample obtained as a timepoint via continuous representation
	if (not (dayFirstSample in timepointsReferenceSample)) or timepointsReferenceSample[0] != dayFirstSample:
		timepointsReferenceSample.insert(0, dayFirstSample)
	#Add day of last sample obtained as a timepoint via continuous representation
	if (not (dayLastSample in timepointsReferenceSample)) or timepointsReferenceSample[-1] != dayLastSample:
		timepointsReferenceSample.insert(-1, dayLastSample)

	referenceSampleStartTimepoint = math.ceil(timepointsReferenceSample[0])
	referenceSampleEndTimepoint = math.floor(timepointsReferenceSample[-1])
	sampleLength = referenceSampleEndTimepoint - referenceSampleStartTimepoint + 1.0
	timepointsReferenceSample = np.arange(referenceSampleStartTimepoint, (referenceSampleEndTimepoint + 1.0), SAMPLING_RATE)

	for currentSampleID in subjectIDs:
		if currentSampleID == referenceSampleID:
			for currentTimepoint in timepointsReferenceSample:
				daysOfAntibiotics = int(round(interpolate.splev(currentTimepoint, referenceSampleDaysOfAntibioticsSplineParameters)))
				gestationalAgeAtBirth = gestationalAgePerSubject[referenceSampleID]
				postconceptionalAgeSampleObtained = (currentTimepoint/7.0) + gestationalAgeAtBirth
				gender = genderPerSubject[referenceSampleID]
				modeOfBirth = modeOfBirthPerSubject[referenceSampleID]
				periodOfStudy = periodOfStudyPerSubject[referenceSampleID]
				roomType = roomTypePerSubject[referenceSampleID]
				humanMilkUsed = humanMilkUsedPerSubject[referenceSampleID]
				if daysOfAntibiotics < 0:
					daysOfAntibiotics = 0
				outline = referenceSampleID + '\t' + str(currentTimepoint) + '\t' + str(gestationalAgeAtBirth) + '\t' + str(postconceptionalAgeSampleObtained) + '\t' + gender + '\t' + modeOfBirth + '\t' + roomType + '\t' + humanMilkUsed + '\t' + str(daysOfAntibiotics)
				#Save the taxon relative abundance values for this timepoint 
				taxonAbundanceValues = {}
				for taxaName, abundanceLevelPerTimepoint in referenceSample.iteritems():
					abundanceValue = interpolate.splev(currentTimepoint, referenceSampleTaxonSplineParameters[taxaName])
					abundanceValue = truncateAbundanceValue(abundanceValue)
					taxonAbundanceValues[taxaName] = abundanceValue
				#Output normalized relative abundance values for the selected taxa
				for taxaName in selectedTaxa:
					normalizedAbundanceValue = taxonAbundanceValues[taxaName] / sum(taxonAbundanceValues.values())
					outline += '\t' + str(normalizedAbundanceValue)
				outline += '\n'
				outfileAligned.writelines(outline)
				outfileUnaligned.writelines(outline)
		else:
			if not (currentSampleID in alignmentParameters):
				continue
			a, b, alpha, beta = alignmentParameters[currentSampleID]
			currentSample = samplesPerSubject[currentSampleID]
			currentSampleDaysOfAntibioticsSplineParameters = daysOfAntibioticsSplinesPerSubject[currentSampleID]
			currentSampleTaxonSplineParameters = taxonSplinesPerSubject[currentSampleID]
			dayFirstSample, dayLastSample, numSamples = samplesPerSubjectInfo[currentSampleID]
			abundanceLevelPerTimepoint = currentSample[taxaNames[0]]
			timepointsCurrentSample = []
			for timepoint, abundance in abundanceLevelPerTimepoint:
				timepointsCurrentSample.append(timepoint)
			#Add day of first sample obtained as a timepoint via continuous representation
			if (not (dayFirstSample in timepointsCurrentSample)) or timepointsCurrentSample[0] != dayFirstSample:
				timepointsCurrentSample.insert(0, dayFirstSample)
			#Add day of last sample obtained as a timepoint via continuous representation
			if (not (dayLastSample in timepointsCurrentSample)) or timepointsCurrentSample[-1] != dayLastSample:
				timepointsCurrentSample.insert(-1, dayLastSample)

			sampleStartTimepoint = math.ceil(timepointsCurrentSample[0])
			sampleEndTimepoint = math.floor(timepointsCurrentSample[-1])
			sampleLength = sampleEndTimepoint - sampleStartTimepoint + 1.0
			timepointsCurrentSample = np.arange(sampleStartTimepoint, (sampleEndTimepoint + 1.0), SAMPLING_RATE)

			#Output unaligned samples data set file
			for currentTimepoint in timepointsCurrentSample:
				daysOfAntibiotics = int(round(interpolate.splev(currentTimepoint, currentSampleDaysOfAntibioticsSplineParameters)))
				gestationalAgeAtBirth = gestationalAgePerSubject[currentSampleID]
				postconceptionalAgeSampleObtained = (currentTimepoint/7.0) + gestationalAgeAtBirth
				gender = genderPerSubject[currentSampleID]
				modeOfBirth = modeOfBirthPerSubject[currentSampleID]
				periodOfStudy = periodOfStudyPerSubject[currentSampleID]
				roomType = roomTypePerSubject[currentSampleID]
				humanMilkUsed = humanMilkUsedPerSubject[currentSampleID]
				if daysOfAntibiotics < 0:
					daysOfAntibiotics = 0
				outline = currentSampleID + '\t' + str(currentTimepoint) + '\t' + str(gestationalAgeAtBirth) + '\t' + str(postconceptionalAgeSampleObtained) + '\t' + gender + '\t' + modeOfBirth + '\t' + roomType + '\t' + humanMilkUsed + '\t' + str(daysOfAntibiotics)
				#Save the taxon relative abundance values for this timepoint 
				taxonAbundanceValues = {}
				for taxaName, abundanceLevelPerTimepoint in currentSample.iteritems():
					abundanceValue = interpolate.splev(currentTimepoint, currentSampleTaxonSplineParameters[taxaName])
					abundanceValue = truncateAbundanceValue(abundanceValue)
					taxonAbundanceValues[taxaName] = abundanceValue
				#Output normalized relative abundance values for the selected taxa
				for taxaName in selectedTaxa:
					normalizedAbundanceValue = taxonAbundanceValues[taxaName] / sum(taxonAbundanceValues.values())
					outline += '\t' + str(normalizedAbundanceValue)
				outline += '\n'
				outfileUnaligned.writelines(outline)

			timepointsAlignedInterval = np.arange(math.ceil(alpha), (math.floor(beta) + 1.0), SAMPLING_RATE)

			#Output aligned samples data set file
			for referenceTimepoint in timepointsAlignedInterval:
				currentTimepoint = warpFunction(a, b, referenceTimepoint)
				daysOfAntibiotics = int(round(interpolate.splev(currentTimepoint, currentSampleDaysOfAntibioticsSplineParameters)))
				gestationalAgeAtBirth = gestationalAgePerSubject[currentSampleID]
				postconceptionalAgeSampleObtained = (currentTimepoint/7.0) + gestationalAgeAtBirth
				gender = genderPerSubject[currentSampleID]
				modeOfBirth = modeOfBirthPerSubject[currentSampleID]
				periodOfStudy = periodOfStudyPerSubject[currentSampleID]
				roomType = roomTypePerSubject[currentSampleID]
				humanMilkUsed = humanMilkUsedPerSubject[currentSampleID]
				if daysOfAntibiotics < 0:
					daysOfAntibiotics = 0
				outline = currentSampleID + '\t' + str(referenceTimepoint) + '\t' + str(gestationalAgeAtBirth) + '\t' + str(postconceptionalAgeSampleObtained) + '\t' + gender + '\t' + modeOfBirth + '\t' + roomType + '\t' + humanMilkUsed + '\t' + str(daysOfAntibiotics)
				#Save the taxon relative abundance values for this timepoint 
				taxonAbundanceValues = {}
				for taxaName, abundanceLevelPerTimepoint in currentSample.iteritems():
					abundanceValue = interpolate.splev(currentTimepoint, currentSampleTaxonSplineParameters[taxaName])
					abundanceValue = truncateAbundanceValue(abundanceValue)
					taxonAbundanceValues[taxaName] = abundanceValue
				#Output normalized relative abundance values for the selected taxa
				for taxaName in selectedTaxa:
					normalizedAbundanceValue = taxonAbundanceValues[taxaName] / sum(taxonAbundanceValues.values())
					outline += '\t' + str(normalizedAbundanceValue)
				outline += '\n'
				outfileAligned.writelines(outline)
	#Close output files
	outfileAligned.close()
	outfileUnaligned.close()

	return 

def getMean(sampleValues):
	if len(sampleValues) < 1:
		return 0.0

	return (sum(sampleValues) / float(len(sampleValues)))

def getVariance(sampleValues):
	variance = 0.0
	if len(sampleValues) < 1:
		return variance
	meanValue = sum(sampleValues) / float(len(sampleValues))
	for currentValue in sampleValues:
		variance += (currentValue - meanValue)**2
	variance = variance / float(len(sampleValues))

	return variance

def truncateAbundanceValue(sampleAbundanceValue):
	if sampleAbundanceValue < LOWER_BOUND:
		sampleAbundanceValue = LOWER_BOUND
	elif sampleAbundanceValue > UPPER_BOUND:
		sampleAbundanceValue = UPPER_BOUND
		
	return sampleAbundanceValue

def main(argv):
	if (len(argv) == 6):
		dataFilename = argv[1]
		taxaListFilename = argv[2]
		alignmentFilename = argv[3]
		outfilenameDBNnoAlignment = argv[4]
		outfilenameDBNalignment = argv[5]
	else:
		print "<<ERROR>> Invalid number of parameters!"
		return
	
	#Get taxa of interest from user-specified list
	selectedTaxa = getTaxaList(taxaListFilename)
	
	#Read dataset and prepare corresponding data structures
	getSamples(dataFilename, selectedTaxa, alignmentFilename, outfilenameDBNnoAlignment, outfilenameDBNalignment)

if __name__ == '__main__':
	main(sys.argv)
	#Alternate call, if you are using IDLE.
#	main(['generateCGBayesNetsDataset.py', './infant_gut/infant_gut_microbiota_annotated.tsv', './infant_gut/selectedTaxa.tsv', './infant_gut/infant_gut_microbiota_optimal_alignment.tsv', 'infant_gut_microbiota_noalignment', 'infant_gut_microbiota_alignment'])
