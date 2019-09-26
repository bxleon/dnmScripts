#! /usr/bin/python3

import sys

inpF = sys.argv[1]
outF = sys.argv[2]

with open(inpF, 'r') as iFh, open(outF,'w') as out:
	for rLine in iFh:
		line = rLine.rstrip()
		if(line[0] == "#"):
			print(line,file = out)
		else:
			# Parse the VCF line
			cols = line.split("\t")
			fmtVals = cols[8].split(":")
			idxGT = fmtVals.index("GT")
			idxPL = fmtVals.index("PL")
			samples = []
			for ix in range(9,12):
				samples.append(cols[ix].split(":")) 
			gtsAll = []
			
			# For all the samples: remove GT alleles, add to a list with all the allele for the samples
			# and check if it's hom ref, hom alt or het to change to 0/1 values.
			for smp in samples:
				gtsVals = smp[idxGT].split("/") 
				gtsAll += gtsVals 
				if(gtsVals[0] == gtsVals[1] and gtsVals[0] == "0"):
					if(gtsVals[0] == "0"):
						smp[idxGT] = "0/0"
					else:
						smp[idxGT] = "1/1"
				elif(gtsVals[0] != gtsVals[1] and gtsVals[0] == "0"):
					smp[idxGT] = "0/1"
					ogHet = int(gtsVals[1])-1
				else:
					print("The genotype for position :",cols[0],cols[1]," don't have a reference allele.\nPlease exclude or verify.")
			
			# Once a list of all alleles for the samples is generated, iterate over samples again 
			# and change the PL scores to match the three possibles PLs for the correct allele.
			# If there is no reference allele along with a 1 or 2 allele value in the list, send an error message 
			for smp in samples:
				allPls = smp[idxPL].split(",")
				if("1" in gtsAll and "0" in gtsAll):
					newPls = ",".join([allPls[0],allPls[1],allPls[2]])
				elif("2" in gtsAll and "0" in gtsAll):
					newPls = ",".join([allPls[0],allPls[3],allPls[5]])
				elif("3" in gtsAll and "0" in gtsAll):
					newPls = ",".join([allPls[0],allPls[6],allPls[9]])
				elif("0" not in gtsAll):
					print("The genotype for position :",cols[0],cols[1]," don't have a reference allele.\nPlease exclude or verify.")
				smp[idxPL] = newPls

			# Change the info filds according to the original genotype of the alternate
			infVals = cols[7].split(";")
			newVals = []
			for inf in infVals:
				if( "," in inf):
					fld,vals = inf.split("=")
					valLst = vals.split(",")
					newInf = fld+"="+valLst[ogHet]
					newVals.append(newInf)				
				else:
					newVals.append(inf)
			newInfoField = ";".join(newVals)

			# Make the new line to output by adding all columns from CHROM to FORMAT and then add the modded samples
			newLine = "\t".join(cols[0:7])
			newLine += "\t" + "\t".join([newInfoField,cols[8]])
			for smp in samples:
				newLine += "\t"+ ":".join(smp)
			print(newLine,file=out)
