#! /usr/bin/python3

import sys
import datetime

# Variables:
inp	= sys.argv[1]
outF = sys.argv[2]
added = False
GQcut = 20
DPcut = 8
ABhetMin = 0.30
ABhetMax = 0.70
ABhomCut = 0.80
debug = False

# Find non mendelian inheratence
def nonMendelian(kid,mom,dad):
	flag = False
	# No uncalled
	if(kid[0] == "." or mom[0] == "." or dad[0] == "."):
		var = False	
	# Parents must have the same GT
	elif(mom[0:3] != dad[0:3]):
		var = False
	# The GT must be homozygous, only 1 possible allele to pass
	elif(mom[0] != mom[2]):
		var = False
	# Kid must not be homozygous
	elif(kid[0] == kid[2]):
		var = False
	# At least one allele must be same as parental
	elif(not((kid[0] == mom[0]) or (kid[2] == mom[0]))):
		var = False
		flag = True
	else:
		var = True
	return var,flag

# Check that the samples being used are the same family and in the correct role
def checkFam(kid,mom,dad):
	if(kid[0:4] != mom[0:4] or kid[0:4] != dad[0:4]):
		fam = False
	elif(kid[4:7] != "001"):
		fam = False
		print("Wrong ID ",kid)
	elif(mom[4:7] != "002"):
		fam = False
		print("Wrong ID ",mom)
	elif(dad[4:7] != "003"):
		fam = False
		print("Wrong ID ",dad)
	else:
		fam = True
	return fam

# Check the quality of each trio member's genotype call
def quality(var,fmtList,kid,mom,dad):
	passes = True
	if("GT" in fmtList):
		gtIdx = fmtList.index("GT")
	if("GQ" in fmtList):
		gqIdx = fmtList.index("GQ")
	else:
		gqIdx = ""
		print("Did not filter for GQ")
	if("DP" in fmtList): 
		dpIdx = fmtList.index("DP")
	else:
		dpIdx = ""
		print("Did not filter for DP")
	if("AD" in fmtList): 
		adIdx = fmtList.index("AD")
	else:
		adIdx = ""
		print("Did not filter for AD")
	famMems = [kid,mom,dad]
	if(len(var[0]) > 1 or len(var[1]) > 1 or var[0] == "*" or var[1] == "*"):
		passes = False
	if(passes):	
		for mem in famMems:
			if(gqIdx):
				GQ = mem.split(":")[gqIdx]
				if( not GQ.isnumeric()):
					passes = False
				elif(int(GQ) < GQcut):
					passes = False
			else:
				passes = False
				print("No GQ here, check out:\n",kid,mom,dad)
			if(dpIdx):
				DP = mem.split(":")[dpIdx]
				if( not DP.isnumeric()):
					passes = False
				elif(int(DP) < DPcut):
					passes = False
			else:
				passes = False
				print("No DP here, check out:\n",kid,mom,dad)
			if(passes and adIdx):
				AD = mem.split(":")[adIdx].split(",")
				GTal = mem.split(":")[gtIdx].split("/")
				adSum = 0
				for ads in AD:
					adSum += int(ads)
				al1 = int(GTal[0])
				al2 = int(GTal[1])
				al2dp = int(AD[al2])
				if(adSum == 0):
					print("\nFLAG02: AD SUM IS 0")
					print("adSum is 0, check out:")
					print("AD is: ",",".join(AD))
					print("al1,al2,al2dp,adSum = ",al1,al2,al2dp,adSum)
					print("Var1,Var2 = ",var[0],var[1]) 
					print("mem= ",mem,"\n")
					print(kid,mom,dad)
					ab = 0
					al1 = 1
					al2 = 2
				else:
					ab = round(al2dp/adSum,4)
				if(al1 != al2 and (ab < ABhetMin or ab > ABhetMax )):
					passes = False
				elif(al1 == al2 and ab < ABhomCut):
					print("\nFLAG03: HOM WITH AB < 1")
					print("Hom with AB less than 1, check out:\n")
					print("al1,al2,al2dp,adSum = ",al1,al2,al2dp,adSum)
					print("AB= ",ab)
					print("Var1,Var2 = ",var[0],var[1]) 
					print("mem= ",mem,"\n")
					print(kid,mom,dad)
	return passes

# Returns the allele balance when given the allele depth values and genotype
def getAB(inAD,inGT,altAlle):
	allDP = inAD.split(",")
	gtAl = int(inGT.split("/")[1])
	dpSum = sum([int(dp) for dp in allDP])
	if( dpSum != 0):
		outAB = round( (int(allDP[gtAl]) / dpSum) ,3)
	else:
		outAB = "Error"
	return(outAB)

# Return the QD for specific child sample
def getQD(inFmt,inKid):
	fmtVals = inFmt.split(":")
	kidVals = inKid.split(":")
	gtIx    = fmtVals.index("GT")
	plIx    = fmtVals.index("PL")
	dpIx    = fmtVals.index("DP")
	altAlle = int(kidVals[gtIx].split("/")[1])
	plVals  = kidVals[plIx].split(",")
	if(altAlle == 1):
		gtPL = sorted([int(plVals[0]),int(plVals[1]),int(plVals[2])])
	elif(altAlle == 2):
		gtPL = sorted([int(plVals[0]),int(plVals[3]),int(plVals[5])])
	elif(altAlle == 3):
		gtPL = sorted([int(plVals[0]),int(plVals[6]),int(plVals[9])])
	elif(altAlle == 4):
		gtPL = sorted([int(plVals[0]),int(plVals[10]),int(plVals[14])])
	elif(altAlle == 5):
		gtPL = sorted([int(plVals[0]),int(plVals[15]),int(plVals[20])])
	outQD = int(gtPL[1])/int(kidVals[dpIx])
	return(outQD)

# Remove values from INFO fields that are not from the alleles in this trio
def removeExtraInfo(inInfo, altGT):
	infVals = inInfo.split(";")
	newVals = []
	# Values in the info will be ordered allele first for 0/1, second for 0/2 etc.
	# Substracting 1 makes it compatible with 0 starting indexes
	altGtInt = int(altGT) - 1
	for inf in infVals:
		if( "," in inf):
			fld,vals = inf.split("=")
			valLst = vals.split(",")
			newInf = fld+"="+valLst[altGtInt]
			newVals.append(newInf)              
		else:
			newVals.append(inf)
	newInfoField = ";".join(newVals)
	return(newInfoField)

# Change the samples fields to have:
#  - AB value for both hets and homs
#  - GT values of 0 and 1, no more multiple alleles per line
#  - AD from different genotypes removed
#  - PL from different genotypes removed
def removeExtraFormat(inFormat,inKid,inMom,inDad):
	outAll = [[],[],[],[]]
	inAll = []
	flag = False
	for iStr in [inFormat,inKid,inMom,inDad]:
		inAll.append(iStr.split(":"))
	altAllele = int(inAll[1][0][2])
	inGTix = inAll[0].index("GT")
	inADix = inAll[0].index("AD") 
	for fIx in range(len(inAll[0])):
		if(fIx == 1):
			outAll[0].append("AB")
			for mIx in [1,2,3]:
				mAllele = inAll[mIx][inGTix]
				mA1 = int(mAllele[0]) 
				mA2 = int(mAllele[2])
				denom = int(inAll[mIx][inADix].split(",")[mA1]) + int(inAll[mIx][inADix].split(",")[mA2])
				if(denom == 0):
					print("Kid in filt: \n",inKid)
					print("Mom in filt: \n",inMom)
					print("Dad in filt: \n",inDad)
				newAB = getAB(inAll[mIx][inADix],mAllele,altAllele)
				if(newAB == "Error"):
					flag = True
					print("\nFLAG04: AB value error")
					print("Kid in filt: \n",inKid)
					print("Mom in filt: \n",inMom)
					print("Dad in filt: \n",inDad)
					newAB = inAll[mIx][inADix]+",ERR"
				outAll[mIx].append(str(newAB))
		if(inAll[0][fIx] == "GT"):
			outAll[0].append("GT")
			for mIx in [1,2,3]:
				mAlle1,mAlle2 = inAll[mIx][fIx].split("/")
				if(mAlle1 != "0"):
					mAlle1 = "1"
				if(mAlle2 != "0"):
					mAlle2 = "1"
				outAll[mIx].append(str(mAlle1+"/"+mAlle2))
		elif(inAll[0][fIx] == "PL"):
			outAll[0].append("PL")
			for mIx in [1,2,3]:
				plVals = inAll[mIx][fIx].split(",")
				if(str(altAllele) == "1"):
					newPL = plVals[0]+","+plVals[1]+","+plVals[2]	
				elif(str(altAllele) == "2"):
					newPL = plVals[0]+","+plVals[3]+","+plVals[5]	
				elif(str(altAllele) == "3"):
					newPL = plVals[0]+","+plVals[6]+","+plVals[9]	
				else:
					newPL = inAll[mIx][fIx]
				outAll[mIx].append(str(newPL))
		elif(inAll[0][fIx] == "AD"):
			outAll[0].append("AD")
			for mIx in [1,2,3]:
				adVals = inAll[mIx][fIx].split(",")
				newAD = adVals[0]+","+adVals[altAllele]
				outAll[mIx].append(str(newAD))
		elif(inAll[0][fIx] not in ["AD","GT","AB","PL"]):
			outAll[0].append(inAll[0][fIx])
			for mIx in [1,2,3]:
				outAll[mIx].append(inAll[mIx][fIx])
	outFmt = ":".join(outAll[0])
	outKid = ":".join(outAll[1]) 
	outMom = ":".join(outAll[2]) 
	outDad = ":".join(outAll[3])
	return(outFmt,outKid,outMom,outDad,flag)



# Main code:
with open(inp, "r") as iFh,open(outF,"w") as out:
	for iLine in iFh:
		if (iLine[0:2] == "##"):
			oLine = iLine.rstrip()
			if("##INFO=" in oLine and not added):
				oLine = "##INFO=<ID=TrioDP,Number=1,Type=Integer,Description=\"Total depth for the trio from trioID\">\n##INFO=<ID=TrioID,Number=1,Type=String,Description=\"Family ID of trio where SNV is found\">\n"+oLine
				added = True
			print(oLine,file=out)
			continue
		elif(iLine[0:2] == "#C"):
			oLine = iLine.rstrip()
			headEls = oLine.split("\t")
			newHeader = "\t".join(headEls[0:9]+["child","mother","father"])
			now = datetime.datetime.now()
			newHeader="##python_filter="+sys.argv[0]+" "+sys.argv[1]+" "+sys.argv[2]+"; Date="+str(now)+"\n"+newHeader
			print(newHeader,file=out)
			continue
		else:
			oLine = iLine.rstrip()
			varLine = oLine.split("\t")
			varList = [varLine[3]]
			# If there are more than 1 alternate allele, make them into a list
			# otherwise make a list with the single alt allele
			if("," in varLine[4]):
				varList += varLine[4].split(",")
			else:
				varList += [varLine[4]]
			# For all indices past 9 (FORMAT) process them in groups of 3
			# If the formatting is correct, these should be ordered child, mother and father
			for ind in range(9,len(headEls),3):
				kidID,momID,dadID = headEls[ind],headEls[ind+1],headEls[ind+2]
				kid,mom,dad = varLine[ind],varLine[ind+1],varLine[ind+2]
				# If the IDs are not compatible with being a trio, exit the script and ask for
				# verification of the sample IDs. If DEBUG mode is true, print the step the script is at.
				if(debug):
					print("Testing family IDs")
				if(not checkFam(kidID,momID,dadID)):
					print("The IDs for this set of three are:",kidID,momID,dadID,sep="\n")
					print("Their first 4 characters, where the Family ID is supposed to be, don't match")
					print("Please check header format and make sure there are 3 samples per trio in the correct order")
					sys.exit()
				# If the kid's genotype has a missing character, '.', then continue to
				# the next variant, otherwise fetch the correct variant from the variant list
				# Since variants are ordered starting at 0, the GT value should correspond to
				# the correct variant identity for the child
				if(kid[0] == "." or kid[2] == "."):
					continue
				else:
					varCheck = [varList[int(kid[0])],varList[int(kid[2])]]
				# If the reference or variant alleles are not a single nucleotide, continue to next
				# If any have the asterisk character representing a deletion, continue to next
				if(len(varCheck[0]) > 1 or len(varCheck[1]) > 1 or "*" in [varCheck]):
					continue
				# If debug is on, print the current IDs with the SNV ID and location being processed
				if(debug):
					print("\n\n\n")
					print("Processing "+kidID+", "+momID+", and "+dadID+" at position: "+varLine[0]+":"+varLine[1])
					print("kid:\n",kid)
					print("mom:\n",mom)
					print("dad:\n",dad)
					print("Testing for non mendelian inheritance")
				# Test if the proband (i.e. child,kid, sample001) has a genotype incompatible with the
				# parental genotypes. Flag variants with double DNMs (almost 99.99% are errors)
				nonMend,flag  = nonMendelian(kid,mom,dad) 
				# If the test returns a True flag, then report in case we want to double check the supposed
				# double DNM. (or incorrect parent files?)
				if(flag):
					print("\nFLAG01: DOUBLE DNM")
					print("Position: ",varLine[0],varLine[1])
					print("Ref: ",varLine[3])
					print("Alts: ",varLine[4])
					print(kidID,kid)
					print(momID,mom)
					print(dadID,dad)
				# If the variant passes as non mendelian inheritance, proceed to test quality.
				if(nonMend):
					# If DEBUG mode, print what step the script is at.
					if(debug):
						print("Passed famCheck and nonMend tests")
					# Remove any previously added TrioID, TrioDP and TrioQD values from the current  INFO field
					# For cases where there are multple samples in the file that have different alternate allele
					# that could be de novo SNVs
					cleanInfoFields = [(iV) for iV in varLine[7].split(";") if ("TrioID" not in iV and "TrioDP" not in iV and "TrioQD" not in iV)]
					# Add current trio ID to INFO field (varLine[7]) and append it to ID field (varLine[2]).
					# In single trio testing, the Family ID (or TrioID) doesn't need to be added to the ID field.
					# We modify it anyways for consistency.
					cleanInfoFields.append("TrioID="+kidID[0:4])
					varLine[7] = ";".join(cleanInfoFields)
					idField = varLine[2]+","+kidID[0:4]
					# Verify depth and quality, if those values pass the minimum requirement set at the start of 
					# the scrip, calculate trio specific depth and qualByDepth values.
					# Also remove values from other genotypes if present.
					varCheck = [varList[int(kid[0])],varList[int(kid[2])]]
					fmt = varLine[8].split(":")
					# Print current algorithm step if on DEBUG mode
					if(debug):
						print("Testing quality")
					if(quality(varCheck,fmt,kid,mom,dad)):
						# Print current algorithm step if on DEBUG mode
						if(debug):
							print("Passed quality test")
						# Find variant specific indeces for GT,DP, and PL values
						gtIdx = fmt.index("GT")
						dpIdx = fmt.index("DP")
						plIdx = fmt.index("PL")
						# Split kid sample info
						kidSampVals = kid.split(":")
						# Obtain depth value from each trio member
						cDp = int(kidSampVals[dpIdx])
						mDp = int(mom.split(":")[dpIdx])
						fDp = int(dad.split(":")[dpIdx])
						# Add it into a single TrioDepth value
						trioDP = cDp + mDp + fDp
						# Calculate trio specific Quality by Depth stat
						# Our heuristic is the childs GQ score divided by Depth.
						trioQD = round(getQD(varLine[8],kid),2)
						# Add the two new INFO fields: TrioDP and TrioQD
						cleanInfoFields.append("TrioDP="+str(trioDP))
						cleanInfoFields.append("TrioQD="+str(trioQD))
						# Use the removeExtraInfo function to extract the values that belong to other
						# alternate alleles that are not the current trios allele.
						# Only here for multi sample VCF files where all alternate alleles of the same
						# genomic coordinate are in one line.
						fixedInfo = removeExtraInfo(";".join(cleanInfoFields),kidSampVals[gtIdx][2])
						# Remove the value for other genotypes that are not in the trio. Same as above but
						# for the individual depth and phenotype likelihood scores. Also, adds AB values.
						fixedFmt, fixedKid, fixedMom, fixedDad, flag2 = removeExtraFormat( varLine[8], kid, mom, dad )
						# If the calculation of AB is impossible return an error flag to report.
						# This could happen if for some reason the Allele Depth for one of the trio member is 0 for both alleles.
						# Usually indicates the quality filters are too low.
						if(flag2):
							print("Position: ",varLine[0],varLine[1])
							print("FMT flag")
						# If the kid's genotype doesn't have the reference allele (e.g. kid genotype 1/2 or 2/3)
						# then but both alternate alleles in the REF field of the VCF. 
						# Otherwise make the ALT field the alternate allele in the kid genotype.
						# Example 1: If reference is A and kid is G/T then REF field is A and ALT field is G,T
						# Example 2: If reference is A and kid is A/T then REF field is A and ALT field is T
						if(varLine[3] in varCheck):
							subVars = varCheck[1]
						else:
							subVars = varCheck[0]+","+varCheck[1]
						if(debug):
							print("PASSED ALL")
						print("\t".join(varLine[0:2]+[idField,varLine[3],subVars]+varLine[5:7]+[fixedInfo,fixedFmt,fixedKid,fixedMom,fixedDad]),file=out)
					else:
						if(debug):
							print("QUALITY FAILED")
				else:
					if(debug):
						print("This variant follows mendelian inheritance")
