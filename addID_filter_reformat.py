#! /usr/bin/python3

import sys
import datetime

# Variables:
inp	= sys.argv[1]
outF = sys.argv[2]
added = False
GQcut = 20
DPcut = 8
ABhetMin = 0.25
ABhetMax = 0.75
ABhomCut = 0.80


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
		print("Double DNM? Check:\n",kid,mom,dad)
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
				adSum = 0
				for ads in AD:
					adSum += int(ads)
				al1 = int(mem[0])
				al2 = int(mem[2])
				al2dp = int(AD[al2])
				if(adSum == 0):
					print("FLAG02: AD SUM IS 0")
					print("adSum is 0, check out:\n")
					print("al1,al2,al2dp,adSum = ",al1,al2,al2dp,adSum)
					print("AB= ",ab)
					print("Var1,Var2 = ",var[0],var[1]) 
					print("mem= ",mem,"\n")
					print(kid,mom,dad)
					ab = 0
					al1 = 1
					al2 = 2
				else:
					ab = round(al2dp/adSum,4)
				if(al1 != al2 and (ab < ABhetMin or ABhetMax >0.75)):
					passes = False
				elif(al1 == al2 and ab < ABhomCut):
					print("FLAG03: HOM WITH AB < 1")
					print("Hom with AB less than 1, check out:\n")
					print("al1,al2,al2dp,adSum = ",al1,al2,al2dp,adSum)
					print("AB= ",ab)
					print("Var1,Var2 = ",var[0],var[1]) 
					print("mem= ",mem,"\n")
					print(kid,mom,dad)
	return passes


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
			if("," in varLine[4]):
				varList += varLine[4].split(",")
			else:
				varList += [varLine[4]]
			for ind in range(9,len(headEls),3):
				kidID,momID,dadID = headEls[ind],headEls[ind+1],headEls[ind+2]
				kid,mom,dad = varLine[ind],varLine[ind+1],varLine[ind+2]
				#Check they are correct IDs
				nonMend,flag  = nonMendelian(kid,mom,dad) 
				if(flag):
					print("FLAG01: DOUBLE DNM")
					print(iLine)
				if(checkFam(kidID,momID,dadID) and nonMend):
					# Add trio ID
					varLine[7] += ";TrioID="+kidID[0:4]
					idField = varLine[2]+","+kidID[0:4]
					# Verify depth and quality, if pass calculate trio specific depth and print
					fmt = varLine[8].split(":")
					varCheck = [varList[int(kid[0])],varList[int(kid[2])]]
					if(quality(varCheck,fmt,kid,mom,dad)):
						dpIdx = fmt.index("DP")
						cDp = int(kid.split(":")[dpIdx])
						mDp = int(mom.split(":")[dpIdx])
						fDp = int(dad.split(":")[dpIdx])
						trioDP = cDp + mDp + fDp
						varLine[7] += ";TrioDP="+str(trioDP)
						if(varLine[3] in varCheck):
							subVars = varCheck[1]
						else:
							subVars = varCheck[0],varCheck[1]
						print("\t".join(varLine[0:2]+[idField,varLine[3],subVars]+varLine[5:9]+[kid,mom,dad]),file=out)    


