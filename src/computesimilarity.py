from __future__ import division



def substitute_annotation(profile,replaced,annotationpool):
	
	available=set.difference(set(range(0,len(profile))),replaced)
	replaceindex=random.sample(available,1)[0]
	replaced.add(replaceindex)
	replacementindex=random.randint(0,len(annotationpool)-1)
	profile[replaceindex]=annotationpool[replacementindex]
	return profile,replaced


def bp_test(profile1,profile2,icdict,ancestors):
	finalsim=0
	bestmatchic=[]
	termmatchic=[]
	matchdata=[]
	for term1 in profile1:
		termmatchic=[]
		for term2 in profile2:
			micaic,mica=getmicaic(term1,term2,ancestors,icdict)
			print term1,"\t",term2,"\t",mica,"\t",micaic
			
			
		





def calculate_bestpairs(profile1,profile2,icdict,ancestors,similaritydict):
	finalsim=0
	bestmatchic=[]
	termmatchic=[]
	matchdata=[]
	for term1 in profile1:
		termmatchic=[]
		for term2 in profile2:
			termtuple=tuple(sorted((term1,term2)))
			if termtuple in similaritydict:
				termmatchic.append(similaritydict[termtuple])
			
			else:
				micaic,mica=getmicaic(term1,term2,ancestors,icdict)
				termmatchic.append((term2,micaic,mica))
				similaritydict[termtuple]=(term2,micaic,mica)
		bestmatch,bestic,bestmica=getmax(termmatchic)
		bestmatchic.append(bestic)
		matchdata.append((term1,bestmatch,bestmica,bestic))	
	mediansim=np.median(bestmatchic)
	return mediansim,matchdata,similaritydict,bestmatchic

def getmax(lcslist):
	micaic=0
	mica=""
	match=""
	for tup in lcslist:
		if tup[1]>micaic:
			micaic=tup[1]
			mica=tup[2]
			match=tup[0]
	return match,micaic,mica
			
def getmicaic(term1,term2,ancestors,icdict):
	micaic=0
	mica=""
	commonancestors=set.intersection(ancestors[term1],ancestors[term2])
	lcslist=[(term2,icdict[anc],anc) for anc in commonancestors]
	match,micaic,mica=getmax(lcslist)

	if len(lcslist)>0:
		return micaic,mica
	else:
		return 0,"None"

def compute(randomprofiles,annotationpool,ancestors,icdict,queryprofilesize,numqueryprofiles):
	selfmatch=dict()
	out=open("../results/PilotResults_ExAttr(Q)_OldProfiles_Test.tsv",'w')
	out.write("Query profile ID\tNumber of annotations replaced\tBest taxon match\tMedian similarity\tIC list\t\n")
	queryprofileids=set()
	bestmatchsimilarity=dict()
	similaritydict=dict()
	for profileid in randomprofiles:
		if len(randomprofiles[profileid])==queryprofilesize:
			queryprofileids.add(profileid)
	randomqueryprofileids=random.sample(queryprofileids, numqueryprofiles)
	for profileid in randomqueryprofileids:
		queryprofile=deepcopy(randomprofiles[profileid])
		queryprofilesize=len(queryprofile)
		print "Conducting experiment on ",profileid," Size",queryprofilesize
		print "Profile",randomprofiles[profileid]
		

		replaced=set()
		for i in range(0,queryprofilesize+1):
			results=[]
			for databaseprofileid in randomprofiles:
				bpmediansim=0
				dbprofile=randomprofiles[databaseprofileid]

				bpmediansim,matchdata,similaritydict,bestmatchic=calculate_bestpairs(queryprofile,dbprofile,icdict,ancestors,similaritydict)
				resulttuple=(bpmediansim,matchdata,bestmatchic,databaseprofileid)
				results.append(resulttuple)
				
			
			print "Number replaced, ", i
			bestmatch=max(results,key=itemgetter(0))[3]
			similarity= max(results,key=itemgetter(0))[0]
			iclist=",".join(str(x) for x in max(results,key=itemgetter(0))[2])
			bestmatchdata=max(results,key=itemgetter(0))[1]
			# if i==0:
			# 	for x in bestmatchdata:
			# 		#term1,bestmatch,bestmica,bestic
			# 		if x[0] == x[1] == x[2]:
			# 			selfmatch[x[0]] = x[3]
			# else:
			# 	for x in bestmatchdata:
			# 		print x[0],"\t\t",x[1],"\t\t",x[2]
			# 		if x[0] != x[1]:
			# 			print "inside"
			# 			if x[3] > selfmatch[x[1]]:
			# 				print x[0],"\t","best match ","\t",x[1],"\t","MICA\t",x[2],x[3]
			# 				if x[2] not in ancestors[x[1]] or x[2] not in ancestors[x[2]]:
			# 					print "MICA not in ancestors"


			out.write(profileid+"\t"+str(i)+"\t"+bestmatch+"\t"+str(similarity)+"\t"+iclist+"\n")
			if len(replaced)<len(queryprofile):
				queryprofile,replaced=substitute_annotation(queryprofile,replaced,annotationpool)
	out.close()

def check_ic(annotationpool,ancestors,icdict,frequency):
	for annotation in set(annotationpool):
		annic=icdict[annotation]
		for anc in ancestors[annotation]:
			ancic=icdict[anc]
			if ancic > annic:
				print annotation,annic,anc,ancic
				print "Ann freq",frequency[annotation],"Anc freq",frequency[anc]

def compute_ic(profiles,ancestors):
	annotationcorpus=open("../AnnotationCorpus.txt",'w')
	annotationlist=[]
	icdict=dict()
	frequency=dict()
	annotationandanclist=dict()
	for profileid in profiles:
		for annotation in profiles[profileid]:
			annotationcorpus.write(annotation+"\n")
			annotationlist.append(annotation)
			
			for anc in ancestors[annotation]:
				annotationcorpus.write(anc+"\n")
				if anc not in frequency:
					frequency[anc]=0
				frequency[anc]+=1


	corpussize=len(annotationlist)
	print "Corpus Size",corpussize
	maxic=round(-math.log(1/corpussize),2)
	for term in frequency:
		freq=frequency[term]
		ic=round((-math.log(freq/corpussize))/maxic,2)
		icdict[term]=ic
	annotationcorpus.close()
	
	return icdict,frequency

def load_randomprofiles():
	randomprofiles=dict()
	annotationpool=[]
	infile=open("../data/RandomProfilesOld.txt")
	for line in infile:
		profileid,annotation=line.strip().split("\t")
		if profileid not in randomprofiles:
			randomprofiles[profileid]=[]
		randomprofiles[profileid].append(annotation)
		annotationpool.append(annotation)
	infile.close()
	return randomprofiles,annotationpool

def load_ancestors():
	ancestors=dict()
	infile=open("../data/Subsumers_EAttr(Q)_OldRealProfiles.txt")
	for line in infile:
		term,subsumer=line.replace(">","").replace("<","").strip().split("\t")
		if term not in ancestors:
			ancestors[term]=set()
			ancestors[term].add(term)
		if subsumer !="owl:Thing":		
			ancestors[term].add(subsumer)
	infile.close()
	return ancestors

def main():
	queryprofilesize=10
	numqueryprofiles=5
	randomprofiles,annotationpool=load_randomprofiles()
	print "Loaded randomprofiles"
	ancestors=load_ancestors() 
	print "Loaded ancestors"
	icdict,frequency=compute_ic(randomprofiles,ancestors)
	print "Computed IC"
	#check_ic(annotationpool,ancestors,icdict,frequency)
	similaritydict=dict()
	compute(randomprofiles,annotationpool,ancestors,icdict,queryprofilesize,numqueryprofiles)
	#bp_test(randomprofiles['VTO_0035450'],randomprofiles['VTO_0035450'],icdict,ancestors)


	

if __name__ == "__main__":
	import random
	from copy import deepcopy
	from operator import itemgetter
	import numpy as np
	import math
	main()
