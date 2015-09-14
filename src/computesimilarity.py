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
			
			
		

def create_random_profiles(size,annotationpool,numprofiles):
	queryprofiles=[]
	randomindices=[]
	randomprofileindices=dict()
	selected=set()		
	poolsize=len(annotationpool)
	realindices=list(range(0,len(annotationpool)))
	while (len(queryprofiles)<numprofiles):
		profilesize=0
		selected=set()
		tempprofiles=[]
		while profilesize<size:
			randomindex=random.randint(0,poolsize-1)
			if (randomindex not in selected):
				tempprofiles.append(annotationpool[randomindex])
				selected.add(randomindex)
				profilesize+=1
		queryprofiles.append(tempprofiles)
	return queryprofiles	




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


def compute_without_substitution(queryprofiles,dbprofiles,ancestors,icdict,similaritydict):
	similarityscores=[]
	querysizes=[]
	dbprofilesizes=[]
	allsimilarityscores=[]
	scores=[]
	for query in queryprofiles:
		results=[]
		for databaseprofileid in dbprofiles:
			querysizes.append(len(query))
			dbprofilesizes.append(len(dbprofiles[databaseprofileid]))
			bpmediansim=0
			dbprofile=dbprofiles[databaseprofileid]
			bpmediansim,matchdata,similaritydict,bestmatchic=calculate_bestpairs(query,dbprofile,icdict,ancestors,similaritydict)
			resulttuple=(bpmediansim,matchdata,bestmatchic,databaseprofileid)
			results.append(resulttuple)
			allsimilarityscores.append(bpmediansim)
		bestmatch=max(results,key=itemgetter(0))[3]
		similarity= max(results,key=itemgetter(0))[0]
		iclist=",".join(str(x) for x in max(results,key=itemgetter(0))[2])
		bestmatchdata=max(results,key=itemgetter(0))[1]
		similarityscores.append(similarity)
	return similarityscores,allsimilarityscores,querysizes,dbprofilesizes,similaritydict




def compute_with_substitution(randomprofiles,annotationpool,ancestors,icdict,queryprofilesize,numqueryprofiles,similaritydict):
	selfmatch=dict()
	bestsimilaritylist=[]
	quantiles=dict()
	queryprofileids=set()
	bestmatchsimilarity=dict()
	for profileid in randomprofiles:
		if len(randomprofiles[profileid])==queryprofilesize:
			queryprofileids.add(profileid)
	randomqueryprofileids=random.sample(queryprofileids, numqueryprofiles)
	for profileid in randomqueryprofileids:
		queryprofile=deepcopy(randomprofiles[profileid])
		queryprofilesize=len(queryprofile)
		replaced=set()
		for i in range(0,queryprofilesize+1):
			
			if i not in quantiles:
				quantiles[i]=dict()
			results=[]
			for databaseprofileid in randomprofiles:
				bpmediansim=0
				dbprofile=randomprofiles[databaseprofileid]

				bpmediansim,matchdata,similaritydict,bestmatchic=calculate_bestpairs(queryprofile,dbprofile,icdict,ancestors,similaritydict)
				resulttuple=(bpmediansim,matchdata,bestmatchic,databaseprofileid)
				results.append(resulttuple)
				
			allsimilarityscores=[]
			for res in results:
				allsimilarityscores.append(res[0])
			for j in [99,99.1,99.2,99.3,99.4,99.5,99.6,99.7,99.8,99.9,100]:
				if j not in quantiles[i]:
					quantiles[i][j]=[]
				quantiles[i][j].append(np.percentile(allsimilarityscores,j))	
			bestmatch=max(results,key=itemgetter(0))[3]
			similarity= max(results,key=itemgetter(0))[0]
			iclist=",".join(str(x) for x in max(results,key=itemgetter(0))[2])
			bestmatchdata=max(results,key=itemgetter(0))[1]
			bestsimilaritylist.append([profileid,i,bestmatch,similarity,iclist])
		

			if len(replaced)<len(queryprofile):
				queryprofile,replaced=substitute_annotation(queryprofile,replaced,annotationpool)

	
	return similaritydict,bestsimilaritylist,quantiles

def plot(figuredir):
	fig = plt.figure()
	for filename in os.listdir(figuredir):
		infile=open(figuredir+filename)
		scoredict=dict()
		expectdict=dict()
		numreplacedset=set()
		for line in infile:
			if "Number" not in line:
				data=line.strip().split("\t")
				profileid,numreplaced,sim,expect=data[0],int(data[1]),float(data[3]),float(data[4])
				numreplacedset.add(numreplaced)
				if numreplaced not in scoredict:
					scoredict[numreplaced]=[]
				scoredict[numreplaced].append(sim)
				if numreplaced not in expectdict:
					expectdict[numreplaced]=[]
				expectdict[numreplaced].append(sim)
		meansim=[]
		for numreplaced in scoredict:
			meansim.append(np.mean(scoredict[numreplaced]))
			meanexpect.append(np.mean(expectdict[numreplaced]))

		plt.plot(list(numreplacedset),meansim)

	plt.xlabel('Number of annotations replaced')
	plt.ylabel('Median best pair similarity')
	fig.savefig("../results/DecaySimilarity.png", dpi=1200)


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
	#print "Corpus Size",corpussize
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

def regression(similarityscores, queryprofilesizes,dbprofilesizes):
	box_query, querylambda_ = boxcox(np.array(queryprofilesizes) + 1)
	box_db, dblambda_ = boxcox(np.array(dbprofilesizes) + 1)
	sizes = [box_query,box_db]

	ones = np.ones(len(sizes[0]))
	X = sm.add_constant(np.column_stack((sizes[0], ones)))
	for ele in sizes[1:]:
		X = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(similarityscores, X).fit()
	


	predicted = results.fittedvalues
	
	return results,box_query,box_db,querylambda_,dblambda_

def prepare_expect_inputs(actualsimilarityscores,box_query,box_db,results):
	dbcoeff=results.params[0]
	querycoeff=results.params[1]
	constant=results.params[2]


	residuals=[]
	studentizedresiduals=[]
	expectscores=[]
	

	for i in range(0,len(actualsimilarityscores)):
		predictedscore=constant+(querycoeff*box_query[i])+(dbcoeff*box_db[i])
		residuals.append(actualsimilarityscores[i]-predictedscore)

	stddev=np.std(residuals)
	return stddev

def compute_expect_scores(results,actualsimilarity,boxquerysize,boxdbsize,resstddev,numofdbprofiles):
	# now compute expect scores for the individual comparison

	dbcoeff=results.params[0]
	querycoeff=results.params[1]
	constant=results.params[2]




	predictedscore=constant+(querycoeff*boxquerysize)+(dbcoeff*boxdbsize)
	residual=actualsimilarity-predictedscore
	stdres=residual/resstddev
	pvalue=1-math.exp(-math.exp(-stdres*math.pi/math.sqrt(6)+ 0.5772156649))
	expect=pvalue*numofdbprofiles
	return expect

def load():
	infile=open("../results/Scores.tsv")
	allsimilarityscores=[]
	querysizes=[]
	dbprofilesizes=[]
	for line in infile:
		score,querysize,dbsize=line.strip().split("\t")
		allsimilarityscores.append(float(score))
		querysizes.append(float(querysize))
		dbprofilesizes.append(float(dbsize))
	return allsimilarityscores,querysizes,dbprofilesizes

def quant_results(infile,size):
	scores=dict()
	initialscores=[]
	finalscores=[]
	decay=[]
	numreplacedlist=[]
	scorelist=[]
	for line in infile:
		if "Number" not in line:
			data=line.strip().split("\t")
			profileid,numreplaced,sim=data[0],int(data[1]),float(data[3])
			scorelist.append(sim)
			numreplacedlist.append(numreplaced)
			if numreplaced==0:
				initialscores.append(sim)
			if numreplaced ==size:
				finalscores.append(sim)
			if profileid not in scores:
				scores[profileid]=dict()
			scores[profileid][numreplaced]=sim
	infile.close()
	
	for i in range (0,len(initialscores)):
		decay.append(((initialscores[i]-finalscores[i])/finalscores[i])*100)
	print "Mean decay",np.mean(decay)
	print "Mean decayed similarity",np.mean(finalscores)
	fig = plt.figure()
	for profileid in scores:
		numreplacedlist=[]
		scorelist=[]
		for replaced in scores[profileid]:
			numreplacedlist.append(replaced)
			scorelist.append(scores[profileid][replaced])
		plt.plot
	plt.xlabel('Number of annotations replaced')
	plt.ylabel('Median best pair similarity')
	plt.legend(['Query 1', 'Query 2','Query 3', 'Query 4','Query 5'], loc='upper right')
	fig.suptitle("Query Profile Size "+str(size))
	fig.savefig("../results/Decay_Profilesize_"+str(size)+".png", dpi=1200)






def get_summary():
	queryprofilesize=40
	infile=open("../results/ProfileSize40_Results.tsv")
	quant_results(infile,queryprofilesize)

	queryprofilesize=20
	infile=open("../results/ProfileSize20_Results.tsv")
	quant_results(infile,queryprofilesize)

	queryprofilesize=10
	infile=open("../results/ProfileSize10_Results.tsv")
	quant_results(infile,queryprofilesize)




def plot_decay():
	figuredir="../results/Decay/"
	fig = plt.figure()
	fig1 = plt.figure()
	legend=[]
	for filename in os.listdir(figuredir):
		infile=open(figuredir+filename)
		scoredict=dict()
		expectdict=dict()
		numreplacedset=set()
		meanexpect=[]
		meansim=[]
		error=[]
		for line in infile:
			if "Number" not in line:
				data=line.strip().split("\t")
				profileid,numreplaced,sim,expect=data[0],int(data[1]),float(data[3]),float(data[4])
				numreplacedset.add(numreplaced)
				if numreplaced not in scoredict:
					scoredict[numreplaced]=[]
				scoredict[numreplaced].append(sim)
				if numreplaced not in expectdict:
					expectdict[numreplaced]=[]
				expectdict[numreplaced].append(sim)
		legend.append("profile size: "+str(np.max(numreplaced)))		
		for numreplaced in scoredict:
			meansim.append(np.mean(scoredict[numreplaced]))
			meanexpect.append(np.mean(expectdict[numreplaced]))
			error.append(2*(np.std(scoredict[numreplaced])/math.sqrt(len(scoredict[numreplaced]))))
		plt.errorbar(list(numreplacedset),meansim,yerr=error)

	plt.ylim(0,1)
	plt.legend(legend, loc='upper right')
	plt.xlabel('Number of annotations replaced')
	plt.ylabel('Median best pair similarity')
	fig.savefig("../results/DecaySimilarity.png", dpi=1200)

def plot_expect():
	fig= plt.figure()
	x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	infile=open("../results/ProfileSize10_Results.tsv")
	scoredict=dict()
	expectdict=dict()
	for line in infile:
		if "Query" not in line:
			data=line.strip().split("\t")
			query,numreplaced,match,sim,expect=data[0],int(data[1]),data[2],float(data[3]),float(data[4])
			if numreplaced not in scoredict:
				scoredict[numreplaced]=[]
			scoredict[numreplaced].append(sim)

			if numreplaced not in expectdict:
				expectdict[numreplaced]=[]
			expectdict[numreplaced].append(expect)

	expectlist=[]
	scorelist=[]

	for numreplaced in scoredict:
		print numreplaced
		expectlist.append(math.log(np.mean(expectdict[numreplaced])))
		scorelist.append(np.mean(scoredict[numreplaced]))
	plt.scatter(x,scorelist)
	plt.plot(x,expectlist)
	plt.xlim(0,10.5)
	plt.legend(['log(Expect score)','Median IC'], loc='lower right')
	plt.xlabel('Num of annotations replaced')
	plt.ylabel('Score')
	fig.savefig("../results/Similarity_Expect.png", dpi=1200)




	

def plotquantiles(noisequantiles,quantiles):
	fig, ax = plt.subplots()


	numreplaced=[]
	x=list(range(0,12))
	for i in [99,99.5,100]:
		quantilelist=[]
		for j in range(0,11):
			quantilelist.append(np.mean(quantiles[j][i]))
		quantilelist.append(noisequantiles[i])
		plt.plot(x,quantilelist)
		print x
		print quantilelist
	legend=['99th percentile','99.5th percentile','100th percentile']		
	# quantilex=[99,99.1,99.2,99.3,99.4,99.5,99.6,99.7,99.8,99.9,100]
	# numreplacedlist=[]
	# legend=[]
	# legend.append('noise')
	# plt.plot(quantilex,noisequantiles,'--')
	# print "Noise ",noisequantiles
	# for numreplaced in quantiles:
	# 	quantilelist=[]
	# 	legend.append(numreplaced)
	# 	for j in quantilex:
	# 		quantilelist.append(np.mean(quantiles[numreplaced][j]))
	# 	plt.plot(quantilex,quantilelist)
	# 	print numreplaced,quantilelist
	#plt.xlim(99,100)
	#plt.ylim(0,1)
	# plt.ylabel('Similarity score')
	# plt.legend(legend, loc='lower left',ncol=3)
	# plt.xlabel('Quantiles')
	# plt.show()
	plt.ylim(0,1)
	ax.set_xticklabels(['0','1','2','3','4','5','6','7','8','9','10','Noise'])
	plt.ylabel('Score')
	plt.legend(legend, loc='lower left',ncol=2)
	plt.xlabel('Num of annotations replaced')
	plt.show()	

def do(queryprofilesize,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,similaritydict,numofdbprofiles):
	filename="../results/Decay/ProfileSize"+str(queryprofilesize)+"_Results.tsv"	
	out=open(filename,'w')
	print "Profile size",queryprofilesize
	queryprofiles=create_random_profiles(queryprofilesize,annotationpool,5)
	
	similaritydict,bestsimilaritylist,quantiles=compute_with_substitution(dbprofiles,annotationpool,ancestors,icdict,queryprofilesize,numqueryprofiles,similaritydict)

	similarityscores,allsimilarityscores,querysizes,dbprofilesizes,similaritydict=compute_without_substitution(queryprofiles,dbprofiles,ancestors,icdict,similaritydict)
	
	noisequantiles=dict()
	for j in [99,99.1,99.2,99.3,99.4,99.5,99.6,99.7,99.8,99.9,100]:
		noisequantiles[j]=np.percentile(allsimilarityscores,j)

	#plotquantiles(noisequantiles,quantiles)

	print "Mean noise to noise similarity",np.mean(similarityscores)

	results,box_query,box_db,querylambda_,dblambda_=regression(allsimilarityscores,querysizes,dbprofilesizes)
	resstddev=np.std(results.resid)

	#resstddev=prepare_expect_inputs(allsimilarityscores,box_query,box_db,results)

	out.write("Query profile ID\tNumber of annotations replaced\tBest taxon match\tMedian similarity\tExpect score\n")
	for tup in bestsimilaritylist:
		query=tup[0]
		numreplaced=tup[1]
		bestmatch=tup[2]
		similarity=tup[3]
		iclist=tup[4]
		dbprofilesize=len(dbprofiles[bestmatch])
		boxquerysize = boxcox(np.array([queryprofilesize]) + 1,lmbda=querylambda_)[0]
		boxdbsize = boxcox(np.array([dbprofilesize]) + 1,lmbda=dblambda_)[0]
		expectscore=compute_expect_scores(results,similarity,boxquerysize,boxdbsize,resstddev,numofdbprofiles)
		out.write(query+"\t"+str(numreplaced)+"\t"+bestmatch+"\t"+str(similarity)+"\t"+str(expectscore)+"\n")
	out.close()
	return similaritydict


def main():
	plot_decay()
	sys.exit()
	queryprofilesize=sys.argv[1]
	numqueryprofiles=int(sys.argv[2])
	numofdbprofiles=659

	similaritydict=dict()
	dbprofiles,annotationpool=load_randomprofiles()
	print "Loaded randomprofiles"
	ancestors=load_ancestors() 
	print "Loaded ancestors"
	icdict,frequency=compute_ic(dbprofiles,ancestors)
	print "Computed IC"
	
	similaritydict=do(10,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,similaritydict,numofdbprofiles)
	print "Done with 10"
	similaritydict=do(20,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,similaritydict,numofdbprofiles)
	print "Done with 20"
	similaritydict=do(40,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,similaritydict,numofdbprofiles)
	print "Done with 40"

	



if __name__ == "__main__":
	import random
	import matplotlib.pyplot as plt
	from copy import deepcopy
	from operator import itemgetter
	import numpy as np
	import math
	import statsmodels.api as sm
	import statsmodels.stats.api as sms
	from scipy import stats
	from scipy.stats import boxcox
	from statsmodels.stats import stattools
	import sys
	main()
