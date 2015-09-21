from __future__ import division



def substitute_annotation(profile,replaced,annotationpool):
	
	available=set.difference(set(range(0,len(profile))),replaced)
	replaceindex=random.sample(available,1)[0]
	replaced.add(replaceindex)
	replacementindex=random.randint(0,len(annotationpool)-1)
	profile[replaceindex]=annotationpool[replacementindex]
	return profile,replaced

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


def get_allpairs_medianic(profile1,profile2,icdict,ancestors,bpsimilaritydict,apsimilaritydict):
	medianic=0
	termmatchic=[]
	for term1 in profile1:
		for term2 in profile2:
			termtuple=tuple(sorted((term1,term2)))

			if termtuple in bpsimilaritydict:
				termmatchic.append(bpsimilaritydict[termtuple][1])
				apsimilaritydict.append((bpsimilaritydict[termtuple][1],bpsimilaritydict[termtuple][2]))
			elif termtuple in apsimilaritydict:
				termmatchic.append(apsimilaritydict[termtuple][0])
			else:
				micaic,mica=getmicaic(term1,term2,ancestors,icdict)
				termmatchic.append(micaic)
				apsimilaritydict[termtuple]=(micaic,mica)	
	return np.median(termmatchic),apsimilaritydict



def calculate_bestpairs_symmetric(profile1,profile2,icdict,ancestors,bpsimilaritydict):
	finalsim=0
	bestmatchiclist=[]
	termmatchic=[]
	matchdata=[]
	for term1 in profile1:
		termmatchic=[]
		for term2 in profile2:
			termtuple=tuple(sorted((term1,term2)))
			if termtuple in bpsimilaritydict:
				termmatchic.append(bpsimilaritydict[termtuple])
			
			else:
				micaic,mica=getmicaic(term1,term2,ancestors,icdict)
				termmatchic.append((term2,micaic,mica))
				bpsimilaritydict[termtuple]=(term2,micaic,mica)
		bestmatch,bestic,bestmica=getmax(termmatchic)
	mediansim1=np.median(bestmatchiclist)
	

	finalsim=0
	bestmatchiclist=[]
	termmatchic=[]
	matchdata=[]
	for term1 in profile2:
		termmatchic=[]
		for term2 in profile1:
			termtuple=tuple(sorted((term1,term2)))
			termmatchic.append(bpsimilaritydict[termtuple])
		bestmatch,bestic,bestmica=getmax(termmatchic)
		bestmatchiclist.append(bestic)
			
	mediansim2=np.median(bestmatchiclist)
	return np.mean(mediansim1,mediansim2),bpsimilaritydict


def calculate_bestpairs_asymmetric(profile1,profile2,icdict,ancestors,bpsimilaritydict):
	finalsim=0
	bestmatchiclist=[]
	termmatchic=[]
	matchdata=[]
	for term1 in profile1:
		termmatchic=[]
		for term2 in profile2:
			termtuple=tuple(sorted((term1,term2)))
			if termtuple in bpsimilaritydict:
				termmatchic.append(bpsimilaritydict[termtuple])
			
			else:
				micaic,mica=getmicaic(term1,term2,ancestors,icdict)
				termmatchic.append((term2,micaic,mica))
				bpsimilaritydict[termtuple]=(term2,micaic,mica)
		bestmatch,bestic,bestmica=getmax(termmatchic)
		bestmatchiclist.append(bestic)
		matchdata.append((term1,bestmatch,bestmica,bestic))	
	mediansim=np.median(bestmatchiclist)
	return mediansim,matchdata,bpsimilaritydict,bestmatchiclist

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


def compute_without_substitution(queryprofiles,dbprofiles,ancestors,icdict,bpsimilaritydict):
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
			bpmediansim,matchdata,bpsimilaritydict,bestmatchiclist=calculate_bestpairs(query,dbprofile,icdict,ancestors,bpsimilaritydict)
			resulttuple=(bpmediansim,matchdata,bestmatchiclist,databaseprofileid)
			results.append(resulttuple)
			allsimilarityscores.append(bpmediansim)
		bestmatch=max(results,key=itemgetter(0))[3]
		similarity= max(results,key=itemgetter(0))[0]
		iclist=",".join(str(x) for x in max(results,key=itemgetter(0))[2])
		bestmatchdata=max(results,key=itemgetter(0))[1]
		similarityscores.append(similarity)
	return similarityscores,allsimilarityscores,querysizes,dbprofilesizes,bpsimilaritydict


def compute_with_substitution(queryprofileids,randomprofiles,annotationpool,ancestors,icdict,queryprofilesize,numqueryprofiles,bpsimilaritydict):
	selfmatch=dict()
	bestsimilaritylist=[]
	quantiles=dict()
	bestmatchsimilarity=dict()
	
	
	for profileid in queryprofileids:
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
				bpmediansim,matchdata,bpsimilaritydict,bestmatchiclist=calculate_bestpairs(queryprofile,dbprofile,icdict,ancestors,bpsimilaritydict)
				resulttuple=(bpmediansim,matchdata,bestmatchiclist,databaseprofileid)
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

	
	return bpsimilaritydict,bestsimilaritylist,quantiles

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


def check_ic(annotationpool,ancestors,icdict):
	print "checking"
	for annotation in set(annotationpool):
		annic=icdict[annotation]
		for anc in ancestors[annotation]:
			ancic=icdict[anc]
			if ancic > annic:
				print annotation,annic,anc,ancic

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


def compute_ic_profile_frequency(infile):
	icdict=dict()
	frequency=dict()
	corpussize=0
	for line in infile:
		profileid,profileancestors=line.strip().split("\t")
		profileancestors=profileancestors.split(",")
		for term in profileancestors:
			if term not in frequency:
				frequency[term]=0
			frequency[term]+=1
		corpussize+=1
	infile.close()
	maxic=round(-math.log(1/corpussize),2)
	for term in frequency:
		ic=round((-math.log(frequency[term]/corpussize))/maxic,2)
		icdict[term]=ic

	return icdict,frequency	
		
def prepare_profile_corpus(dbprofiles,ancestors):
	outfile=open("../data/Corpus_Profile_Frequency.txt",'w')
	for profileid in dbprofiles:
		profileancestors=set()
		for annotation in set(dbprofiles[profileid]):
			ancset=ancestors[annotation]
			profileancestors.add(annotation)
			profileancestors=set.union(profileancestors,ancset)
		outfile.write(profileid+"\t"+",".join(profileancestors)+"\n")
	outfile.close()





	
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
	
	colors=['b','r','g']
	i=0
	for filename in [x for x in os.listdir(figuredir) if ".tsv" in x]:
		fig = plt.figure()
		infile=open(figuredir+filename)
		scoredict=dict()
		expectdict=dict()
		numreplacedset=set()
		meanexpect=[]
		meansim=[]
		error=[]
		iclistdict=dict()
		meanfirstquartile=[]
		meanthirdquartile=[]
		for line in infile:
			if "Number" not in line:
				data=line.strip().split("\t")
				profileid,numreplaced,sim,iclist,expect=data[0],int(data[1]),float(data[3]),data[4],float(data[5])
				numreplacedset.add(numreplaced)
				if numreplaced not in scoredict:
					scoredict[numreplaced]=[]
				scoredict[numreplaced].append(sim)
				if numreplaced not in expectdict:
					expectdict[numreplaced]=[]
				expectdict[numreplaced].append(sim)

				if numreplaced not in iclistdict:
					iclistdict[numreplaced]=[]
				iclistdict[numreplaced].append([float(x) for x in iclist.split(",")])
		

		
		for numreplaced in scoredict:
			meansim.append(np.mean(scoredict[numreplaced]))
			meanexpect.append(np.mean(expectdict[numreplaced]))
			error.append(2*(np.std(scoredict[numreplaced])/math.sqrt(len(scoredict[numreplaced]))))
			firstquartilelist=[]
			thirdquartilelist=[]
			for iclist in iclistdict[numreplaced]:
				firstquartilelist.append(np.percentile(iclist,25))
				thirdquartilelist.append(np.percentile(iclist,75))
			meanfirstquartile.append(np.mean(firstquartilelist))
			meanthirdquartile.append(np.mean(thirdquartilelist))


		legend=[]	
		plt.errorbar(list(numreplacedset),meansim,yerr=error,color=colors[i])
		plt.plot(list(numreplacedset),meanfirstquartile,'--',color=colors[i])
		plt.plot(list(numreplacedset),meanthirdquartile,':',color=colors[i])
		#legend.append("profile size: "+str(numreplaced)+" median")
		#legend.append("profile size: "+str(numreplaced)+" first quartiles")
		#legend.append("profile size: "+str(numreplaced)+" thirdquartiles")

		legend.append( "median")
		legend.append("first quartiles")
		legend.append("thirdquartiles")
		i+=1


		plt.ylim(0,1)
		plt.legend(legend, loc='lower right')

		
		plt.xlabel('Number of annotations replaced')
		plt.ylabel('Median best pair similarity')
		fig.savefig("../results/DecaySimilarity_"+str(numreplaced)+".png", dpi=1200)

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




def plot_ic_comparison():
	fig,ax=plt.subplots()
	colors=['b','g','r','b','g','r']
	#Query profile ID	Number of annotations replaced	Best taxon match	Median similarity	IC list	Expect score
	i=0
	for filename in [x for x in os.listdir("../results/Decay/") if "Based" in x]:
		infile=open("../results/Decay/"+filename)
		size=int(filename.split("ProfileSize")[1].split("_")[0])
		print filename,size
		quant_results(infile,size)
		infile=open("../results/Decay/ "+filename)
		numreplaced=set()
		expect=[]
		expectdict=dict()
		for line in infile:
			if "Query" not in line:
				data=line.strip().split()
				numreplaced.add(int(data[1]))
				if int(data[1]) not in expectdict:
					expectdict[int(data[1])]=[]
				expectdict[int(data[1])].append(float(data[3]))

		error=[]
	
		for num in expectdict:
			expect.append(np.mean(expectdict[num]))

			error.append(2*(np.std(expectdict[num])/math.sqrt(len(expectdict[num]))))
	

		
		if "ProfileBased" in filename:
			#eb1=plt.errorbar(list(numreplaced),expect,yerr=error,color=colors[i],ls='--')
			#eb1[-1][0].set_linestyle('--')
			plt.plot(list(numreplaced),expect,"--",color=colors[i])
		else:
			#plt.errorbar(list(numreplaced),expect,yerr=error,color=colors[i])
			plt.plot(list(numreplaced),expect,color=colors[i])
		i+=1
	plt.ylim(0,1)
	
	#legend=['Annotation frequency IC', 'Profile frequency IC']
	#plt.legend([red_dot, (red_dot, white_cross)], ["Attr A", "Attr A+B"])
	#plt.legend(legend, loc='upper right')

	dotted_line = mlines.Line2D([], [], linewidth=1, linestyle="--", dashes=(3.7, 2), color='black')
	line = mlines.Line2D([], [], linewidth=1, color='black')
	blueline = mlines.Line2D([], [], linewidth=1, color='blue')
	greenline = mlines.Line2D([], [], linewidth=1, color='green')
	redline = mlines.Line2D([], [], linewidth=1, color='red')

	plt.legend([redline,greenline,blueline,line, dotted_line], ["Size 40","Size 20","Size 10","Annotation based frequency", "Profile based frequency"],loc="lower right")


	plt.xlabel('Num of annotations replaced')
	plt.ylabel('Median Best pairs similarity')
	plt.savefig("../results/ICComparison.png",dpi=1200)

def do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,bpsimilaritydict,numofdbprofiles,filename):
		
	out=open(filename,'w')
	print "Profile size",queryprofilesize
	
	
	bpsimilaritydict,bestsimilaritylist,quantiles=compute_with_substitution(queryprofileids,dbprofiles,annotationpool,ancestors,icdict,queryprofilesize,numqueryprofiles,bpsimilaritydict)

	noisequeryprofiles=create_random_profiles(queryprofilesize,annotationpool,5)
	similarityscores,allsimilarityscores,querysizes,dbprofilesizes,bpsimilaritydict=compute_without_substitution(noisequeryprofiles,dbprofiles,ancestors,icdict,bpsimilaritydict)

	print "Mean noise to noise similarity",np.mean(similarityscores)

	#results,box_query,box_db,querylambda_,dblambda_=regression(allsimilarityscores,querysizes,dbprofilesizes)

	results=load_regression()


	resstddev=np.std(results.resid)

	out.write("Query profile ID\tNumber of annotations replaced\tBest taxon match\tMedian similarity\tIC list\tExpect score\n")
	for tup in bestsimilaritylist:
		query=tup[0]
		numreplaced=tup[1]
		bestmatch=tup[2]
		similarity=tup[3]
		iclist=tup[4]
		dbprofilesize=len(dbprofiles[bestmatch])
		#boxquerysize = boxcox(np.array([queryprofilesize]) + 1,lmbda=querylambda_)[0]
		#boxdbsize = boxcox(np.array([dbprofilesize]) + 1,lmbda=dblambda_)[0]
		
		boxquerysize = np.log(queryprofilesize)
		boxdbsize = np.log(dbprofilesize)

		expectscore=compute_expect_scores(results,similarity,boxquerysize,boxdbsize,resstddev,numofdbprofiles)
		out.write(query+"\t"+str(numreplaced)+"\t"+bestmatch+"\t"+str(similarity)+"\t"+ iclist+"\t"+str(expectscore)+"\n")
	out.close()
	return bpsimilaritydict

def load_regression():
	infile=open("../results/test_lm.txt")
	querysizes=[]
	dbsizes=[]
	scores=[]
	for line in infile:
		querysize,dbsize,score=line.strip().split("\t")
		querysizes.append(float(querysize))
		dbsizes.append(float(dbsize))
		scores.append(float(score))

	sizes = [querysizes,dbsizes]

	ones = np.ones(len(sizes[0]))
	X = sm.add_constant(np.column_stack((sizes[0], ones)))
	for ele in sizes[1:]:
		X = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(scores, X).fit()
	


	predicted = results.fittedvalues
	infile.close()
	return results

def compare_ic(dbprofiles,ancestors,numqueryprofiles,numofdbprofiles,annotationpool):
		bpsimilaritydict_profile=dict()
		bpsimilaritydict_annotation=dict()
		profileicdict=dict()
		prepare_profile_corpus(dbprofiles,ancestors)
		infile=open("../data/Corpus_Profile_Frequency.txt")
		profileicdict,profilefrequency=compute_ic_profile_frequency(infile)
		infile.close()
		print "Computed profile based IC"
		annotationicdict=dict()
		annotationicdict,annotationfrequency=compute_ic(dbprofiles,ancestors)
		print "Computed annotation based IC"

		icfile=open("IC.tsv",'w')
		for term in annotationicdict:
			icfile.write(str(annotationicdict[term])+"\t"+str(profileicdict[term])+"\n")
		icfile.close()
		

		
		

		queryprofilesize=10
		profileids=set()
		for profileid in dbprofiles:
			if len(dbprofiles[profileid])==queryprofilesize:
				profileids.add(profileid)
		queryprofileids=random.sample(profileids, numqueryprofiles)
		filename="../results/Decay/ProfileBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
		bpsimilaritydict_profile=do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,profileicdict,numqueryprofiles,bpsimilaritydict_profile,numofdbprofiles,filename)


		filename="../results/Decay/AnnotationBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
		bpsimilaritydict_annotation=do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,annotationicdict,numqueryprofiles,bpsimilaritydict_annotation,numofdbprofiles,filename)

		print "Done", queryprofilesize

		queryprofilesize=20
		profileids=set()
		for profileid in dbprofiles:
			if len(dbprofiles[profileid])==queryprofilesize:
				profileids.add(profileid)
		queryprofileids=random.sample(profileids, numqueryprofiles)
		filename="../results/Decay/ProfileBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
		bpsimilaritydict_profile=do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,profileicdict,numqueryprofiles,bpsimilaritydict_profile,numofdbprofiles,filename)

		filename="../results/Decay/AnnotationBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
		bpsimilaritydict_annotation=do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,annotationicdict,numqueryprofiles,bpsimilaritydict_annotation,numofdbprofiles,filename)

		print "Done", queryprofilesize

		queryprofilesize=40
		profileids=set()
		for profileid in dbprofiles:
			if len(dbprofiles[profileid])==queryprofilesize:
				profileids.add(profileid)
		queryprofileids=random.sample(profileids, numqueryprofiles)
		filename="../results/Decay/ProfileBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
		bpsimilaritydict_profile=do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,profileicdict,numqueryprofiles,bpsimilaritydict_profile,numofdbprofiles,filename)

		filename="../results/Decay/AnnotationBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
		bpsimilaritydict_annotation=do(queryprofileids,queryprofilesize,dbprofiles,annotationpool,ancestors,annotationicdict,numqueryprofiles,bpsimilaritydict_annotation,numofdbprofiles,filename)
		print "Done", queryprofilesize


def main():
	
	
	queryprofilesize=sys.argv[1]
	numqueryprofiles=int(sys.argv[2])
	icflag=sys.argv[3]
	numofdbprofiles=659
	plot_ic_comparison()
	sys.exit()


	bpsimilaritydict=dict()
	dbprofiles,annotationpool=load_randomprofiles()
	print "Loaded randomprofiles"
	ancestors=load_ancestors() 
	print "Loaded ancestors"
	compare_ic(dbprofiles,ancestors,numqueryprofiles,numofdbprofiles,annotationpool)
	sys.exit()

	




	if icflag =="profilebased":
		prepare_profile_corpus(dbprofiles,ancestors)
		infile=open("../data/Corpus_Profile_Frequency.txt")
		icdict=compute_ic_profile_frequency(infile)
		infile.close()
		print "Computed profile based IC"
	else:
		icdict,frequency=compute_ic(dbprofiles,ancestors)
		print "Computed annotation based IC"
	
	bpsimilaritydict=do(10,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,bpsimilaritydict,numofdbprofiles)
	sys.exit()
	print "Done with 10"
	bpsimilaritydict=do(20,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,bpsimilaritydict,numofdbprofiles)
	print "Done with 20"
	bpsimilaritydict=do(40,dbprofiles,annotationpool,ancestors,icdict,numqueryprofiles,bpsimilaritydict,numofdbprofiles)
	print "Done with 40"

	plot_decay()
	sys.exit()


	



if __name__ == "__main__":
	import random
	import matplotlib.lines as mlines
	import os
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
