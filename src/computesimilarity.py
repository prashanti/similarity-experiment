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
				apsimilaritydict[termtuple]=(bpsimilaritydict[termtuple][1],bpsimilaritydict[termtuple][2])
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
		bestmatchiclist.append(bestic)
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
	return np.mean([mediansim1,mediansim2]),bpsimilaritydict


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
	
	return icdict

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


def load_subclasses(profileterms):
	subclasses=dict()
	infile=open("../data/Subclasses_EAttr(Q)_OldRealProfiles.txt")
	for line in infile:
		term,subclass=line.replace(">","").replace("<","").strip().split("\t")
		if term in profileterms:
			if term not in subclasses:
				subclasses[term]=set()
				subclasses[term].add(term)
			if subclass !="owl:Nothing":		
				subclasses[term].add(subclass)
	infile.close()
	return subclasses

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





	
def compute_expect_scores(results,actualsimilarity,queryprofilesize,dbprofilesize, querylambda_,dblambda_,numofdbprofiles):
	# now compute expect scores for the individual comparison

	dbcoeff=results.params[0]
	querycoeff=results.params[1]
	constant=results.params[2]

	resstddev=np.std(results.resid)

	boxquerysize = boxcox(np.array([queryprofilesize]) + 1,lmbda=querylambda_)[0]
	boxdbsize = boxcox(np.array([dbprofilesize]) + 1,lmbda=dblambda_)[0]
	
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


def error(scorelist):
	return 2*(np.std(scorelist)/math.sqrt(len(scorelist)))

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
	
	return results,querylambda_,dblambda_



def regression_allmetrics(dbprofiles):
	reg_results=dict()
	df=pd.read_csv('../results/RegressionInputs.tsv', sep='\t')
	bp_sym_mediansim_pic,bp_asym_mediansim_pic,bp_sym_mediansim_aic,bp_asym_mediansim_aic,ap_mediansim_pic,ap_mediansim_aic,bp_asym_simj,bp_sym_simj,ap_simj=df['bp_sym_mediansim_pic'],df['bp_asym_mediansim_pic'],df['bp_sym_mediansim_aic'],df['bp_asym_mediansim_aic'],df['ap_mediansim_pic'],df['ap_mediansim_aic'],df['bp_asym_simj'],df['bp_sym_simj'],df['ap_simj']
	dbprofilesizes=[len(dbprofiles[x]) for x in df['dbprofileid']]
	queryprofilesizes=[len(dbprofiles[x]) for x in df['queryprofileid']]

	
	reg_results['bp_sym_pic']=(regression(bp_sym_mediansim_pic.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['bp_asym_pic']=(regression(bp_asym_mediansim_pic.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['bp_sym_aic']=(regression(bp_sym_mediansim_aic.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['bp_asym_aic']=(regression(bp_asym_mediansim_aic.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['ap_pic']=(regression(ap_mediansim_pic.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['ap_aic']=(regression(ap_mediansim_aic.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['bp_asym_simj']=(regression(bp_asym_simj.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['bp_sym_simj']=(regression(bp_sym_simj.tolist(),queryprofilesizes,dbprofilesizes))
	reg_results['ap_simj']=(regression(ap_simj.tolist(),queryprofilesizes,dbprofilesizes))

	return reg_results

def compare_probability(similarityscore,noisescores):
	greaterscores=[x for x in noisescores if x>=similarityscore]
	empiricalpvalue=len(greaterscores)/len(noisescores)
	return empiricalpvalue




def compute_expect_allmetrics(resfile,noisefile,reg_results,dbprofiles,numofdbprofiles,profilesize):

	infile=open(resfile)
	out=open(resfile.replace("/results/","/results/Expect_"),'w')
	out.write("queryid\tnumber of annotations replaced\tmatch_bp_sym_pic\tscore_bp_sym_pic\texpect_bp_sym_pic\tmatch_bp_sym_aic\tscore_bp_sym_aic\texpect_bp_sym_aic\tmatch_bp_asym_pic\tscore_bp_asym_pic\texpect_bp_asym_pic\tmatch_bp_asym_aic\tscore_bp_asym_aic\texpect_bp_asym_aic\tmatch_ap_pic\tscore_ap_pic\texpect_ap_pic\tmatch_ap_aic\tscore_ap_aic\texpect_ap_pic\tmatch_bp_asym_simj\tscore_bp_asym_simj\texpect_bp_asym_simj\tmatch_bp_sym_simj\tscore_bp_sym_simj\texpect_bp_sym_simj\tmatch_ap_simj\tscore_ap_simj\texpect_ap_simj\n")
	
	noisedf=pd.read_csv(noisefile, sep='\t')
	metrics=['bp_sym_pic','bp_sym_aic','bp_asym_pic','bp_asym_aic','ap_pic','ap_aic','bp_asym_simj','bp_sym_simj','ap_simj']
	for metric in metrics:
		comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'w')
		comp.write("Empirical p-value\tRegression p-value\n")
		comp.close()

	for line in infile:
		if "QueryID" not in line:
			queryid,numreplaced,match_bp_sym_pic,score_bp_sym_pic,match_bp_sym_aic,score_bp_sym_aic,match_bp_asym_pic,score_bp_asym_pic,match_bp_asym_aic,score_bp_asym_aic,match_ap_pic,score_ap_pic,match_ap_aic,score_ap_aic,match_bp_asym_simj,score_bp_asym_simj,match_bp_sym_simj,score_bp_sym_simj,match_ap_simj,score_ap_simj=line.strip().split("\t")

			queryprofilesize=dbprofilesize=profilesize




			metric='bp_sym_pic'

			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')

			#dbprofilesize=len(dbprofiles[match_bp_sym_pic])
			expect_bp_sym_pic=compute_expect_scores(reg_results[metric][0],float(score_bp_sym_pic),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_bp_sym_pic),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_bp_sym_pic/numofdbprofiles)+"\n")
			comp.close()


			metric='bp_sym_aic'

			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')

			#dbprofilesize=len(dbprofiles[match_bp_sym_aic])
			expect_bp_sym_aic=compute_expect_scores(reg_results[metric][0],float(score_bp_sym_aic),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_bp_sym_aic),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_bp_sym_aic/numofdbprofiles)+"\n")
			comp.close()


			metric='bp_asym_pic'
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_bp_asym_pic])
			expect_bp_asym_pic=compute_expect_scores(reg_results[metric][0],float(score_bp_asym_pic),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_bp_asym_pic),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_bp_asym_pic/numofdbprofiles)+"\n")
			comp.close()


			metric='bp_asym_aic'
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_bp_asym_aic])
			expect_bp_asym_aic=compute_expect_scores(reg_results[metric][0],float(score_bp_asym_aic),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_bp_asym_aic),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_bp_asym_aic/numofdbprofiles)+"\n")
			comp.close()


			metric='ap_pic'
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_ap_pic])
			expect_ap_pic=compute_expect_scores(reg_results[metric][0],float(score_ap_pic),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_ap_pic),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_ap_pic/numofdbprofiles)+"\n")
			comp.close()

			metric='ap_aic'
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_ap_aic])
			expect_ap_aic=compute_expect_scores(reg_results[metric][0],float(score_ap_aic),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_ap_aic),noisedf[metric].tolist())
			
			comp.write(str(empiricalpvalue)+"\t"+str(expect_ap_aic/numofdbprofiles)+"\n")
			comp.close()


			metric='bp_asym_simj'
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_bp_asym_simj])
			expect_bp_asym_simj=compute_expect_scores(reg_results[metric][0],float(score_bp_asym_simj),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_bp_asym_simj),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_bp_asym_simj/numofdbprofiles)+"\n")
			comp.close()

			metric='bp_sym_simj'
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_bp_sym_simj])
			expect_bp_sym_simj=compute_expect_scores(reg_results[metric][0],float(score_bp_sym_simj),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_bp_sym_simj),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_bp_sym_simj/numofdbprofiles)+"\n")
			comp.close()


			metric='ap_simj'		
			comp=open("../results/EmpiricalvsRegression_"+metric+".tsv",'a')
			#dbprofilesize=len(dbprofiles[match_ap_simj])
			expect_ap_simj=compute_expect_scores(reg_results[metric][0],float(score_ap_simj),queryprofilesize,dbprofilesize, reg_results[metric][1],reg_results[metric][2],numofdbprofiles)
			empiricalpvalue=compare_probability(float(score_ap_simj),noisedf[metric].tolist())
			comp.write(str(empiricalpvalue)+"\t"+str(expect_ap_simj/numofdbprofiles)+"\n")
			comp.close()

			out.write(queryid+"\t"+numreplaced+"\t"+match_bp_sym_pic+"\t"+score_bp_sym_pic+"\t"+ str(expect_bp_sym_pic)+"\t"+  match_bp_sym_aic+"\t"+score_bp_sym_aic+"\t" +  str(expect_bp_sym_aic)+"\t"    +match_bp_asym_pic+"\t"+score_bp_asym_pic+"\t"+   str(expect_bp_asym_pic)+"\t"+    match_bp_asym_aic+"\t"+score_bp_asym_aic+"\t"+     str(expect_bp_asym_aic)+"\t"+    match_ap_pic+"\t"+score_ap_pic+"\t"+   str(expect_ap_pic)+"\t"+   match_ap_aic+"\t"+score_ap_aic+"\t"+    str(expect_ap_pic)+"\t"+    match_bp_asym_simj+"\t"+score_bp_asym_simj+"\t"+     str(expect_bp_asym_simj)+"\t"+ match_bp_sym_simj+"\t"+score_bp_sym_simj+"\t"+   str(expect_bp_sym_simj)+"\t"+  match_ap_simj+"\t"+score_ap_simj+"\t"+str(expect_ap_simj)+"\n")


	infile.close()
	out.close()
	comp.close()


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
def plot_signal_noise(comparisonfile,profilesize):
	numreplacedlist=[]
	fig = plt.figure()
	noise=dict()
	signal=dict()
	legend=set()
	results=dict()
	colors=['b','b','g','g','r','r']
	
	#VTO_0037490	4	AP_AIC x y difference
	numreplacedlist=list(range(0,profilesize+1))
	
	for line in comparisonfile:
		data=line.strip().split("\t")
		queryid,numreplaced,metric,difference=data[0],data[1],data[2],data[5]
		if queryid not in results:
			results[queryid]=dict()
		if metric not in results[queryid]:
			results[queryid][metric]=dict()
		results[queryid][metric][int(numreplaced)]=float(difference)
	
	i=0
	metriclist=['BP_Sym_AIC','BP_Sym_PIC','BP_Asym_AIC','BP_Asym_PIC','AP_AIC','AP_PIC']
	#for queryid in results:
	queryid='VTO_0052581'
	for metric in metriclist:
		differencelist=[]
		for numreplaced in numreplacedlist:
			differencelist.append(results[queryid][metric][numreplaced])
		if "PIC" in metric:
			plt.plot(numreplacedlist,differencelist,"--",color=colors[i])
		else:
			plt.plot(numreplacedlist,differencelist,color=colors[i])
		i+=1

			

	plt.legend(metriclist, loc='lower right')	
	plt.xlabel('Number of annotations replaced')
	plt.ylabel('Best Match Rank Difference')
	fig.savefig("../results/Decay/SignalNoiseDifference.png", dpi=1200)
def compareranks(signalscores,noisescores):
	allscores=sorted(signalscores,reverse=True)+sorted(noisescores,reverse=True)
	allranks=rankdata(allscores)
	signalranks=allranks[:len(signalscores)]
	noiseranks=allranks[-len(noisescores):]
	signalrank=allranks[0]
	noiserank=allranks[len(signalscores)]
	print "Scores,",allscores[0],allscores[len(signalscores)]
	print "Ranks,",signalrank,noiserank
	print allscores
	print "Ranks"
	print allranks
	return signalranks,noiseranks,signalrank-noiserank



def compare_signal_noise(noisefile,signalfile,profilesize):
	fig=plt.figure()
	colors=['b','g','r','m','k','c']
	legend=[]
	out=open("../results/Decay/SignalvsNoise_"+str(profilesize)+".tsv",'w')
	noise=dict()
	signal=dict()
	for line in noisefile:
		data=line.strip().split("\t")
		noise[data[2]]=[float(x) for x in data[3].split(",")]
	for line in signalfile:
		data=line.strip().split("\t")
		if data[0] not in signal:
			signal[data[0]]=dict()
		if int(data[1]) not in signal[data[0]]:
			signal[data[0]][int(data[1])]=dict()
		signal[data[0]][int(data[1])][data[2]] =[float(x) for x in data[3].split(",")]
	i=0
	signalbestscores=dict()
	noisebestscores=dict()
	for queryid in signal:
		for numreplaced in sorted(signal[queryid]):
			for metric in signal[queryid][numreplaced]:
				signalscores=signal[queryid][numreplaced][metric]
				noisescores=noise[metric]
				if metric not in signalbestscores:
					signalbestscores[metric]=[]
				signalbestscores[metric].append(np.max(signalscores))
				if metric not in noisebestscores:
					noisebestscores[metric]=[]
				noisebestscores[metric].append(np.max(noisescores))
				print metric,numreplaced
				signalranks,noiseranks,difference=compareranks(signalscores,noisescores) 
				out.write(queryid+"\t"+str(numreplaced)+"\t"+metric+"\t"+",".join(str(x) for x in signalranks )+"\t"+",".join(str(x) for x in noiseranks)+"\t"+str(difference)+"\n")


		
		# for metric in signalbestscores:
		# 	plt.plot(list(range(0,21)),signalbestscores[metric],color=colors[i])
		# 	i+=1
		# 	legend.append(metric)
		# i=0
		# for metric in signalbestscores:
		# 	plt.plot(list(range(0,21)),noisebestscores[metric],"--",color=colors[i])
		# 	i+=1
		# plt.legend(legend, loc='upper right',ncol=2)
		# plt.xlabel('Number of annotations replaced')
		# plt.ylabel('Similarity of best match')
		# plt.show()



		# break
	out.close()

def signal_vs_noise(signalfile,noisefile):

	noisedf=pd.read_csv(noisefile, sep='\t')



	signalbestscores=dict()
	noisebestscores=dict()
	numreplacedset=set()
	metricset=set()
	signalmeanbestscores=dict()
	signalfile.next()
	for line in signalfile:
		queryid,numreplaced,match_bp_sym_pic,bp_sym_pic,match_bp_sym_aic,bp_sym_aic,match_bp_asym_pic,bp_asym_pic,match_bp_asym_aic,bp_asym_aic,match_ap_pic,ap_pic,match_ap_aic,ap_aic,match_bp_asym_simj,bp_asym_simj,match_bp_sym_simj,bp_sym_simj,match_ap_simj,ap_simj=line.strip().split("\t")
		numreplaced=int(numreplaced)
	

		if numreplaced not in signalbestscores:
			signalbestscores[numreplaced]=dict()
			numreplacedset.add(numreplaced)
		
		metric="bp_sym_pic"

		metricset.add(metric)
		noisebestscores[metric]=noisedf[metric]
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_sym_pic))

		metric="bp_sym_aic"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_sym_aic))


		metric="bp_asym_pic"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_asym_pic))

		metric="bp_asym_aic"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_asym_aic))


		metric="ap_pic"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(ap_pic))

		metric="ap_aic"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(ap_aic))

		metric="bp_sym_simj"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_sym_simj))

		metric="bp_asym_simj"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_asym_simj))

		metric="ap_simj"
		noisebestscores[metric]=noisedf[metric]
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(ap_simj))

	signalfile.close()
	out=open("../results/MetricComparison_SignalvsNoise.tsv",'w')
	out.write("Number of annotations replaced\tMetric\tMean signal best score\t2*std error\tPercentile rank in noise\n")
	fig=plt.figure()
	for metric in sorted(metricset):
		signallist=[]
		errorlist=[]
		percentilelist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
			percentilelist.append(stats.percentileofscore(noisebestscores[metric],np.mean(signalbestscores[numreplaced][metric])))
			out.write(str(numreplaced)+"\t"+metric+"\t"+str(np.mean(signalbestscores[numreplaced][metric]))+"\t"+str(error(signalbestscores[numreplaced][metric]))+"\t"+str(stats.percentileofscore(noisebestscores[metric],np.mean(signalbestscores[numreplaced][metric])))+"\n")

		
	title=dict()
	title['ap_aic']="Annotation IC"
	title['ap_pic']="All Pairs \n\n Profile IC"
	title['ap_simj']="Jaccard"
	title['bp_asym_aic']="Annotation IC"
	title['bp_asym_pic']="Best Pairs Asymmetric\n\nProfile IC"
	title['bp_asym_simj']="Jaccard"
	title['bp_sym_aic']="Annotation IC"
	title['bp_sym_pic']="Best Pairs Symmetric\n\nProfile IC"
	title['bp_sym_simj']="Jaccard"
	f, axarr = plt.subplots(3, 3)
	i=j=0
	for metric in sorted(metricset):
		signallist=[]
		errorlist=[]
		percentilelist=[np.percentile(noisebestscores[metric],99),np.percentile(noisebestscores[metric],99.9)]
		alphalist=['#B6B6B4','#736F6E','#0C090A']
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color='black')
		axarr[i,j].axhline(y=np.percentile(noisebestscores[metric],99),linestyle='--',color='orange')
		axarr[i,j].axhline(y=np.percentile(noisebestscores[metric],99.9),linestyle='--',color='blue',label='x')
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1.1)
		axarr[i, j].set_xlim(0,11)
		
		j+=1
		if j==3:
			i+=1
			j=0
		if i==3:
			i=0	
	plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
	plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	#plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)	
	#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)	
	axarr[2,1].set_xlabel('Number of annotations replaced')
	axarr[1,0].set_ylabel('Similarity score')
	plt.tight_layout()	
	
	plt.savefig("../results/MetricComparison_DecayvsNoise_Quantiles.png")
	
	out.close()



	

def run_noise_to_noise(dbprofiles,ancestors,queryprofilesize,annotationpool,profileicdict,annotationicdict):
	bpsimilaritydict_profile=dict()
	bpsimilaritydict_annotation=dict()
	apsimilaritydict_profile=dict()
	apsimilaritydict_annotation=dict()
	noisequeryprofiles=create_random_profiles(queryprofilesize,annotationpool,10)
	out=open("../results/Decay/NoisetoNoise_ProfileSize"+str(queryprofilesize)+"_Results.tsv",'w')
	#outfile=open("../results/Decay/AllScores_NoisetoNoise_"+str(queryprofilesize)+"_Results.tsv",'w')
	out.write("QueryID\tNumber of annotations replaced\tBest Match Similarity--BP Symmetric Profile IC\tBest Match Similarity--BP Symmetric Annotation IC\tBest Match Similarity--BP Asymmetric Profile IC\tBest Match Similarity--BP Asymmetric Annotation IC\tBest Match Similarity--AP Profile IC\tBest Match Similarity--AP Annotation IC\n")
	i=1
	for queryprofile in noisequeryprofiles:
		print i
		queryid="Query_"+str(i)
		i+=1
		bp_asym_pic_results=[]
		bp_asym_aic_results=[]
		bp_sym_pic_results=[]
		bp_sym_aic_results=[]
		ap_pic_results=[]
		ap_aic_results=[]
		for dbprofileid in dbprofiles:
			profile2=dbprofiles[dbprofileid]

			bp_sym_mediansim_pic,bpsimilaritydict_profile=calculate_bestpairs_symmetric(queryprofile,profile2,profileicdict,ancestors,bpsimilaritydict_profile)
			bp_sym_pic_results.append((bp_sym_mediansim_pic,dbprofileid))
			
			bp_asym_mediansim_pic,matchdata,bpsimilaritydict_profile,bestmatchiclist=calculate_bestpairs_asymmetric(queryprofile,profile2,profileicdict,ancestors,bpsimilaritydict_profile)
			bp_asym_pic_results.append((bp_asym_mediansim_pic,dbprofileid))

			bp_sym_mediansim_aic,bpsimilaritydict_annotation=calculate_bestpairs_symmetric(queryprofile,profile2,annotationicdict,ancestors,bpsimilaritydict_annotation)
			bp_sym_aic_results.append((bp_sym_mediansim_aic,dbprofileid))

			bp_asym_mediansim_aic,matchdata,bpsimilaritydict_annotation,bestmatchiclist=calculate_bestpairs_asymmetric(queryprofile,profile2,annotationicdict,ancestors,bpsimilaritydict_annotation)
			bp_asym_aic_results.append((bp_asym_mediansim_aic,dbprofileid))
			

			ap_mediansim_pic,apsimilaritydict_profile=get_allpairs_medianic(queryprofile,profile2,profileicdict,ancestors,bpsimilaritydict_profile,apsimilaritydict_profile)
			ap_pic_results.append((ap_mediansim_pic,dbprofileid))


			ap_mediansim_aic,apsimilaritydict_annotation=get_allpairs_medianic(queryprofile,profile2,annotationicdict,ancestors,bpsimilaritydict_annotation,apsimilaritydict_annotation)
			ap_aic_results.append((ap_mediansim_aic,dbprofileid))


		
		bestmatch_bp_sym_pic=max(bp_sym_pic_results,key=itemgetter(0))[1]	
		bestmatchsimilarity_bp_sym_pic=max(bp_sym_pic_results,key=itemgetter(0))[0]
		#allresults=[y[0] for y in bp_sym_pic_results]
		# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Sym_PIC\t"+",".join(str(x)for x in allresults)+"\n")
		
		bestmatch_bp_sym_aic=max(bp_sym_aic_results,key=itemgetter(0))[1]	
		bestmatchsimilarity_bp_sym_aic=max(bp_sym_aic_results,key=itemgetter(0))[0]
		# allresults=[y[0] for y in bp_sym_aic_results]
		# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Sym_AIC\t"+",".join(str(x)for x in allresults)+"\n")


		bestmatch_bp_asym_pic=max(bp_asym_pic_results,key=itemgetter(0))[1]	
		bestmatchsimilarity_bp_asym_pic=max(bp_asym_pic_results,key=itemgetter(0))[0]
		# allresults=[y[0] for y in bp_asym_pic_results]
		# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Asym_PIC\t"+",".join(str(x)for x in allresults)+"\n")
		
		bestmatch_bp_asym_aic=max(bp_asym_aic_results,key=itemgetter(0))[1]	
		bestmatchsimilarity_bp_asym_aic=max(bp_asym_aic_results,key=itemgetter(0))[0]
		# allresults=[y[0] for y in bp_asym_aic_results]
		# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Asym_AIC\t"+",".join(str(x)for x in allresults)+"\n")
			
		bestmatch_ap_pic=max(ap_pic_results,key=itemgetter(0))[1]	
		bestmatchsimilarity_ap_pic=max(ap_pic_results,key=itemgetter(0))[0]
		# allresults=[y[0] for y in ap_pic_results]
		# outfile.write(queryid+"\t"+str(i)+"\t"+"AP_PIC\t"+",".join(str(x)for x in allresults)+"\n")

		bestmatch_ap_aic=max(ap_aic_results,key=itemgetter(0))[1]	
		bestmatchsimilarity_ap_aic=max(ap_aic_results,key=itemgetter(0))[0]
		# allresults=[y[0] for y in ap_aic_results]
		# outfile.write(queryid+"\t"+str(i)+"\t"+"AP_AIC\t"+",".join(str(x)for x in allresults)+"\n")
		
		
		
		out.write(queryid+"\t"+str(i)+"\t"+str(bestmatchsimilarity_bp_sym_pic)+"\t"+str(bestmatchsimilarity_bp_sym_aic)+"\t"+str(bestmatchsimilarity_bp_asym_pic)+"\t"+str(bestmatchsimilarity_bp_asym_aic)+"\t"+str(bestmatchsimilarity_ap_pic)+"\t"+str(bestmatchsimilarity_ap_aic)+"\n")


		
	out.close()
	#outfile.close()



def load_hrss_precomputation():
	distancedict=dict()
	micadict=dict()
	infile=open("../data/Distance_HRSS1.tsv")
	for line in infile:
		data=line.split("\t")
		term1,mica,distance=data[0],data[1],data[2]
		if int(distance) !=0:
			if ((term1,mica)) not in distancedict:
				distancedict[(term1,mica)] =distance
	infile.close()
	return distancedict
def experiment_hrss(decayedprofiles,dbprofiles,queryprofilesize,ancestors,icflag,icdict,subclasses):
	distancedict=load_hrss_precomputation()
	print "loaded distance"
	out=open("../results/Decay/HRSS_ProfileSize"+str(queryprofilesize)+"_"+icflag+"_Results.tsv",'w')

	out.write("QueryID\tNumber of annotations replaced\tAP Best Match\tAP Best Match Similarity\tBP Best Match\tBP Best Match Similarity\n")

	
	for queryid in decayedprofiles:
		for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
			ap_hrss_results=[]
			bp_hrss_results=[]
			queryprofile=decayedprofiles[queryid][str(numreplaced)]
			for dbprofileid in dbprofiles:
				dbprofile=dbprofiles[dbprofileid]
				ap_hrss,bp_hrss=hrss(queryprofile,dbprofile,distancedict,icdict,ancestors,subclasses)
				ap_hrss_results.append((ap_hrss,dbprofileid))
				bp_hrss_results.append((bp_hrss,dbprofileid))
			bestmatch_ap_hrss=max(ap_hrss_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_ap_hrss=max(ap_hrss_results,key=itemgetter(0))[0]
			bestmatch_bp_hrss=max(bp_hrss_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_bp_hrss=max(bp_hrss_results,key=itemgetter(0))[0]
			out.write(queryid+"\t"+str(numreplaced)+"\t"+   bestmatch_ap_hrss+"\t"+str(bestmatchsimilarity_ap_hrss)+"\t"+bestmatch_bp_hrss+"\t"+str(bestmatchsimilarity_bp_hrss)+"\n")
			print queryid,numreplaced
	out.close()	
	

def get_profile_terms(decayedprofiles,dbprofiles):
	termset=set()
	for profileid in decayedprofiles:
		profile=decayedprofiles[profileid]
		for term in profile:
			termset.add(term)

	for profileid in dbprofiles:
		profile=dbprofiles[profileid]
		for term in profile:
			termset.add(term)

	return termset

def getmicaic(term1,term2,ancestors,icdict):
	micaic=0
	mica=""
	commonancestors=set.intersection(ancestors[term1],ancestors[term2])
	#for anc in commonancestors:
	#	print term2,icdict[anc],anc
	lcslist=[(term2,icdict[anc],anc) for anc in commonancestors]
	match,micaic,mica=getmax(lcslist)
	
	if len(lcslist)>0:
		 
		return micaic,mica
	else:
		return 0,"None"


def hrss(queryprofile,dbprofile,distancedict,icdict,ancestors,subclasses):
	bp_hrssscores=[]
	ap_hrssscores=[]
	scores=dict()
	for term1 in queryprofile:
		for term2 in dbprofile:	
			micaic,mica=getmicaic(term1,term2,ancestors,icdict)
			
			term1dist=0
			term2dist=0
			if (term1,mica) in distancedict:
				term1dist=distancedict[(term1,mica)]
			if (term2,mica) in distancedict:
				term2dist=distancedict[(term2,mica)]

			gamma=int(term1dist)+int(term2dist)
			mil1ic=mil2ic=0
			leaficlist=[]
			if term1 in subclasses:
				for term in subclasses[term1]:
					if term in icdict:
						leaficlist.append(icdict[term])
			
			
			
			if len(leaficlist)>0:
				mil1ic=max(leaficlist)
			else:
				mil1ic=icdict[term1]
			leaficlist=[]
			
			if term2 in subclasses:
				for term in subclasses[term2]:
					if term in icdict:
						leaficlist.append(icdict[term])

			if len(leaficlist)>0:
				mil2ic=max(leaficlist)
			else:
				mil2ic=icdict[term2]

	
			alpha=micaic
			
			beta=(mil1ic-icdict[term1]+mil2ic-icdict[term2])/2


			if alpha==0:
				hrss=0
			else:
				hrss=(1/(1+gamma))*(alpha/(alpha+beta))
			
			if term1 not in scores:
				scores[term1]=dict()
			if term2 not in scores[term1]:
				scores[term1][term2]=hrss
			ap_hrssscores.append(hrss)
	for term1 in scores:
		bp_hrssscores.append(max(scores[term1].iteritems(), key=itemgetter(1))[1])

			
	return np.median(ap_hrssscores),np.median(bp_hrssscores)
	
def groupwise_jaccard(profile1,profile2,ancestors):
	
	ancestors1=set()
	ancestors2=set()
	for term in profile1:
		ancestors1=set.union(ancestors1,ancestors[term])
	for term in profile2:
		ancestors2=set.union(ancestors2,ancestors[term])
	common=set.intersection(ancestors1,ancestors2)
	
	if len(common) > 0:
		union=set.union(ancestors1,ancestors2)	
		simj=len(common)/len(union)

	else:
		simj=0	
	return simj,union

def calculate_bestpairs_jaccard_symmetric(profile1,profile2,ancestors):
	finalsim=0
	bestmatchsimj=[]
	termmatchsimj=[]
	for term1 in profile1:
		termmatchsimj=[]
		for term2 in profile2:
			simj=getsimj(term1,term2,ancestors)
			termmatchsimj.append(simj)
		bestmatchsimj.append(max(termmatchsimj))
	sim=np.median(bestmatchsimj)

	bestmatchsimj=[]
	termmatchsimj=[]
	for term1 in profile2:
		termmatchsimj=[]
		for term2 in profile1:
			simj=getsimj(term1,term2,ancestors)
			termmatchsimj.append(simj)
		bestmatchsimj.append(max(termmatchsimj))
	revsim=np.median(bestmatchsimj)

	finalsim=np.mean([sim,revsim])
	return finalsim



def getsimj(term1,term2,ancestors):
	if len(set.union(ancestors[term1],ancestors[term2])) >0:
		simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
	else:
		simj=0
	return simj


def calculate_bestpairs_jaccard_asymmetric(profile1,profile2,ancestors):
	finalsim=0
	bestmatchsimj=[]
	termmatchsimj=[]
	for term1 in profile1:
		termmatchsimj=[]
		for term2 in profile2:

			simj=getsimj(term1,term2,ancestors)
			termmatchsimj.append(simj)
		bestmatchsimj.append(max(termmatchsimj))
	sim=np.median(bestmatchsimj)
	return sim



	

def calculate_allpairs_jaccard(profile1,profile2,ancestors):
	
	mediansimj=0
	simj=[]
	for term1 in profile1:
		for term2 in profile2:
			simj.append(getsimj(term1,term2,ancestors))		
	return np.median(simj)


def create_decayed_profiles(dbprofiles,numqueryprofiles,queryprofilesize,annotationpool):
	decayedprofiles=dict()
	
	dictout=open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt",'w')
	profileids=set()
	for profileid in dbprofiles:
		if len(dbprofiles[profileid])==queryprofilesize:
			profileids.add(profileid)
	queryprofileids=random.sample(profileids, numqueryprofiles)
	for queryid in queryprofileids:
		queryprofile=deepcopy(dbprofiles[queryid])
		replaced=set()
		for numreplaced in range(0,queryprofilesize+1):
			if queryid not in decayedprofiles:
				decayedprofiles[queryid]=dict()
			decayedprofiles[queryid][int(numreplaced)]=deepcopy(queryprofile)			
			if len(replaced)<len(queryprofile):	
				queryprofile,replaced=substitute_annotation(queryprofile,replaced,annotationpool)
				
	json.dump(decayedprofiles, dictout)	
	dictout.close()		
	
def experiment_groupwise_jaccard(dbprofiles,queryprofilesize,ancestors,subclassset):
	out=open("../results/Decay/Groupwise_SimJ_ProfileSize"+str(queryprofilesize)+"_Results.tsv",'w')

	out.write("QueryID\tNumber of annotations replaced\tBest Match\tBest Match Similarity\n")

	decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
	for queryid in decayedprofiles:
		for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
			groupwise_simj_results=[]
			queryprofile=decayedprofiles[queryid][str(numreplaced)]
			for dbprofileid in dbprofiles:
				dbprofile=dbprofiles[dbprofileid]
				simj,union=groupwise_jaccard(queryprofile,dbprofile,ancestors)
				subclassset=set.union(subclassset,union)
				groupwise_simj_results.append((simj,dbprofileid))

			bestmatch_groupwise_simj=max(groupwise_simj_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_groupwise_simj=max(groupwise_simj_results,key=itemgetter(0))[0]
			out.write(queryid+"\t"+str(numreplaced)+"\t"+bestmatch_groupwise_simj+"\t"+str(bestmatchsimilarity_groupwise_simj)+"\n")
			print queryid,numreplaced
	out.close()	
	return subclassset				



def compute_ind_profileic(profiles):	
	annotationset=set()
	icdict=dict()
	frequency=dict()
	annotationandanclist=dict()
	for profileid in profiles:
		for annotation in set(profiles[profileid]):
			if annotation not in frequency:
				frequency[annotation]=0
			frequency[annotation]+=1

	corpussize=len(profiles)
	maxic=round(-math.log(1/corpussize),2)
	for term in frequency:
		freq=frequency[term]
		if freq > corpussize:
			print "Frequency is greater than corpus size"
		ic=round((-math.log(freq/corpussize))/maxic,2)
		icdict[term]=ic
	
	return icdict


def compute_ind_ic(profiles):
	annotationcorpus=open("../Independent_AnnotationCorpus.txt",'w')
	annotationlist=[]
	icdict=dict()
	frequency=dict()
	annotationandanclist=dict()
	for profileid in profiles:
		for annotation in profiles[profileid]:
			annotationcorpus.write(annotation+"\n")
			annotationlist.append(annotation)
			if annotation not in frequency:
				frequency[annotation]=0
			frequency[annotation]+=1
			

	corpussize=len(annotationlist)
	maxic=round(-math.log(1/corpussize),2)
	for term in frequency:
		freq=frequency[term]
		ic=round((-math.log(freq/corpussize))/maxic,2)
		icdict[term]=ic
	annotationcorpus.close()
	
	return icdict


def sum_ic(termset,indannic,indprofileic,subclasses):
	print "in sum ic"
	nonredtermset=set()
	annicsum=0
	profileicsum=0
	for term in termset:
		if term in subclasses:
			nonredtermset=set.union(nonredtermset,subclasses[term])

	for term in nonredtermset:
		if term in indannic:
			ic=indannic[term]
		else:
			ic=0
		annicsum+=ic

		if term in indprofileic:
			ic=indprofileic[term]
		else:
			ic=0
		profileicsum+=ic
	print "out of sum_ic"
	return annicsum,profileicsum






def groupwise_related_resnik(profile1,profile2,icdict,profileicdict,ancestors):
	expprofile1=set()
	expprofile2=set()
	for term in profile1:
		expprofile1=set.union(expprofile1,ancestors[term])
	for term in profile2:
		expprofile2=set.union(expprofile2,ancestors[term])
	
	common=set.intersection(expprofile1,expprofile2)
	union=set.union(expprofile1,expprofile2)
	
	commonannic=np.sum([icdict[term] for term in common])
	commonprofileic=np.sum([profileicdict[term] for term in common])

	unionannic=np.sum([icdict[term] for term in union])
	unionprofileic=np.sum([profileicdict[term] for term in union])
	



	annresnik=commonannic/unionannic
	profileresnik=commonprofileic/unionprofileic
	
	return annresnik,profileresnik
	

def experiment_groupwise_related_resnik(dbprofiles,queryprofilesize,ancestors,icdict,profileicdict):
	out=open("../results/Decay/Groupwise_Related_Resnik_ProfileSize"+str(queryprofilesize)+"_Results.tsv",'w')

	out.write("QueryID\tNumber of annotations replaced\tBest Match Groupwise Annotation Resnik\tBest Match Similarity Groupwise Annotation Resnik\t Best Match Groupwise Profile Resnik\tBest Match Similarity Groupwise Profile Resnik\n")



	decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
	for queryid in decayedprofiles:
		for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
			groupwise_annresnik_results=[]
			groupwise_profileresnik_results=[]
			queryprofile=decayedprofiles[queryid][str(numreplaced)]
			for dbprofileid in dbprofiles:
				dbprofile=dbprofiles[dbprofileid]
				annresnik,profileresnik=groupwise_related_resnik(queryprofile,dbprofile,icdict,profileicdict,ancestors)

				groupwise_annresnik_results.append((annresnik,dbprofileid))
				groupwise_profileresnik_results.append((profileresnik,dbprofileid))
			bestmatch_groupwise_annresnik=max(groupwise_annresnik_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_groupwise_annresnik=max(groupwise_annresnik_results,key=itemgetter(0))[0]

			bestmatch_groupwise_profileresnik=max(groupwise_profileresnik_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_groupwise_profileresnik=max(groupwise_profileresnik_results,key=itemgetter(0))[0]


			out.write(queryid+"\t"+str(numreplaced)+"\t"+bestmatch_groupwise_annresnik+"\t"+str(bestmatchsimilarity_groupwise_annresnik)+"\t"+bestmatch_groupwise_profileresnik+"\t"+str(bestmatchsimilarity_groupwise_profileresnik)+"\n")
			print queryid,numreplaced
	out.close()









def groupwise_resnik(profile1,profile2,indannic,indprofileic,subclasses,ancestors):
	print "In grpwise"
	expprofile1=set()
	expprofile2=set()
	for term in profile1:
		expprofile1=set.union(expprofile1,ancestors[term])
	for term in profile2:
		expprofile2=set.union(expprofile2,ancestors[term])
	
	common=set.intersection(expprofile1,expprofile2)
	union=set.union(expprofile1,expprofile2)
	commonannic,commonprofileic=sum_ic(common,indannic,indprofileic,subclasses)
	unionannic,unionprofileic=sum_ic(union,indannic,indprofileic,subclasses)
	annresnik=commonannic/unionannic
	profileresnik=commonprofileic/unionprofileic
	print "out of grpwise"
	return annresnik,profileresnik
	

def experiment_groupwise_resnik(dbprofiles,queryprofilesize,ancestors,indannic,indprofileic,subclasses):
	print "In experiment"
	out=open("../results/Decay/Groupwise_Resnik_ProfileSize"+str(queryprofilesize)+"_Results.tsv",'w')

	out.write("QueryID\tNumber of annotations replaced\tBest Match Groupwise Annotation Resnik\tBest Match Similarity Groupwise Annotation Resnik\t Best Match Groupwise Profile Resnik\tBest Match Similarity Groupwise Profile Resnik\n")



	decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
	for queryid in decayedprofiles:
		for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
			groupwise_annresnik_results=[]
			groupwise_profileresnik_results=[]
			queryprofile=decayedprofiles[queryid][str(numreplaced)]
			for dbprofileid in dbprofiles:
				dbprofile=dbprofiles[dbprofileid]
				annresnik,profileresnik=groupwise_resnik(queryprofile,dbprofile,indannic,indprofileic, subclasses,ancestors)

				groupwise_annresnik_results.append((annresnik,dbprofileid))
				groupwise_profileresnik_results.append((profileresnik,dbprofileid))
			bestmatch_groupwise_annresnik=max(groupwise_annresnik_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_groupwise_annresnik=max(groupwise_annresnik_results,key=itemgetter(0))[0]

			bestmatch_groupwise_profileresnik=max(groupwise_profileresnik_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_groupwise_profileresnik=max(groupwise_profileresnik_results,key=itemgetter(0))[0]


			out.write(queryid+"\t"+str(numreplaced)+"\t"+bestmatch_groupwise_annresnik+"\t"+str(bestmatchsimilarity_groupwise_annresnik)+"\t"+bestmatch_groupwise_profileresnik+"\t"+str(bestmatchsimilarity_groupwise_profileresnik)+"\n")
			print queryid,numreplaced
	out.close()					


def experiment(dbprofiles,ancestors,numqueryprofiles,queryprofilesize,annotationpool,profileicdict,annotationicdict):
	bpsimilaritydict_profile=dict()
	bpsimilaritydict_annotation=dict()
	apsimilaritydict_profile=dict()
	apsimilaritydict_annotation=dict()


	profileids=set()
	for profileid in dbprofiles:
		if len(dbprofiles[profileid])==queryprofilesize:
			profileids.add(profileid)


	queryprofileids=random.sample(profileids, numqueryprofiles)
	filename="../results/Decay/ProfileBased_ProfileSize"+str(queryprofilesize)+"_Results.tsv"
	out=open("../results/Decay/Integrated_ProfileSize_SimJ"+str(queryprofilesize)+"_Results.tsv",'w')
	outfile=open("../results/Decay/AllScores_SimJ"+str(queryprofilesize)+"_Results.tsv",'w')
	#out.write("QueryID\tNumber of annotations replaced\tBest Match--BP Symmetric Profile IC\tBest Match Similarity--BP Symmetric Profile IC\tBest Match--BP Symmetric Annotation IC\tBest Match Similarity--BP Symmetric Annotation IC\tBest Match--BP Asymmetric Profile IC\tBest Match Similarity--BP Asymmetric Profile IC\tBest Match--BP Asymmetric Annotation IC\tBest Match Similarity--BP Asymmetric Annotation IC\tBest Match--AP Profile IC\tBest Match Similarity--AP Profile IC\tBest Match--AP Annotation IC\tBest Match Similarity--AP Annotation IC\tBest Match--BP Asymmetric Simj\tBest Match Similarity--BP Asymmetric Simj\tBest Match--BP Symmetric Simj\tBest Match Similarity--BP Symmetric Simj\tBest Match--AP Simj\tBest Match Similarity--AP Simj\n")
	out.write("QueryID\tNumber of annotations replaced\tBest Match--BP Asymmetric Simj\tBest Match Similarity--BP Asymmetric Simj\tBest Match--BP Symmetric Simj\tBest Match Similarity--BP Symmetric Simj\tBest Match--AP Simj\tBest Match Similarity--AP Simj\n")





		
	for queryid in queryprofileids:
		print queryid

		queryprofile=deepcopy(dbprofiles[queryid])
		replaced=set()
		for i in range(0,queryprofilesize+1):
			bp_asym_pic_results=[]
			bp_asym_aic_results=[]
			bp_sym_pic_results=[]
			bp_sym_aic_results=[]
			ap_pic_results=[]
			ap_aic_results=[]

			bp_asym_simj_results=[]
			bp_sym_simj_results=[]
			ap_simj_results=[]

			for dbprofileid in dbprofiles:


				profile2=dbprofiles[dbprofileid]

				# bp_sym_mediansim_pic,bpsimilaritydict_profile=calculate_bestpairs_symmetric(queryprofile,profile2,profileicdict,ancestors,bpsimilaritydict_profile)
				# bp_sym_pic_results.append((bp_sym_mediansim_pic,dbprofileid))
				
				# bp_asym_mediansim_pic,matchdata,bpsimilaritydict_profile,bestmatchiclist=calculate_bestpairs_asymmetric(queryprofile,profile2,profileicdict,ancestors,bpsimilaritydict_profile)
				# bp_asym_pic_results.append((bp_asym_mediansim_pic,dbprofileid))

				# bp_sym_mediansim_aic,bpsimilaritydict_annotation=calculate_bestpairs_symmetric(queryprofile,profile2,annotationicdict,ancestors,bpsimilaritydict_annotation)
				# bp_sym_aic_results.append((bp_sym_mediansim_aic,dbprofileid))

				# bp_asym_mediansim_aic,matchdata,bpsimilaritydict_annotation,bestmatchiclist=calculate_bestpairs_asymmetric(queryprofile,profile2,annotationicdict,ancestors,bpsimilaritydict_annotation)
				# bp_asym_aic_results.append((bp_asym_mediansim_aic,dbprofileid))
				

				# ap_mediansim_pic,apsimilaritydict_profile=get_allpairs_medianic(queryprofile,profile2,profileicdict,ancestors,bpsimilaritydict_profile,apsimilaritydict_profile)
				# ap_pic_results.append((ap_mediansim_pic,dbprofileid))


				# ap_mediansim_aic,apsimilaritydict_annotation=get_allpairs_medianic(queryprofile,profile2,annotationicdict,ancestors,bpsimilaritydict_annotation,apsimilaritydict_annotation)
				# ap_aic_results.append((ap_mediansim_aic,dbprofileid))

				bp_asym_simj=calculate_bestpairs_jaccard_asymmetric(queryprofile,profile2,ancestors)
				bp_asym_simj_results.append((bp_asym_simj,dbprofileid))

				bp_sym_simj=calculate_bestpairs_jaccard_symmetric(queryprofile,profile2,ancestors)
				bp_sym_simj_results.append((bp_sym_simj,dbprofileid))

				ap_simj=calculate_allpairs_jaccard(queryprofile,profile2,ancestors)
				ap_simj_results.append((ap_simj,dbprofileid))
		


			
			# bestmatch_bp_sym_pic=max(bp_sym_pic_results,key=itemgetter(0))[1]	
			# bestmatchsimilarity_bp_sym_pic=max(bp_sym_pic_results,key=itemgetter(0))[0]
			# allresults=[y[0] for y in bp_sym_pic_results]
			# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Sym_PIC\t"+",".join(str(x)for x in allresults)+"\n")
			
			# bestmatch_bp_sym_aic=max(bp_sym_aic_results,key=itemgetter(0))[1]	
			# bestmatchsimilarity_bp_sym_aic=max(bp_sym_aic_results,key=itemgetter(0))[0]
			# allresults=[y[0] for y in bp_sym_aic_results]
			# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Sym_AIC\t"+",".join(str(x)for x in allresults)+"\n")


			# bestmatch_bp_asym_pic=max(bp_asym_pic_results,key=itemgetter(0))[1]	
			# bestmatchsimilarity_bp_asym_pic=max(bp_asym_pic_results,key=itemgetter(0))[0]
			# allresults=[y[0] for y in bp_asym_pic_results]
			# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Asym_PIC\t"+",".join(str(x)for x in allresults)+"\n")
			
			# bestmatch_bp_asym_aic=max(bp_asym_aic_results,key=itemgetter(0))[1]	
			# bestmatchsimilarity_bp_asym_aic=max(bp_asym_aic_results,key=itemgetter(0))[0]
			# allresults=[y[0] for y in bp_asym_aic_results]
			# outfile.write(queryid+"\t"+str(i)+"\t"+"BP_Asym_AIC\t"+",".join(str(x)for x in allresults)+"\n")
				
			# bestmatch_ap_pic=max(ap_pic_results,key=itemgetter(0))[1]	
			# bestmatchsimilarity_ap_pic=max(ap_pic_results,key=itemgetter(0))[0]
			# allresults=[y[0] for y in ap_pic_results]
			# outfile.write(queryid+"\t"+str(i)+"\t"+"AP_PIC\t"+",".join(str(x)for x in allresults)+"\n")

			# bestmatch_ap_aic=max(ap_aic_results,key=itemgetter(0))[1]	
			# bestmatchsimilarity_ap_aic=max(ap_aic_results,key=itemgetter(0))[0]
			# allresults=[y[0] for y in ap_aic_results]
			# outfile.write(queryid+"\t"+str(i)+"\t"+"AP_AIC\t"+",".join(str(x)for x in allresults)+"\n")
			
			
			bestmatch_bp_asym_simj=max(bp_asym_simj_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_bp_asym_simj=max(bp_asym_simj_results,key=itemgetter(0))[0]
			
			bestmatch_bp_sym_simj=max(bp_sym_simj_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_bp_sym_simj=max(bp_sym_simj_results,key=itemgetter(0))[0]
			
			bestmatch_ap_simj=max(ap_simj_results,key=itemgetter(0))[1]	
			bestmatchsimilarity_ap_simj=max(ap_simj_results,key=itemgetter(0))[0]
			

			out.write(queryid+"\t"+str(i)+"\t"+bestmatch_bp_asym_simj+"\t"+str(bestmatchsimilarity_bp_asym_simj)  +"\t"+ bestmatch_bp_sym_simj+"\t"+str(bestmatchsimilarity_bp_sym_simj)  +"\t"+ bestmatch_ap_simj+"\t"+str(bestmatchsimilarity_ap_simj)+   "\n")


			if len(replaced)<len(queryprofile):	
				queryprofile,replaced=substitute_annotation(queryprofile,replaced,annotationpool)
	out.close()




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



def load_groupwise_resnik_decay(infile,signalbestscores):
	numreplacedset=set()
	infile.next()
	for line in infile:
		queryid,numreplaced,bestmatch_aic,score_aic,bestmatch_pic,score_pic=line.split("\t")
		numreplaced=int(numreplaced)
	

		if numreplaced not in signalbestscores:
			signalbestscores[numreplaced]=dict()
			numreplacedset.add(numreplaced)

		if 'groupwise_aic_resnik' not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced]['groupwise_aic_resnik']=[]
		signalbestscores[numreplaced]['groupwise_aic_resnik'].append(float(score_aic))

		if 'groupwise_pic_resnik' not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced]['groupwise_pic_resnik']=[]
		signalbestscores[numreplaced]['groupwise_pic_resnik'].append(float(score_pic))


	return signalbestscores



def load_groupwise_simj_decay(infile,metric,signalbestscores):
	numreplacedset=set()
	infile.next()
	for line in infile:
		queryid,numreplaced,bestmatch,score=line.split("\t")
		numreplaced=int(numreplaced)
	

		if numreplaced not in signalbestscores:
			signalbestscores[numreplaced]=dict()
			numreplacedset.add(numreplaced)

		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(score))

	return signalbestscores

def load_decay_results(infile):
	metricset=set()
	numreplacedset=set()
	infile.next()
	signalbestscores=dict()
	for line in infile:
		queryid,numreplaced,match_bp_sym_pic,bp_sym_pic,match_bp_sym_aic,bp_sym_aic,match_bp_asym_pic,bp_asym_pic,match_bp_asym_aic,bp_asym_aic,match_ap_pic,ap_pic,match_ap_aic,ap_aic,match_bp_asym_simj,bp_asym_simj,match_bp_sym_simj,bp_sym_simj,match_ap_simj,ap_simj=line.split("\t")
		numreplaced=int(numreplaced)
	

		if numreplaced not in signalbestscores:
			signalbestscores[numreplaced]=dict()
			numreplacedset.add(numreplaced)


		metric="bp_sym_pic"
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_sym_pic))

		metric="bp_sym_aic"
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_sym_aic))


		metric="bp_asym_pic"
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_asym_pic))

		metric="bp_asym_aic"
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(bp_asym_aic))


		metric="ap_pic"
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(ap_pic))

		metric="ap_aic"
		metricset.add(metric)
		if metric not in signalbestscores[numreplaced]:
			signalbestscores[numreplaced][metric]=[]
		signalbestscores[numreplaced][metric].append(float(ap_aic))

		metric="bp_sym_simj"
		metricset.add(metric)
		if bp_sym_simj.strip() !="":
			if metric not in signalbestscores[numreplaced]:
				signalbestscores[numreplaced][metric]=[]
			
			signalbestscores[numreplaced][metric].append(float(bp_sym_simj))

		metric="bp_asym_simj"
		metricset.add(metric)
		if bp_asym_simj.strip() !="":
			if metric not in signalbestscores[numreplaced]:
				signalbestscores[numreplaced][metric]=[]
			signalbestscores[numreplaced][metric].append(float(bp_asym_simj))

		metric="ap_simj"
		metricset.add(metric)
		if ap_simj.strip() !="":
			if metric not in signalbestscores[numreplaced]:
				signalbestscores[numreplaced][metric]=[]
			signalbestscores[numreplaced][metric].append(float(ap_simj))

	return metricset,numreplacedset,signalbestscores
def plot_decay_allsizes():
	f, axarr = plt.subplots(3, 3)
	colors=['black','red','green']
	legend=[]
	colorindex=0


	title=dict()
	title['ap_aic']="Annotation IC"
	title['ap_pic']="All Pairs \n\n Profile IC"
	title['ap_simj']="Jaccard"
	title['bp_asym_aic']="Annotation IC"
	title['bp_asym_pic']="Best Pairs Asymmetric\n\nProfile IC"
	title['bp_asym_simj']="Jaccard"
	title['bp_sym_aic']="Annotation IC"
	title['bp_sym_pic']="Best Pairs Symmetric\n\nProfile IC"
	title['bp_sym_simj']="Jaccard"
	infile=open("../results/Decay/Integrated_ProfileSize10_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	print numreplacedset
	i=j=0
	for metric in sorted(metricset):
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		#axarr[i, j].set_xlim(0,max(numreplacedset))
		
		j+=1
		if j==3:
			i+=1
			j=0
		if i==3:
			i=0	
	


	infile=open("../results/Decay/Integrated_ProfileSize20_Results.tsv")
	colorindex+=1
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	print numreplacedset
	i=j=0
	for metric in sorted(metricset):
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		#axarr[i, j].set_xlim(0,max(numreplacedset))
		
		j+=1
		if j==3:
			i+=1
			j=0
		if i==3:
			i=0	
	


	infile=open("../results/Decay/Integrated_ProfileSize40_Results.tsv")
	colorindex+=1
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	print numreplacedset
	i=j=0
	for metric in sorted(metricset):
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		#axarr[i, j].set_xlim(0,max(numreplacedset))
		
		j+=1
		if j==3:
			i+=1
			j=0
		if i==3:
			i=0	



	plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
	plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	#plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)	
	#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)	
	axarr[2,1].set_xlabel('Number of annotations replaced')
	axarr[1,0].set_ylabel('Similarity score')	
	#plt.show()
	plt.tight_layout()
	plt.savefig("../results/Decay_Nov.png")
	
def compare_simj():
	title=dict()
	title['ap_simj']="All Pairs Jaccard"
	title['bp_asym_simj']="Best Pairs Asymmetric\nJaccard"
	title['bp_sym_simj']="Best Pairs Symmetric\nJaccard"
	title['groupwise_jaccard']="Groupwise Jaccard"
	f, axarr = plt.subplots(2, 2)
	colors=['black','red','green']
	legend=[]
	colorindex=0

	infile=open("../results/Decay/Integrated_ProfileSize10_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	infile.close()
	infile=open("../results/Decay/Groupwise_SimJ_ProfileSize10_Results.tsv")
	signalbestscores=load_groupwise_simj_decay(infile,"groupwise_jaccard",signalbestscores)
	
	i=j=0
	for metric in ['ap_simj','bp_asym_simj','bp_sym_simj','groupwise_jaccard']:
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		axarr[i, j].set_xlim(0,40)
		j+=1
		if j==2:
			i+=1
			j=0
		if i==2:
			i=0


	infile=open("../results/Decay/Integrated_ProfileSize20_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	infile.close()
	infile=open("../results/Decay/Groupwise_SimJ_ProfileSize20_Results.tsv")
	signalbestscores=load_groupwise_simj_decay(infile,"groupwise_jaccard",signalbestscores)
	colorindex+=1
	i=j=0
	for metric in ['ap_simj','bp_asym_simj','bp_sym_simj','groupwise_jaccard']:
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		axarr[i, j].set_xlim(0,40)
		j+=1
		if j==2:
			i+=1
			j=0
		if i==2:
			i=0

	infile=open("../results/Decay/Integrated_ProfileSize40_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	infile.close()
	infile=open("../results/Decay/Groupwise_SimJ_ProfileSize40_Results.tsv")
	signalbestscores=load_groupwise_simj_decay(infile,"groupwise_jaccard",signalbestscores)
	
	colorindex+=1
	i=j=0
	for metric in ['ap_simj','bp_asym_simj','bp_sym_simj','groupwise_jaccard']:
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		axarr[i, j].set_xlim(0,40)
		print i,j
		j+=1
		if j==2:
			i+=1
			j=0
		if i==2:
			i=0




	#plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	#plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	#plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)	
	#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)	
	axarr[1,1].set_xlabel('Number of annotations replaced')
	axarr[1,0].set_xlabel('Number of annotations replaced')
	axarr[0,0].set_ylabel('Similarity score')
	axarr[1,0].set_ylabel('Similarity score')
	plt.tight_layout()
	plt.savefig("../results/JaccardComparison.eps",dpi=1200)


def compare_resnik():
	title=dict()
	title['groupwise_aic_resnik']="Groupwise Resnik Annotation IC"
	title['groupwise_pic_resnik']="Groupwise Resnik Profile IC"
	title['ap_aic']="All pairs Annotation IC"
	title['ap_pic']="All Pairs Profile IC"
	title['bp_asym_aic']="Best Pairs Asymmetric Annotation IC"
	title['bp_asym_pic']="Best Pairs Asymmetric Profile IC"
	title['bp_sym_aic']="Best Pairs Symmetric Annotation IC"
	title['bp_sym_pic']="Best Pairs Symmetric Profile IC"
	
	f, axarr = plt.subplots(4, 2)
	colors=['black','red','green']
	legend=[]
	colorindex=0

	infile=open("../results/Decay/Integrated_ProfileSize10_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	infile.close()
	infile=open("../results/Decay/Groupwise_Related_Resnik_ProfileSize10_Results.tsv")
	signalbestscores=load_groupwise_resnik_decay(infile,signalbestscores)
	infile.close()
	i=j=0
	for metric in ['ap_aic','ap_pic','bp_sym_aic','bp_sym_pic','bp_asym_aic','bp_asym_pic','groupwise_aic_resnik','groupwise_pic_resnik']:
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		axarr[i, j].set_xlim(0,40)
		j+=1
		if j==2:
			i+=1
			j=0
		if i==4:
			i=0


	infile=open("../results/Decay/Integrated_ProfileSize20_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	infile.close()
	infile=open("../results/Decay/Groupwise_Related_Resnik_ProfileSize20_Results.tsv")
	signalbestscores=load_groupwise_resnik_decay(infile,signalbestscores)
	infile.close()
	colorindex+=1
	i=j=0
	for metric in ['ap_aic','ap_pic','bp_sym_aic','bp_sym_pic','bp_asym_aic','bp_asym_pic','groupwise_aic_resnik','groupwise_pic_resnik']:
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		axarr[i, j].set_xlim(0,40)
		j+=1
		if j==2:
			i+=1
			j=0
		if i==4:
			i=0

	infile=open("../results/Decay/Integrated_ProfileSize40_Results.tsv")
	metricset,numreplacedset,signalbestscores=load_decay_results(infile)
	infile.close()
	infile=open("../results/Decay/Groupwise_Related_Resnik_ProfileSize40_Results.tsv")
	signalbestscores=load_groupwise_resnik_decay(infile,signalbestscores)
	infile.close()
	print signalbestscores
	colorindex+=1
	i=j=0
	for metric in ['ap_aic','ap_pic','bp_sym_aic','bp_sym_pic','bp_asym_aic','bp_asym_pic','groupwise_aic_resnik','groupwise_pic_resnik']:
		signallist=[]
		errorlist=[]
		for numreplaced in sorted(numreplacedset):
			signallist.append(np.mean(signalbestscores[numreplaced][metric]))
			errorlist.append(error(signalbestscores[numreplaced][metric]))
		axarr[i, j].errorbar(list(sorted(numreplacedset)),signallist,yerr=errorlist,color=colors[colorindex])
		axarr[i, j].set_title(title[metric])
		axarr[i, j].set_ylim(0,1)
		axarr[i, j].set_xlim(0,40)
		print i,j
		j+=1
		if j==2:
			i+=1
			j=0
		if i==4:
			i=0




	#plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	#plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	#plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)	
	#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)	
	axarr[3,1].set_xlabel('Number of annotations replaced')
	axarr[3,0].set_xlabel('Number of annotations replaced')
	axarr[0,0].set_ylabel('Similarity score')
	axarr[1,0].set_ylabel('Similarity score')
	axarr[2,0].set_ylabel('Similarity score')
	axarr[3,0].set_ylabel('Similarity score')
	f.set_tight_layout(True)
	plt.savefig("../results/ResnikComparison.eps",dpi=1200)

def plot_expect_difference():
	i=j=0
	f, axarr = plt.subplots(3, 3)
	metrics=['bp_sym_pic','bp_sym_aic','bp_asym_pic','bp_asym_aic','ap_pic','ap_aic','bp_asym_simj','bp_sym_simj','ap_simj']
	for metric in metrics:
		df=pd.read_csv('../results/EmpiricalvsRegression_'+metric+'.tsv', sep='\t')
		empirical=np.log10(df['Empirical p-value'].tolist())
		regression=np.log10(df['Regression p-value'].tolist())
		axarr[i,j].scatter(empirical,regression,color='black',alpha=0.2)
		axarr[i, j].set_title(metric)
		axarr[i, j].set_xlim(-5,0)
		axarr[i, j].set_ylim(-5,0)
		j+=1
		if j==3:
			i+=1
			j=0
		if i==3:
			i=0	

	plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
	plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)	
	plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)	
	axarr[2,1].set_xlabel('log10(Empirical p-value)')
	axarr[1,0].set_ylabel('log10(Regression p-value)')	

	plt.savefig("../results/pvalueDifference.png", dpi=1200)


def main():
	#plot_expect_difference()
	#sys.exit()
	#plot_decay_allsizes()
	#sys.exit()
	
	################################
	# Compare Groupwise metrics to regular metrics
	# queryprofilesize=int(sys.argv[1])
	# numqueryprofiles=int(sys.argv[2])
	# numofdbprofiles=659
	
	# compare_simj()
	# compare_resnik()
	# sys.exit()

	################################
	
	

	################################
	# Run AP and BP HRSS
	queryprofilesize=int(sys.argv[1])
	dbprofiles,annotationpool=load_randomprofiles()
	decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
	ancestors=load_ancestors()
	print "loaded ancestors"
	profileterms=get_profile_terms(decayedprofiles,dbprofiles)
	subclasses=load_subclasses(profileterms)
	print "loaded subclasses"
	annotationicdict=json.load(open("../data/AnnotationIC.txt"))
	experiment_hrss(decayedprofiles,dbprofiles,queryprofilesize,ancestors,"Annotation IC",annotationicdict,subclasses)
	sys.exit()
	################################




	# noisefile="../results/Noise_Size10_BestMatches.tsv"
	# signalfile=open("../results/Decay/Integrated_ProfileSize10_Results.tsv")
	# signal_vs_noise(signalfile,noisefile)
	# sys.exit()

	# noisefile=open("../results/Decay/AllScores_NoisetoNoise_20_Results.tsv")
	# signalfile=open("../results/Decay/AllScores20_Results.tsv")
	# comparisonfile=open("../results/Decay/SignalvsNoise_20.tsv")
	
	#sys.exit()
	# compare_signal_noise(noisefile,signalfile,queryprofilesize)
	# plot_signal_noise(comparisonfile, queryprofilesize)
	# sys.exit()
	# plot_ic_comparison()
	# sys.exit()


	bpsimilaritydict=dict()
	dbprofiles,annotationpool=load_randomprofiles()
	#create_decayed_profiles(dbprofiles,numqueryprofiles,queryprofilesize,annotationpool)
	#sys.exit()


	print "Loaded randomprofiles"
	ancestors=load_ancestors() 
	print "Loaded ancestors"
	
	#indprofileic=compute_ind_profileic(dbprofiles)
	#indannic=compute_ind_ic(dbprofiles)
	#subclasses=load_subclasses()
	#print "Loaded subclasses"
	
	#experiment_groupwise_related_resnik(dbprofiles,20,ancestors,indannic,indprofileic,subclasses)
	#experiment_groupwise_related_resnik(dbprofiles,40,ancestors,indannic,indprofileic,subclasses)
	
	#subclassset=set()
	#subclassset=experiment_groupwise_jaccard(dbprofiles,10,ancestors,subclassset)
	#subclassset=experiment_groupwise_jaccard(dbprofiles,20,ancestors,subclassset)
	#subclassset=experiment_groupwise_jaccard(dbprofiles,40,ancestors,subclassset)
	#out=open("../data/SubclassesNeeded.txt",'w')
	#for term in subclassset:
	#	out.write(term+"\n")
	#out.close()
	#sys.exit()

	prepare_profile_corpus(dbprofiles,ancestors)
	infile=open("../data/Corpus_Profile_Frequency.txt")
	profileicdict,profilefrequency=compute_ic_profile_frequency(infile)
	infile.close()
	print "Computed profile based IC"

	

	annotationicdict=compute_ic(dbprofiles,ancestors)
	print "Computed annotation based IC"


	#experiment_groupwise_related_resnik(dbprofiles,10,ancestors,annotationicdict,profileicdict)
	experiment_groupwise_related_resnik(dbprofiles,20,ancestors,annotationicdict,profileicdict)
	#experiment_groupwise_related_resnik(dbprofiles,40,ancestors,annotationicdict,profileicdict)
	sys.exit()


	# resfile="../results/54_54_DecayResults.tsv"
	# noisefile="../results/NoiseResults_54_54.tsv"
	# reg_results=regression_allmetrics(dbprofiles)
	# compute_expect_allmetrics(resfile,noisefile,reg_results,dbprofiles,numofdbprofiles,queryprofilesize)
	# plot_expect_difference()
	# sys.exit()


	# plot_expect_difference()

	# sys.exit()
	# compare_ic(dbprofiles,ancestors,numqueryprofiles,numofdbprofiles,annotationpool)
	# sys.exit()

	




	
	
	#run_noise_to_noise(dbprofiles,ancestors,queryprofilesize,annotationpool,profileicdict,annotationicdict)
	#sys.exit()
	

	
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
	import json
	import pandas as pd
	import matplotlib.pyplot as plt
	from copy import deepcopy
	from operator import itemgetter
	from scipy.stats import rankdata
	import numpy as np
	import math
	import statsmodels.api as sm
	import statsmodels.stats.api as sms
	from scipy import stats
	from scipy.stats import boxcox
	from statsmodels.stats import stattools
	import sys
	main()
