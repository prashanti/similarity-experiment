from __future__ import division
def getmicaic(term1,term2,ancestors,icdict):
    micaic=0
    mica=""
    commonancestors=set.intersection(ancestors[term1],ancestors[term2])
    lcslist=[icdict[anc] for anc in commonancestors]
    
    
    if len(lcslist)>0:
        micaic=np.max(lcslist)     
        return micaic
    else:
        return 0

def load_ancestors(granularity):
    ancestors=dict()
    if granularity =='E':
        infile=open("../data/ESubsumers.txt")
    elif granularity=='EA':
        infile=open("../data/EASubsumers.txt")
    elif granularity =='EQ':
        infile=open("../data/EQSubsumers.txt")
    for line in infile:
        term,subsumer=line.strip().split("\t")
        if term not in ancestors:
            ancestors[term]=set()
            ancestors[term].add(term)
        if subsumer !="owl:Thing":      
            ancestors[term].add(subsumer)
    infile.close()
    return ancestors

def load_randomprofiles(granularity):
    randomprofiles=dict()
    annotationpool=[]
    if granularity=='E':
        infile=open("../data/RandomProfiles2016_E.txt")
    else:
        infile=open("../data/RandomProfiles2016.txt")
    for line in infile:
        profileid,annotation=line.strip().split("\t")
        if profileid not in randomprofiles:
            randomprofiles[profileid]=[]
        randomprofiles[profileid].append(annotation)
        annotationpool.append(annotation)
    infile.close()
    return randomprofiles,annotationpool

def calculate_bestpairs_symmetric_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
    finalsim=0
    bestmatchiclist1=[]
    bestmatchiclist2=[]
    termmatchic=[]
    matchdata=[]
    for term1 in profile1:
        termmatchic=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in bpsimilaritydict:
                termmatchic.append(bpsimilaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                termmatchic.append(micaic)
                bpsimilaritydict[termtuple]=micaic
        
        bestmatchiclist1.append(np.max(termmatchic))
    

    termmatchic=[]
    matchdata=[]
    for term1 in profile2:
        termmatchic=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            termmatchic.append(bpsimilaritydict[termtuple])
        bestmatchiclist2.append(np.max(termmatchic))
            
    return bestmatchiclist1,bestmatchiclist2,bpsimilaritydict


def calculate_bestpairs_asymmetric_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
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
                micaic=getmicaic(term1,term2,ancestors,icdict)
                termmatchic.append(micaic)
                bpsimilaritydict[termtuple]=micaic
        bestmatchiclist.append(np.max(termmatchic))
    return bpsimilaritydict,bestmatchiclist

def calculate_bestpairs_asymmetric_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    bestmatchsimj=[]
    for term1 in profile1:
        termmatchsimj=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                simj=similaritydict[termtuple]
                termmatchsimj.append(simj)
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
            termmatchsimj.append(simj)
        bestmatchsimj.append(max(termmatchsimj))
    return similaritydict,bestmatchsimj
    
    



def calculate_bestpairs_symmetric_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    bestmatchsimj1=[]
    bestmatchsimj2=[]
    termmatchsimj=[]
    for term1 in profile1:
        termmatchsimj=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                simj=similaritydict[termtuple]
                termmatchsimj.append(simj)
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
            termmatchsimj.append(simj)
        bestmatchsimj1.append(max(termmatchsimj))
    
    
    termmatchsimj=[]
    for term1 in profile2:
        termmatchsimj=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            simj=similaritydict[termtuple]
            termmatchsimj.append(simj)
        bestmatchsimj2.append(max(termmatchsimj))
    
    
    return bestmatchsimj1,bestmatchsimj2,similaritydict



def getsimj(term1,term2,ancestors):
    if len(set.union(ancestors[term1],ancestors[term2])) >0:
        simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
    else:
        simj=0
    return simj



def experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,icflag,metric,granularity):
    bpsimilaritydict=dict()
    decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
    out=open("../results/FullDistribution/AnnotationReplacement/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPAsym_"+icflag+ "_"+metric+"_Results.tsv",'w')
    out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List\n")

    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    if metric == "Resnik":
                        bpsimilaritydict,bestmatchlist=calculate_bestpairs_asymmetric_resnik(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)
                    elif metric == "Jaccard":
                        bpsimilaritydict,bestmatchlist=calculate_bestpairs_asymmetric_jaccard(queryprofile,profile2,ancestors,bpsimilaritydict)
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in bestmatchlist)+"\n")
    out.close()

def experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,icflag,metric,granularity):
    bpsimilaritydict=dict()
    decayedprofiles=json.load(open("../data/DecayedProfilesDict_Size"+str(queryprofilesize)+".txt"))
    
    out=open("../results/FullDistribution/AnnotationReplacement/"+granularity+"_ProfileSize"+str(queryprofilesize)+"_BPSym_"+icflag+ "_"+metric+"_Results.tsv",'w')

    out.write("Number of annotations replaced\tQuery ID\tDatabase ID\tScore List1\tScore list2\n")

    for queryid in decayedprofiles:
            for numreplaced in sorted([int(x) for x in decayedprofiles[queryid]]):
                print queryid,numreplaced
                queryprofile=decayedprofiles[queryid][str(numreplaced)]        
                for dbprofileid in dbprofiles:
                    profile2=dbprofiles[dbprofileid]
                    if metric =="Resnik":
                        bestmatchlist1,bestmatchlist2,bpsimilaritydict=calculate_bestpairs_symmetric_resnik(queryprofile,profile2,icdict,ancestors,bpsimilaritydict)
                    elif metric=="Jaccard":
                        bestmatchlist1,bestmatchlist2,bpsimilaritydict=calculate_bestpairs_symmetric_jaccard(queryprofile,profile2,ancestors,bpsimilaritydict)

                    
                    out.write(str(numreplaced)+"\t"+queryid+"\t"+dbprofileid+"\t"+','.join(str(x) for x in bestmatchlist1)+"\t"+','.join(str(x) for x in bestmatchlist2)+"\n")
    out.close()

def compute_annotation_ic(profiles,ancestors,granularity):
    
    outfile=open("../data/"+granularity+"_AnnotationIC.txt",'w')
    annotationlist=[]
    icdict=dict()
    frequency=dict()
    annotationandanclist=dict()
    for profileid in profiles:
        for annotation in profiles[profileid]:
            annotationlist.append(annotation)
            
            for anc in ancestors[annotation]:
                if anc not in frequency:
                    frequency[anc]=0
                frequency[anc]+=1


    corpussize=len(annotationlist)
    maxic=round(-math.log(1/corpussize),2)
    for term in frequency:
        freq=frequency[term]
        ic=round((-math.log(freq/corpussize))/maxic,2)
        icdict[term]=ic
    json.dump(icdict,outfile)
    outfile.close()
    
def run_bp_asym_resnik(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=load_ancestors(granularity)
    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        print "Beginning experiment"
        experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Resnik",granularity)
    sys.exit()

def run_bp_sym_resnik(icflag):
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=load_ancestors(granularity)
    if icflag=="PIC":
        icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
        print "Beginning experiment"
        experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,icdict,"PIC","Resnik",granularity)
    sys.exit()

def run_bp_sym_jaccard():
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=load_ancestors(granularity)
    print "Beginning experiment"
    experiment_bestpairs_symmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,dict(),"","Jaccard",granularity)
    sys.exit()

def run_bp_asym_jaccard():    
    queryprofilesize=int(sys.argv[1])
    granularity=sys.argv[2]
    dbprofiles,annotationpool=load_randomprofiles(granularity)
    ancestors=load_ancestors(granularity)
    print "Beginning experiment"
    experiment_bestpairs_asymmetric(dbprofiles,ancestors,queryprofilesize,annotationpool,dict(),"","Jaccard",granularity)
    
    sys.exit()
def precompute_annotation_ic():
    dbprofiles,annotationpool=load_randomprofiles('E')
    compute_annotation_ic(dbprofiles,ancestors,'E')
    
    dbprofiles,annotationpool=load_randomprofiles('EA')
    compute_annotation_ic(dbprofiles,ancestors,'EA')

    dbprofiles,annotationpool=load_randomprofiles('EQ')
    ancestors=load_ancestors('EQ')
    compute_annotation_ic(dbprofiles,ancestors,'EQ')


def prepare_profile_corpus(dbprofiles,ancestors):
    outfile=open("../data/Corpus_Profile_Frequency.txt",'w')
    for profileid in dbprofiles:
        profileancestors=set()
        for annotation in set(dbprofiles[profileid]):
            ancset=ancestors[annotation]
            profileancestors.add(annotation)
            profileancestors=set.union(profileancestors,ancset)
        outfile.write(profileid+"\t"+";".join(profileancestors)+"\n")
    outfile.close()


def compute_ic_profile_frequency(dbprofiles,ancestors,granularity):
    prepare_profile_corpus(dbprofiles,ancestors)
    infile=open("../data/Corpus_Profile_Frequency.txt")

    icdict=dict()
    frequency=dict()
    corpussize=0
    for line in infile:
        profileid,profileancestors=line.strip().split("\t")
        profileancestors=profileancestors.split(";")
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
    outfile=open("../data/"+granularity+"_ProfileIC.txt",'w')
    json.dump(icdict,outfile)    
    
def precompute_profile_ic():
    dbprofiles,annotationpool=load_randomprofiles('E')
    ancestors=load_ancestors('E')
    compute_ic_profile_frequency(dbprofiles,ancestors,'E')
    
    dbprofiles,annotationpool=load_randomprofiles('EA')
    ancestors=load_ancestors('EA')
    compute_ic_profile_frequency(dbprofiles,ancestors,'EA')

    dbprofiles,annotationpool=load_randomprofiles('EQ')
    ancestors=load_ancestors('EQ')
    compute_ic_profile_frequency(dbprofiles,ancestors,'EQ')

def main():
    precompute_profile_ic()
    precompute_ic()
    
    metriccombo=sys.argv[3] 
    if metriccombo=="bp_asym_pic_resnik":
        run_bp_asym_resnik("PIC")  
    elif metriccombo=="bp_sym_pic_resnik":
        run_bp_sym_resnik("PIC")
    elif metriccombo=="bp_asym_aic_resnik":
        run_bp_asym_resnik("AIC")  
    elif metriccombo=="bp_sym_aic_resnik":
        run_bp_sym_resnik("AIC")    
    elif metriccombo=="bp_sym_jaccard":
        run_bp_sym_jaccard()
    elif metriccombo=="bp_asym_jaccard":    
        run_bp_asym_jaccard()
    
    


if __name__ == "__main__":
    import os
    import json
    import pandas as pd
    from copy import deepcopy
    from operator import itemgetter
    from scipy.stats import rankdata
    import numpy as np
    import math
    from scipy import stats
    import sys
    main()    