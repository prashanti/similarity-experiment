from __future__ import division
def get_profiles():
	realprofiles=dict()
	realprofileindices=dict()
	annotationpool=[]
	out=open("../results/EvolutionaryProfileSizes.tsv",'w')
	register('text/rdf+n3', Parser, 'rdflib.plugins.parsers.notation3', 'N3Parser')
	g = rdflib.Graph()
	result = g.parse('../data/realprofilesold.ttl', format='n3')
	index=0
	for stmt in g:
		profileid=stmt[0]
		annotation=stmt[2]
		if "VTO" in profileid:
			profileid=profileid.replace("http://purl.obolibrary.org/obo/","").replace("#profile","")
			if profileid not in realprofiles:
				realprofiles[profileid]=set()
			realprofiles[profileid].add(annotation)
			annotationpool.append(annotation)

			if profileid not in realprofileindices:
				realprofileindices[profileid]=set()	
			realprofileindices[profileid].add(index)
			index+=1	
	print "Created real profiles"
	for profileid in realprofiles:
		out.write(profileid+"\t"+str(len(realprofiles[profileid]))+"\n")
	return realprofiles,annotationpool,realprofileindices

def create_random_profiles(realprofiles,annotationpool):
	randomprofiles=dict()
	randomindices=[]
	randomprofileindices=dict()

	selected=set()		
	poolsize=len(annotationpool)
	realindices=list(range(0,len(annotationpool)))
	for profileid in realprofiles:
		realprofilesize=len(realprofiles[profileid])
		if profileid not in randomprofiles:
			randomprofiles[profileid]=set()
		while (len(randomprofiles[profileid])<realprofilesize):
			randomindex=random.randint(0,poolsize-1)
			if (annotationpool[randomindex] not in randomprofiles[profileid]) and (randomindex not in selected):
				randomprofiles[profileid].add(annotationpool[randomindex])
				if profileid not in randomprofileindices:
					randomprofileindices[profileid]=set()
				randomprofileindices[profileid].add(randomindex)
				randomindices.append(randomindex)
				selected.add(randomindex)
	return randomprofiles,realindices,randomindices,randomprofileindices	

def checkrandomization(realindices,randomindices):
	coeff,p=scipy.stats.pearsonr(realindices, randomindices)
	
	#0.000562644964 0.893749763777
	
	plt.scatter(realindices[40:120],randomindices[40:120],c='blue')
	plt.scatter(realindices[40:120],realindices[40:120],c='red')
	plt.show()
	return coeff,p


def checkoverlap(realprofileindices,randomprofileindices):
	
	realcooccurrence=dict()
	randomcooccurrence=dict()
	for profileid in realprofileindices:
		for index in realprofileindices[profileid]:
			if index not in realcooccurrence:
				realcooccurrence[index]=set()

			remaining=set(realprofileindices[profileid]).difference({index})

			realcooccurrence[index]=remaining


	for profileid in randomprofileindices:
		for index in randomprofileindices[profileid]:
			if index not in randomcooccurrence:
				randomcooccurrence[index]=set()
			remaining=set(randomprofileindices[profileid]).difference({index})
			randomcooccurrence[index]=remaining
                 



	overlaplist=[]
	for index in realcooccurrence:
		realset=realcooccurrence[index]
		randomset=randomcooccurrence[index]
		intersection=set.intersection(realset,randomset)
		union=set.union(realset,randomset)
		overlap=len(intersection)/len(union)
		overlaplist.append(overlap)

	return np.median(overlaplist)	


def main():
	realprofiles,annotationpool,realprofileindices=get_profiles()
	
	coefflist=[]
	plist=[]
	# out=open("../results/AnnotationSetOld.txt",'w')
	# for annotation in set(annotationpool):
	# 	out.write(annotation+"\n")
	# out.close()	
	randomprofiles,realindices,randomindices,randomprofileindices=create_random_profiles(realprofiles,annotationpool)
	out=open("../data/RandomProfilesOld.txt",'w')
	for profileid in randomprofiles:
		for annotation in randomprofiles[profileid]:
			out.write(profileid+"\t"+annotation+"\n")
	out.close()
	coeff,p=checkrandomization(realindices,randomindices)	
	coefflist.append(coeff)
	plist.append(p)

	medianoverlap=checkoverlap(realprofileindices,randomprofileindices)	


	print "Median Pearson's coeff",np.median(coefflist)
	print "Median p-value",np.median(p)
	print "Median overlap",medianoverlap




if __name__ == "__main__":
	import json
	import os
	import rdflib
	import matplotlib.pyplot as plt
	import sys
	import numpy as np
	from rdflib.plugin import register, Serializer, Parser
	import random
	import scipy.stats
	main()