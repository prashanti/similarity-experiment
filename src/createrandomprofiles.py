
def get_profiles():
	profiles=dict()
	annotationpool=[]
	register('text/rdf+n3', Parser, 'rdflib.plugins.parsers.notation3', 'N3Parser')
	g = rdflib.Graph()
	result = g.parse('../data/realprofiles.ttl', format='n3')
	for stmt in g:
		profileid=stmt[0]
		annotation=stmt[2]
		if "VTO" in profileid:
			profileid=profileid.replace("http://purl.obolibrary.org/obo/","").replace("#profile","")
			if profileid not in profiles:
				profiles[profileid]=set()
			profiles[profileid].add(annotation)
			annotationpool.append(annotation)
	return profiles

def create_random_profiles(realprofiles):
	randomprofiles=dict()
	annotationpool=[]
	for profileid in realprofiles:
		profile=realprofiles[profileid]
		for ann in profile:
			annotationpool.append(ann)
	selected=set()		
	poolsize=len(annotationpool)
	for profileid in realprofiles:
		realprofilesize=len(realprofiles[profileid])
		if profileid not in randomprofiles:
			randomprofiles[profileid]=set()
		while (len(randomprofiles[profileid])<realprofilesize):
			randomindex=random.randint(0,poolsize-1)
			if (annotationpool[randomindex] not in randomprofiles[profileid]) and (randomindex not in selected):
				randomprofiles[profileid].add(annotationpool[randomindex])
				selected.add(randomindex)


def main():
	realprofiles=get_profiles()
	randomprofiles=create_random_profiles(realprofiles)





if __name__ == "__main__":
	import os
	import rdflib
	import sys
	from rdflib.plugin import register, Serializer, Parser
	import random
	main()