
def check_num_of_profiles():
	infile=open("../data/realprofiles.ttl")
	outfile=open("../data/VTOprofiles.tsv",'w')
	for line in infile:
		if "VTO" in line:
			outfile.write(line.replace(" a ","\t"))
	infile.close()
	outfile.close()

def get_profiles():
	profiles=dict()
	annotationpool=[]
	infile=open("../data/VTOprofiles.tsv")
	for line in infile:
		profileid,annotations=line.replace(" .","").strip().split("\t")
		profileid=profileid.replace("<http://purl.obolibrary.org/obo/","").replace("#profile>","")

		annotationlist=annotations.split(" , ")
		
		if profileid not in profiles:
			profiles[profileid]=set()
		for annotation in annotationlist:
			profiles[profileid].add(annotation.strip())
	infile.close()
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
	check_num_of_profiles()
	realprofiles=get_profiles()
	randomprofiles=create_random_profiles(realprofiles)





if __name__ == "__main__":
	import os
	import random
	main()