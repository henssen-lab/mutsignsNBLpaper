import sys
import gzip
import re

#declare arguments
filename= sys.argv[1] #delly output to transform to palipmpsest
id_sample= sys.argv[2] #patient id

#Open file to read mutations
delly_o = open(filename,'r')

#Open the file to write in
fin_file = open(filename +'.palimpsest.vcf', 'w')

#Select mutations for Yapsa
for line in delly_o:
	if not line.startswith("#"): # exclude header lines
		chro = line.split()[0]
		pos = line.split()[1]
		qual =line.split()[6]
		info = line.split()[7]
		info_typeof = info.split(";")[1]
		word_type = info_typeof.split("=")[1]
		info_chr2 = info.split(";")[3]
		chr2 = info_chr2.split("=")[1]
		info_pos2 = info.split(";")[4]
		pos2 = info_pos2.split("=")[1]
		if str(qual) == str("PASS"):
			fin_file.write(str(id_sample) + "\t" + str(word_type) + "\t" + str(chro) + "\t" + str(pos) + "\t" + str(chr2) + "\t" + str(pos2) + "\t" + str(".") + "\t" + str(".") + "\t" + str(".") + "\n")
		else:
			pass
fin_file.close()
