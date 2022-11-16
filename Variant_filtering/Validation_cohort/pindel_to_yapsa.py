import sys
import gzip
import re

#declare arguments
filename= sys.argv[1] #pindel output to transform to vcf
id_sample= sys.argv[2] #patient id

#Open file to read mutations
pindel_o = open(filename,'r')


#Open the file to write in
fin_file = open(filename +'.yapsa.vcf', 'w')

#Select mutations for Yapsa
for line in pindel_o:
	if not line.startswith("#"): # exclude header lines
		chro = line.split()[0]
		pos = line.split()[1]
		ref = line.split()[3]
		alt = line.split()[4]
		info = line.split()[7]
		info_typeof = info.split(";")[3]
		word_type = info_typeof.split("=")[0]
		if str(word_type) == str("SVLEN"):
			info_typeof = info.split(";")[4]
			type_of = info_typeof.split("=")[1]
			info_svlen = info.split(";")[3]
			if type_of == str("RPL"):
				info_ntlen = info.split(";")[5]
				ntlen = info_ntlen.split("=")[1]
			else:
				pass
		else:
			info_svlen = info.split(";")[2]
			info_typeof = info.split(";")[3]
			type_of = info_typeof.split("=")[1]
			if type_of == str("RPL"):
				info_ntlen = info.split(";")[4]
				ntlen = info_ntlen.split("=")[1]
			else:
				pass
		#type_of = info_typeof.split("=")[1]
		#info_svlen = info.split(";")[3]
		#svlen = info_svlen.split("-")[1]
		#svlen_sis = info_svlen.split("=")[1]
		#First filter by type Deletions
		if type_of == str("DEL"):
			svlen = info_svlen.split("-")[1]
			if int(svlen) <= int(50):
				fin_file.write(str(id_sample) + "\t" + str(chro) + "\t" + str(pos) + "\t" + str(ref) + "\t" + str(alt) + "\n")
			else:
				pass
		elif type_of == str("INS"):
			svlen_sis = info_svlen.split("=")[1]
			if int(svlen_sis) <= int(50) :
				fin_file.write(str(id_sample) + "\t" + str(chro) + "\t" + str(pos) + "\t" + str(ref) + "\t" + str(alt) + "\n")
			else:
				pass
		elif type_of == str("RPL"):
			svlen = info_svlen.split("-")[1]
			if int(svlen) <= int(50) and int(ntlen) <= int(50):
				fin_file.write(str(id_sample) + "\t" + str(chro) + "\t" + str(pos) + "\t" + str(ref) + "\t" + str(alt) + "\n")
			else:
				pass
		else:
			pass
fin_file.close()


