#Libraries
import sys
import subprocess
import re
import operator
from itertools import chain
import gzip
from collections import Counter

#Declare lists
FINAL_vcfs_results = []

#Declare arguments
input_file_vcfs = sys.argv[1]
patient_id = sys.argv[2]

##########

#Open VCFs collapsed results file
file_vcfs = open(input_file_vcfs,'r')
vcfs_results = re.split('\n',file_vcfs.read())
for res in range(0,len(vcfs_results)-1):
 res_line = re.split('\t',vcfs_results[res])
 info = res_line[5]
 if len(re.findall(r'PASS',info))>0:
  info_list = [info] #list with all the info field for each caller
  type_list = [info] #list with all the info field for each caller
  length_list = [info] #list with all the info field for each caller
  strand_list = [info] #list with all the info field for each caller
  ref_list = [info] #list with all the info field for each caller
  alt_list = [info] #list with all the info field for each caller
  for i in range(0,len(info_list)):
   caller = re.split('=',info_list[i])[0]
   sv_type = re.split('=',info_list[i])[1].split(',')[5]
   length = re.split('=',info_list[i])[1].split(',')[6]
   strand = re.split('=',info_list[i])[1].split(',')[4]
   ref = re.split('=',info_list[i])[1].split(',')[2]
   alt = re.split('=',info_list[i])[1].split(',')[3]
   info_list[i] = caller
   type_list[i] = sv_type
   length_list[i] = length
   strand_list[i] = strand
   ref_list[i] = ref
   alt_list[i] = alt
   if sv_type == 'insertion':
    length_list[i] = length
   else:
    length_list[i] = abs(int(res_line[2])-int(res_line[4]))
  first_svtype = list(Counter(type_list))[0]
  first_strand = list(Counter(strand_list))[0]
  first_ref = list(Counter(ref_list))[0]
  first_alt = list(Counter(alt_list))[0]
  res_line.append(first_svtype)
  res_line.append(first_strand)
  res_line.append(length_list[0])
  res_line.append(first_ref)
  res_line.append(first_alt)
  FINAL_vcfs_results.append(res_line)
 else:
  pass

#Open the file to write in
results_patient = open(input_file_vcfs+'.FILTERPASS.svtype','w')

#Print results
for r in range(0,len(FINAL_vcfs_results)):
 results_patient.write('\t'.join(str(v) for v in FINAL_vcfs_results[r]) + '\n')

results_patient.close()

print(patient_id + ' completed')
