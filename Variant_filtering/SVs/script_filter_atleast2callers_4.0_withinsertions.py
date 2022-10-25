#python

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
 info_list = list(re.split(';',info)) #list with all the info field for each caller
 type_list = list(re.split(';',info)) #list with all the info field for each caller
 length_list = list(re.split(';',info)) #list with all the info field for each caller
 strand_list = list(re.split(';',info)) #list with all the info field for each caller
 for i in range(0,len(info_list)):
  caller = re.split('=',info_list[i])[0]
  sv_type = re.split('=',info_list[i])[1].split(',')[4]
  length = re.split('=',info_list[i])[1].split(',')[5]
  strand = re.split('=',info_list[i])[1].split(',')[2]
  info_list[i] = caller
  type_list[i] = sv_type
  strand_list[i] = strand
  if sv_type == 'insertion':
   length_list[i] = length
  else:
    length_list[i] = abs(int(res_line[2])-int(res_line[4]))
 first_svtype = list(Counter(type_list))[0]
 first_strand = list(Counter(strand_list))[0]
 if len(set(type_list))>1 and first_svtype == 'unknown' and len(set(strand_list))>1:
  second_svtype = list(Counter(type_list))[1]
  second_strand = list(Counter(strand_list))[1]
  first_svtype = second_svtype
  first_strand = second_strand
  res_line.append(first_svtype)
  res_line.append(first_strand)
 elif len(set(type_list))>1 and first_svtype == 'unknown' and len(set(strand_list))==1:
  second_svtype = list(Counter(type_list))[1]
  first_svtype = second_svtype
  res_line.append(first_svtype)
  res_line.append(first_strand)
 else:
  res_line.append(first_svtype)
  res_line.append(first_strand)
 res_line.append(length_list[0])
 unique_list = set(info_list) #unique list with the callers names
 if len(unique_list)>1: #IF SV SUPPORTED FOR 2 DIFFERENT CALLERS I KEEP IT
  FINAL_vcfs_results.append(res_line)
 elif len(unique_list)==1 and first_svtype == 'insertion' and len(re.findall(r'PASS',info))>0:
  FINAL_vcfs_results.append(res_line)
 else:
  pass

#Open the file to write in
results_patient = open(input_file_vcfs+'.FILTER2CALLERS.svtype','w')

#Print results
for r in range(0,len(FINAL_vcfs_results)):
 results_patient.write('\t'.join(str(v) for v in FINAL_vcfs_results[r]) + '\n')

results_patient.close()

print(patient_id + ' completed')
