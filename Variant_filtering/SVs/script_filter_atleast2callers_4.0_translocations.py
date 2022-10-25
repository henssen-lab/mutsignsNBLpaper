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
 strand_list = list(re.split(';',info))
 for i in range(0,len(info_list)):
  caller = re.split('=',info_list[i])[0]
  sv_type = "translocation"
  strand = re.split('=',info_list[i])[1].split(',')[2]
  info_list[i] = caller
  type_list[i] = sv_type
  strand_list[i] = strand
 first_svtype = list(Counter(type_list))[0]
 first_strand = list(Counter(strand_list))[0]
 if first_strand == '..' and len(set(strand_list))>1:
  second_strand = list(Counter(strand_list))[1]
  first_strand = second_strand
  res_line.append(first_svtype)
  res_line.append(first_strand)
 elif first_strand == '..' and len(set(strand_list))==1:
  res_line.append(first_svtype)
  res_line.append(first_strand)
 else:
  res_line.append(first_svtype)
  res_line.append(first_strand)
 res_line.append("0")
 unique_list = set(info_list) #unique list with the callers names
 if len(unique_list)>1: #IF SV SUPPORTED FOR 2 DIFFERENT CALLERS I KEEP IT
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
