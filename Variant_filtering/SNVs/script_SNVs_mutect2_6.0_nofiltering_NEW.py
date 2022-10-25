#python

#Libraries
import sys
import subprocess
import re
import operator
import statistics
from itertools import chain
import gzip

#Declare lists
results_pointmut_mutect_snv = []
results_pointmut_TOTAL = []
results_POINTMUT_FINAL = []
pointmut_R = []
POINTMUTSV_uniq = []
POINTMUTSV_same = []
rec_block = []
rec_block_final = []
VCs_used = []
results_POINTMUT_FINAL_2 = []
results_POINTMUT_FINAL_3 = []

#Declare arguments
input_file_mutect_snv = sys.argv[1]
patient_id = sys.argv[2]


########## RETRIEVE SNVS FROM  THE CALLERS ##########

#Open and read the different variant callers results files

#mutect - using the PASS file of snv, all snv are PASS (we don't have to re-filter)
#Different format than the usual SV vcf from mutect
if input_file_mutect_snv != ".":
 with gzip.open(input_file_mutect_snv, 'rb') as f:
  mutect_file_snv = f.read().decode('utf-8')
  mutect_results_snv = re.split('\n',mutect_file_snv)
  for res5 in range(0,len(mutect_results_snv)-1):
   if len(re.findall(r'#',mutect_results_snv[res5]))>0:
    pass
   else:
    res5_line = re.split('\t',mutect_results_snv[res5])
    #pass vcf to smufin style
    res5_line_chrom2 = res5_line[0]
    if len(res5_line[3])==1 and len(res5_line[4])==1 and len(res5_line[3])==len(res5_line[4]):
     res5_line_pos2 = res5_line[1]
     sv_type = "snv"
     length = 1
     strands = ".."
     res5_line_pointmut = list(res5_line[0:2])
     res5_line_pointmut.append(res5_line_chrom2)
     res5_line_pointmut.append(res5_line_pos2)
     res5_line_pointmut.append('MUTECT_snv=' + res5_line[6] + ',.,' + res5_line[3] + ',' + res5_line[4])
     res5_line_pointmut.append(strands)
     res5_line_pointmut.append(sv_type)
     res5_line_pointmut.append(length)
     res5_line_pointmut.insert(0, patient_id)
     results_pointmut_mutect_snv.append(res5_line_pointmut)
    else:
     pass
else:
 pass


if len(results_pointmut_mutect_snv)>0:
 results_pointmut_TOTAL.append(list(results_pointmut_mutect_snv))

#Unnest the lists
results_POINTMUT_FINAL = list(chain.from_iterable(results_pointmut_TOTAL))
results_POINTMUT_FINAL = sorted(results_POINTMUT_FINAL, key=operator.itemgetter(1,2,4),  reverse=True)

########## COLLAPSE DUPLICATED SNVS FROM DIFFERENT/SAME CALLERS ##########

#Change chrX and chrY for 23 and 24 to order them as integers
for T in range(0,len(results_POINTMUT_FINAL)-1):
 if results_POINTMUT_FINAL[T][1]=='X' and results_POINTMUT_FINAL[T][3]=='Y':
  results_POINTMUT_FINAL[T][1] = '23'
  results_POINTMUT_FINAL[T][3] = '24'
 elif results_POINTMUT_FINAL[T][1]=='Y' and results_POINTMUT_FINAL[T][3]=='X':
  results_POINTMUT_FINAL[T][1] = '24'
  results_POINTMUT_FINAL[T][3] = '23'
 elif results_POINTMUT_FINAL[T][1]=='X':
  results_POINTMUT_FINAL[T][1] = '23'
 elif results_POINTMUT_FINAL[T][1]=='Y':
  results_POINTMUT_FINAL[T][1] = '24'
 elif results_POINTMUT_FINAL[T][3]=='X':
  results_POINTMUT_FINAL[T][3] = '23'
 elif results_POINTMUT_FINAL[T][3]=='Y':
  results_POINTMUT_FINAL[T][3] = '24'
 else:
  pass


#Order the SNVS this way: first lowest position ex. 2 67587578 2 13231331 (original: 2 13231331 2 67587578)
for R in range(0,len(results_POINTMUT_FINAL)):
 #I remove the chromosomes that are not 1:22,X,Y and remove results with length<50bp
 if len(re.findall(r'hs37d5',results_POINTMUT_FINAL[R][1]))>0 or len(re.findall(r'_random',results_POINTMUT_FINAL[R][1]))>0 or len(re.findall(r'Un_gl',results_POINTMUT_FINAL[R][1]))>0 or results_POINTMUT_FINAL[R][1] == 'M' or len(re.findall(r'hs37d5',results_POINTMUT_FINAL[R][3]))>0 or len(re.findall(r'_random',results_POINTMUT_FINAL[R][3]))>0 or len(re.findall(r'Un_gl',results_POINTMUT_FINAL[R][3]))>0 or results_POINTMUT_FINAL[R][3] == 'M' or len(re.findall(r'GL',results_POINTMUT_FINAL[R][1]))>0 or len(re.findall(r'GL',results_POINTMUT_FINAL[R][3]))>0 or len(re.findall(r'MT',results_POINTMUT_FINAL[R][1]))>0 or len(re.findall(r'MT',results_POINTMUT_FINAL[R][3]))>0 or int(results_POINTMUT_FINAL[R][8])>1:
  pass
 else:
  results_POINTMUT_FINAL_2.append(results_POINTMUT_FINAL[R])


#results_POINTMUT_FINAL without rare chromosomes
results_POINTMUT_FINAL = []
results_POINTMUT_FINAL = results_POINTMUT_FINAL_2

#Re-sort list
results_POINTMUT_FINAL = sorted(results_POINTMUT_FINAL, key=operator.itemgetter(1,2,4),  reverse=False)

#insert strand info and svtype on info field
for s in range(0,len(results_POINTMUT_FINAL)):
  collapse_info_strand = results_POINTMUT_FINAL[s][5] + ',' + results_POINTMUT_FINAL[s][6] + ',' + results_POINTMUT_FINAL[s][7] + ',' + str(results_POINTMUT_FINAL[s][8])
  results_POINTMUT_FINAL[s][5] = collapse_info_strand
  results_POINTMUT_FINAL_3.append(results_POINTMUT_FINAL[s][0:6])

#results_POINTMUT_FINAL without rare chromosomes and with strands in info
results_POINTMUT_FINAL = []
results_POINTMUT_FINAL = results_POINTMUT_FINAL_3

POINTMUTSV_uniq = results_POINTMUT_FINAL

########## SAVE RESULTS IN FILE ##########

#Print callers names for each patient in the file name

if input_file_mutect_snv != ".":
 f_name_mutect_snv = "mutect_snv"
else:
 f_name_mutect_snv = "X"

#Open the file to write in
results_patient = open('/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/SNV_' + f_name_mutect_snv + '_' + 'patient_' + patient_id + '_nofilters.txt','w')

#Re-sort list
POINTMUTSV_uniq = sorted(POINTMUTSV_uniq, key=operator.itemgetter(1,3,2,4),  reverse=False)

#Change chr 23 and 24 for X and Y again
for T in range(0,len(POINTMUTSV_uniq)):
 if POINTMUTSV_uniq[T][1]=='23' and POINTMUTSV_uniq[T][3]=='24':
  POINTMUTSV_uniq[T][1] = 'X'
  POINTMUTSV_uniq[T][3] = 'Y'
 elif POINTMUTSV_uniq[T][1]=='24' and results_POINTMUT_FINAL[T][3]=='23':
  POINTMUTSV_uniq[T][1] = 'Y'
  POINTMUTSV_uniq[T][3] = 'X'
 elif POINTMUTSV_uniq[T][1]=='23':
  POINTMUTSV_uniq[T][1] = 'X'
 elif POINTMUTSV_uniq[T][1]=='24':
  POINTMUTSV_uniq[T][1] = 'Y'
 elif POINTMUTSV_uniq[T][3]=='23':
  POINTMUTSV_uniq[T][3] = 'X'
 elif POINTMUTSV_uniq[T][3]=='24':
  POINTMUTSV_uniq[T][3] = 'Y'
 else:
  pass


#Print results
for r in range(0,len(POINTMUTSV_uniq)):
 results_patient.write('\t'.join(str(v) for v in POINTMUTSV_uniq[r][0:6]) + '\n')

results_patient.close()

print(patient_id + ' completed')
