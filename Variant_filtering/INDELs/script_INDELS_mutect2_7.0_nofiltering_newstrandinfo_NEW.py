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
results_indel_mutect_indels = []
results_indel_TOTAL = []
results_INDEL_FINAL = []
indel_R = []
INDELmut_uniq = []
INDELmut_same = []
rec_block = []
rec_block_final = []
VCs_used = []
results_INDEL_FINAL_2 = []
results_INDEL_FINAL_3 = []

#Declare arguments
input_file_mutect_indels = sys.argv[1]
patient_id = sys.argv[2]


########## RETRIEVE INDELS FROM ALL THE CALLERS ##########

#Open and read the different variant callers results files

#mutect - using the PASS file of indels, all indels are PASS (we don't have to re-filter)
#Different format than the usual mut vcf from mutect
if input_file_mutect_indels != ".":
 with gzip.open(input_file_mutect_indels, 'rb') as f:
  mutect_file_indels = f.read().decode('utf-8')
  mutect_results_indels = re.split('\n',mutect_file_indels)
  for res5 in range(0,len(mutect_results_indels)-1):
   if len(re.findall(r'#',mutect_results_indels[res5]))>0:
    pass
   else:
    res5_line = re.split('\t',mutect_results_indels[res5])
    res5_line_chrom2 = res5_line[0]
    if len(res5_line[3])>len(res5_line[4]):
     res5_line_pos2 = str(int(res5_line[1]) + int(len(res5_line[3]))-1)
     sv_type = "deletion"
     length = abs(int(res5_line[1])-int(res5_line_pos2))
     strands = "+-"
     res5_line_indel = list(res5_line[0:2])
     res5_line_indel.append(res5_line_chrom2)
     res5_line_indel.append(res5_line_pos2)
     res5_line_indel.append('MUTECT_indels=' + res5_line[6] + ',.,' + res5_line[3] + ',' + res5_line[4])
     res5_line_indel.append(strands)
     res5_line_indel.append(sv_type)
     res5_line_indel.append(length)
     res5_line_indel.insert(0, patient_id)
     results_indel_mutect_indels.append(res5_line_indel)
    elif len(res5_line[3])<len(res5_line[4]):
     res5_line_pos2 = str(int(res5_line[1]) + 1)
     sv_type = "insertion"
     length = int(len(res5_line[4])) - 1 
     strands = "++"
     res5_line_indel = list(res5_line[0:2])
     res5_line_indel.append(res5_line_chrom2)
     res5_line_indel.append(res5_line_pos2)
     res5_line_indel.append('MUTECT_indels=' + res5_line[6] + ',.,' + res5_line[3] + ',' + res5_line[4])
     res5_line_indel.append(strands)
     res5_line_indel.append(sv_type)
     res5_line_indel.append(length)
     res5_line_indel.insert(0, patient_id)
     results_indel_mutect_indels.append(res5_line_indel)
    else:
     pass
else:
 pass


if len(results_indel_mutect_indels)>0:
 results_indel_TOTAL.append(list(results_indel_mutect_indels))

#Unnest the lists
results_INDEL_FINAL = list(chain.from_iterable(results_indel_TOTAL))
results_INDEL_FINAL = sorted(results_INDEL_FINAL, key=operator.itemgetter(1,2,4),  reverse=True)

########## COLLAPSE DUPLICATED INDELS FROM DIFFERENT/SAME CALLERS ##########

#Change chrX and chrY for 23 and 24 to order them as integers
for T in range(0,len(results_INDEL_FINAL)-1):
 if results_INDEL_FINAL[T][1]=='X' and results_INDEL_FINAL[T][3]=='Y':
  results_INDEL_FINAL[T][1] = '23'
  results_INDEL_FINAL[T][3] = '24'
 elif results_INDEL_FINAL[T][1]=='Y' and results_INDEL_FINAL[T][3]=='X':
  results_INDEL_FINAL[T][1] = '24'
  results_INDEL_FINAL[T][3] = '23'
 elif results_INDEL_FINAL[T][1]=='X':
  results_INDEL_FINAL[T][1] = '23'
 elif results_INDEL_FINAL[T][1]=='Y':
  results_INDEL_FINAL[T][1] = '24'
 elif results_INDEL_FINAL[T][3]=='X':
  results_INDEL_FINAL[T][3] = '23'
 elif results_INDEL_FINAL[T][3]=='Y':
  results_INDEL_FINAL[T][3] = '24'
 else:
  pass

#Order the indelmuts this way: first lowest position ex. 2 67587578 2 13231331 (original: 2 13231331 2 67587578)
for R in range(0,len(results_INDEL_FINAL)):
 #I remove the chromosomes that are not 1:22,X,Y and remove results with length<50bp
 if len(re.findall(r'hs37d5',results_INDEL_FINAL[R][1]))>0 or len(re.findall(r'_random',results_INDEL_FINAL[R][1]))>0 or len(re.findall(r'Un_gl',results_INDEL_FINAL[R][1]))>0 or results_INDEL_FINAL[R][1] == 'M' or len(re.findall(r'hs37d5',results_INDEL_FINAL[R][3]))>0 or len(re.findall(r'_random',results_INDEL_FINAL[R][3]))>0 or len(re.findall(r'Un_gl',results_INDEL_FINAL[R][3]))>0 or results_INDEL_FINAL[R][3] == 'M' or len(re.findall(r'GL',results_INDEL_FINAL[R][1]))>0 or len(re.findall(r'GL',results_INDEL_FINAL[R][3]))>0 or len(re.findall(r'MT',results_INDEL_FINAL[R][1]))>0 or len(re.findall(r'MT',results_INDEL_FINAL[R][3]))>0 or int(results_INDEL_FINAL[R][8])>50:
  pass
 else:
  results_INDEL_FINAL_2.append(results_INDEL_FINAL[R])


#results_INDEL_FINAL without rare chromosomes
results_INDEL_FINAL = []
results_INDEL_FINAL = results_INDEL_FINAL_2

for R in range(0,len(results_INDEL_FINAL)-1):
 #Change the order of chr-pos
 if int(results_INDEL_FINAL[R][2])>int(results_INDEL_FINAL[R][4]):
  new_pos_1 = results_INDEL_FINAL[R][4]
  new_pos_2 = results_INDEL_FINAL[R][2]
  if results_INDEL_FINAL[R][6]=='++':
   new_strand = '--'
  elif results_INDEL_FINAL[R][6]=='--':
   new_strand = '++'
  else:
   new_strand = results_INDEL_FINAL[R][6][::-1]
  results_INDEL_FINAL[R][2] = new_pos_1
  results_INDEL_FINAL[R][4] = new_pos_2
  results_INDEL_FINAL[R][6] = new_strand
 else:
  pass

#Re-sort list
results_INDEL_FINAL = sorted(results_INDEL_FINAL, key=operator.itemgetter(1,2,4),  reverse=False)

#insert strand info and svtype on info field
for s in range(0,len(results_INDEL_FINAL)):
  collapse_info_strand = results_INDEL_FINAL[s][5] + ',' + results_INDEL_FINAL[s][6] + ',' + results_INDEL_FINAL[s][7] + ',' + str(results_INDEL_FINAL[s][8])
  results_INDEL_FINAL[s][5] = collapse_info_strand
  results_INDEL_FINAL_3.append(results_INDEL_FINAL[s][0:6])

#results_INDEL_FINAL without rare chromosomes and with strands in info
results_INDEL_FINAL = []
results_INDEL_FINAL = results_INDEL_FINAL_3

INDELmut_uniq = results_INDEL_FINAL


########## SAVE RESULTS IN FILE ##########

#Print callers names for each patient in the file name

if input_file_mutect_indels != ".":
 f_name_mutect_indels = "mutect_indels"
else:
 f_name_mutect_indels = "X"

#Open the file to write in
results_patient = open('/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/INDELS_' + f_name_mutect_indels + '_patient_' + patient_id + '_nofilters_newinfostrands.txt','w')

#Re-sort list
INDELmut_uniq = sorted(INDELmut_uniq, key=operator.itemgetter(1,3,2,4),  reverse=False)

#Change chr 23 and 24 for X and Y again
for T in range(0,len(INDELmut_uniq)):
 if INDELmut_uniq[T][1]=='23' and INDELmut_uniq[T][3]=='24':
  INDELmut_uniq[T][1] = 'X'
  INDELmut_uniq[T][3] = 'Y'
 elif INDELmut_uniq[T][1]=='24' and results_INDEL_FINAL[T][3]=='23':
  INDELmut_uniq[T][1] = 'Y'
  INDELmut_uniq[T][3] = 'X'
 elif INDELmut_uniq[T][1]=='23':
  INDELmut_uniq[T][1] = 'X'
 elif INDELmut_uniq[T][1]=='24':
  INDELmut_uniq[T][1] = 'Y'
 elif INDELmut_uniq[T][3]=='23':
  INDELmut_uniq[T][3] = 'X'
 elif INDELmut_uniq[T][3]=='24':
  INDELmut_uniq[T][3] = 'Y'
 else:
  pass


#Print results
for r in range(0,len(INDELmut_uniq)):
 results_patient.write('\t'.join(str(v) for v in INDELmut_uniq[r][0:6]) + '\n')

results_patient.close()

print(patient_id + ' completed')
