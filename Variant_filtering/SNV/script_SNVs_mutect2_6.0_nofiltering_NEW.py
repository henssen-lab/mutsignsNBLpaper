#python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/aux.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/CB2009/mutect/mutect-sv.PASS.vcf.gz /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/CB2009/delly2/delly2-svs.pre.vcf /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/novobreak_results_hg19_joern/MERGE_BerlinCohort_PeiferWGS_NovobreakResults_filter_hg19_goodformat.txt CB2009

#Libraries
import sys
import subprocess
import re
import operator
import statistics
from itertools import chain
import gzip

#Declare lists
results_intra_mutect_snv = []
results_intra_TOTAL = []
results_INTRA_FINAL = []
intra_R = []
INTRASV_uniq = []
INTRASV_same = []
rec_block = []
rec_block_final = []
VCs_used = []
results_INTRA_FINAL_2 = []
results_INTRA_FINAL_3 = []

#Declare arguments
input_file_mutect_snv = sys.argv[1]
patient_id = sys.argv[2]



########## RETRIEVE TRANSLOCATIONS FROM ALL THE CALLERS ##########

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
     res5_line_intra = list(res5_line[0:2])
     res5_line_intra.append(res5_line_chrom2)
     res5_line_intra.append(res5_line_pos2)
     res5_line_intra.append('MUTECT_snv=' + res5_line[6] + ',.,' + res5_line[3] + ',' + res5_line[4])
     res5_line_intra.append(strands)
     res5_line_intra.append(sv_type)
     res5_line_intra.append(length)
     res5_line_intra.insert(0, patient_id)
     results_intra_mutect_snv.append(res5_line_intra)
    else:
     pass
else:
 pass

#Join all lists and Remove the duplicates using deduplicate function, except smufin and novobreak - we don't have SV ids provided by those callers

if len(results_intra_mutect_snv)>0:
 results_intra_TOTAL.append(list(results_intra_mutect_snv))

#Unnest the lists
results_INTRA_FINAL = list(chain.from_iterable(results_intra_TOTAL))
results_INTRA_FINAL = sorted(results_INTRA_FINAL, key=operator.itemgetter(1,2,4),  reverse=True)

########## COLLAPSE DUPLICATED INTRASVS FROM DIFFERENT/SAME CALLERS ##########

#Change chrX and chrY for 23 and 24 to order them as integers
for T in range(0,len(results_INTRA_FINAL)-1):
 if results_INTRA_FINAL[T][1]=='X' and results_INTRA_FINAL[T][3]=='Y':
  results_INTRA_FINAL[T][1] = '23'
  results_INTRA_FINAL[T][3] = '24'
 elif results_INTRA_FINAL[T][1]=='Y' and results_INTRA_FINAL[T][3]=='X':
  results_INTRA_FINAL[T][1] = '24'
  results_INTRA_FINAL[T][3] = '23'
 elif results_INTRA_FINAL[T][1]=='X':
  results_INTRA_FINAL[T][1] = '23'
 elif results_INTRA_FINAL[T][1]=='Y':
  results_INTRA_FINAL[T][1] = '24'
 elif results_INTRA_FINAL[T][3]=='X':
  results_INTRA_FINAL[T][3] = '23'
 elif results_INTRA_FINAL[T][3]=='Y':
  results_INTRA_FINAL[T][3] = '24'
 else:
  pass


#Order the intraSVs this way: first lowest position ex. 2 67587578 2 13231331 (original: 2 13231331 2 67587578)
for R in range(0,len(results_INTRA_FINAL)):
 #I remove the chromosomes that are not 1:22,X,Y and remove results with length<50bp
 if len(re.findall(r'hs37d5',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'_random',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'Un_gl',results_INTRA_FINAL[R][1]))>0 or results_INTRA_FINAL[R][1] == 'M' or len(re.findall(r'hs37d5',results_INTRA_FINAL[R][3]))>0 or len(re.findall(r'_random',results_INTRA_FINAL[R][3]))>0 or len(re.findall(r'Un_gl',results_INTRA_FINAL[R][3]))>0 or results_INTRA_FINAL[R][3] == 'M' or len(re.findall(r'GL',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'GL',results_INTRA_FINAL[R][3]))>0 or len(re.findall(r'MT',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'MT',results_INTRA_FINAL[R][3]))>0 or int(results_INTRA_FINAL[R][8])>1:
  pass
 else:
  results_INTRA_FINAL_2.append(results_INTRA_FINAL[R])


#results_INTRA_FINAL without rare chromosomes
results_INTRA_FINAL = []
results_INTRA_FINAL = results_INTRA_FINAL_2

#Re-sort list
results_INTRA_FINAL = sorted(results_INTRA_FINAL, key=operator.itemgetter(1,2,4),  reverse=False)

#insert strand info and svtype on info field
for s in range(0,len(results_INTRA_FINAL)):
  collapse_info_strand = results_INTRA_FINAL[s][5] + ',' + results_INTRA_FINAL[s][6] + ',' + results_INTRA_FINAL[s][7] + ',' + str(results_INTRA_FINAL[s][8])
  results_INTRA_FINAL[s][5] = collapse_info_strand
  results_INTRA_FINAL_3.append(results_INTRA_FINAL[s][0:6])

#results_INTRA_FINAL without rare chromosomes and with strands in info
results_INTRA_FINAL = []
results_INTRA_FINAL = results_INTRA_FINAL_3

INTRASV_uniq = results_INTRA_FINAL

########## SAVE RESULTS IN FILE ##########

#Print callers names for each patient in the file name

if input_file_mutect_snv != ".":
 f_name_mutect_snv = "mutect_snv"
else:
 f_name_mutect_snv = "X"

#Open the file to write in
results_patient = open('/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/SNV_' + f_name_mutect_snv + '_' + 'patient_' + patient_id + '_nofilters.txt','w')

#Re-sort list
INTRASV_uniq = sorted(INTRASV_uniq, key=operator.itemgetter(1,3,2,4),  reverse=False)

#Change chr 23 and 24 for X and Y again
for T in range(0,len(INTRASV_uniq)):
 if INTRASV_uniq[T][1]=='23' and INTRASV_uniq[T][3]=='24':
  INTRASV_uniq[T][1] = 'X'
  INTRASV_uniq[T][3] = 'Y'
 elif INTRASV_uniq[T][1]=='24' and results_INTRA_FINAL[T][3]=='23':
  INTRASV_uniq[T][1] = 'Y'
  INTRASV_uniq[T][3] = 'X'
 elif INTRASV_uniq[T][1]=='23':
  INTRASV_uniq[T][1] = 'X'
 elif INTRASV_uniq[T][1]=='24':
  INTRASV_uniq[T][1] = 'Y'
 elif INTRASV_uniq[T][3]=='23':
  INTRASV_uniq[T][3] = 'X'
 elif INTRASV_uniq[T][3]=='24':
  INTRASV_uniq[T][3] = 'Y'
 else:
  pass


#Print results
for r in range(0,len(INTRASV_uniq)):
 results_patient.write('\t'.join(str(v) for v in INTRASV_uniq[r][0:6]) + '\n')

results_patient.close()

print(patient_id + ' completed')
