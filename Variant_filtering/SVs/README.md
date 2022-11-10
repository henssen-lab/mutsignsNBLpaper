# Commands to filter, merge, and format SVs in the discovery cohort

Collapse SVs (>50bp) from different callers (Svaba, Delly and Novobreak)<br/>
Filter them at least 2 callers except insertions<br/>
Including SVs from indels from Svaba and insertions PASS from delly2


## INTRACHROMOSOMAL SVs

### Get the intrachromosomal SVs per patient
python script_merge_intraSVs_svaba_svabafromindels_delly2_novobreak_5.0_collapsedups_nofiltering_newstrandinfo_precisioninfo_NEW.py svaba-sv.vcf.gz svaba-sv-fromindels.vcf.gz delly2-sv.vcf novobreak-sv.txt sampleID

### Format and filter the intrachromosomal SVs per patient
python script_filter_atleast2callers_4.0_withinsertions.py merge_intraSVs_diffvariantcallers_svaba_svaba_fromindels_delly2_novobreak_patient_sampleID_noVCFdups_collapsedups_nofilters_newinfostrands_precisioninfo.txt sampleID


## INTERCHROMOSOMAL SVs

### Get the interchromosomal SVs per patient
python script_merge_translocations_svaba_delly2_novobreak_4.0_collapsedups_nofiltering_newstrandinfo_precisioninfo_NEW.py svaba-sv.PASS.vcf.gz delly2-sv.vcf novobreak-sv.txt sampleID

### Format and filter the interchromosomal SVs per patient
python script_filter_atleast2callers_4.0_translocations.py merge_translocations_diffvariantcallers_svaba_delly2_novobreak_patient_sampleID_noVCFdups_collapsedups_nofilters_newinfostrands_precisioninfo.txt sampleID
