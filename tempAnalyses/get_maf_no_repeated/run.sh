ml anaconda3/5.0

# join pileupts
sample=/scratch/users/paedugar/somaticMutationsProject/pileups/Skin_Sun_Exposed_Lower_leg/SRR662438.txt
blacklist_bed=/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/GRCh37.p13.genome.encode.blacklist.bed
rnaEditBed=/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/Human_AG_all_hg19_v2.txt.bed
rnaEditBed2=/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/hg19_darned.bed
exonBoundariesBed=/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/gencode.v19.genes.v7.patched_contigs.onlyExonBoundaries.bed
sampleAlias=SRR662438
PoN=/scratch/users/paedugar/somaticMutationsProject/mutationCount/panel_of_normal_mutations/n6_0.0_0.7.txt
repeated=/scratch/users/paedugar/somaticMutationsProject/mutationCount/repeated_mutations/n6_0.0_0.7/repeated_muts.txt
genome=/scratch/PI/hbfraser/GRCh37.p13.genome.encode.fa

ml anaconda3/5.0

# find mutations
python3 ../../mutationCalling/python/printMutationMap_fromPileup.py $sample True 6 40 0 1 output.unfiltered_map

# Do filters
        python3 ../../mutationCalling/python/filterAppend_distanceToBed.py output.unfiltered_map $blacklist_bed 1 blacklisted_region > $sampleAlias.fitlerBlacklist
        
        # Append filter: RNA edits
        python3 ../../mutationCalling/python/filterAppend_distanceToBed.py $sampleAlias.fitlerBlacklist $rnaEditBed 1 rna_edit > $sampleAlias.fitlerRnaEdit    
        
        # Append filter: RNA edits 2
        python3 ../../mutationCalling/python/filterAppend_distanceToBed.py $sampleAlias.fitlerRnaEdit $rnaEditBed2 1 rna_edit_darned > $sampleAlias.fitlerRnaEdit2
        
        # Append filter: splicing junction errors
        python3 ../../mutationCalling/python/filterAppend_distanceToBed.py $sampleAlias.fitlerRnaEdit2 $exonBoundariesBed 7 splicing_junction_error > $sampleAlias.fitlerExonBoundary
        
        # Append filer: sequencing error
        python3 ../../mutationCalling/python/filterAppend_sequenceError.py $sampleAlias.fitlerExonBoundary 0.0001 30 sequencing_error > $sampleAlias.fitlerSeqError
        
        # Append filter: clustered mutations
        python3 ../../mutationCalling/python/filterAppend_clusterMutations.py $sampleAlias.fitlerSeqError 100 3 clustered_mutation >  $sampleAlias.clustered 
        
        # Append filter: bcf stats
        # only if the exist in the original pileup file
         
        python3 ../../mutationCalling/python/filterAppend_statFromBCF.py $sampleAlias.clustered $sample bcf 0.05 0.05 0.05 0.05 0.05 > $sampleAlias.filtered_map
        
        # Includes cluster mutations but NOT variant distance bias
        awk -v FS="\t" -v OFS="\t" '{if ($8 == "PASS" || $8 == "clustered_mutation") print $0}' $sampleAlias.filtered_map | cut -f-7 > $sampleAlias.clean_map
        
        # Eliminates PoN
        ../../bin/mutationToBed  $sampleAlias.clean_map | bedtools intersect -v -a stdin -b $PoN |  cut -f 3 --complement > $sampleAlias.out_map
        
        # Eliminates repeated
        bedtools intersect -f 1 -loj -wo -a <(../../bin/mutationToBed $sampleAlias.out_map) -b <(../../bin/mutationToBed <(cut -f 2- $repeated)) | awk -v FS="\\t" -v OFS="\\t" '($4 == $12 && $5 != $13) || $13 == "."' | cut -f 1-2,4-8 > $sampleAlias.context
        
        python3 ../../mutationCalling/python/get_ref_base.py $genome $sampleAlias.context 0 1 4 5 > $sampleAlias.out_map

