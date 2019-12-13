#snakemake --snakefile Snakefile_validationPipeline.smk --cluster-config ../mappingMutationCalling_pipeline/cluster.json --cluster-status jobState --jobs 1900 --cluster "sbatch --parsable --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error}"

##--------------------------------------------------##
## GLOBAL

localrules: download_mutationSignatures, clean_mutationSignatures, create_simulatedMutationMat, plot_Mutationsignatures

NINDS_PER_SIGNATURE = "100"
SIGNATURES = "'Signature.7, Signature.10, Signature.13, Signature.17, Signature.18, Signature.20'"

##--------------------------------------------------##
## PIPELINE

rule all:
    input: 
        'validationPipeline/data/reconstructed_MutationSignatures.pdf',
        'validationPipeline/data/unbiased_reconstructed_W.txt'

rule download_mutationSignatures:
    '''
    Downloads the mutation signatures from Alexandrov et al 2013 Nature
    '''
    
    output:
        'validationPipeline/data/signatures.txt'
    shell:
        'wget -P validationPipeline/data/ ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/signatures.txt'

rule clean_mutationSignatures:
    '''
    Cleans mutation signture table for following step
    '''
    
    input:
        'validationPipeline/data/signatures.txt'
    output:
        'validationPipeline/data/signaturesClean.txt'
    
    shell:
        '''
        module load R/3.4.0
        Rscript validationPipeline/R/cleanMuationSignatures.r {input} {output}
        '''

rule create_simulatedMutationMat:
    '''
    Creates a simulated mutation matrix based on mutation signatures
    '''
    
    input:
        'validationPipeline/data/signaturesClean.txt'
    
    params:
        NINDS_PER_SIGNATURE,
        SIGNATURES
    
    output:
        'validationPipeline/data/simulatedMutationMatrix.txt'
    
    shell:
        '''
        module load R/3.4.0
        Rscript validationPipeline/R/simulateMutationMatrix.r {params} {input} {output}
        '''

rule call_MutationSignatures:
    '''
    Tries to recreate the same number of mutation signatures that were used 
    to construct the simulated mutation matrix
    '''
    
    input:
        'validationPipeline/data/simulatedMutationMatrix.txt'

    params:
        nRuns = 5000,
        nSign = len(SIGNATURES.split(",")),
        discoverNsign = "FALSE",
        contextLength = "3",
        useFreq = "FALSE",
        normilizeContext = "FALSE",
        nCores = 16,
        clusterWorkingDir = "/scratch/users/paedugar/deletemeSuperApply_mutationSignaturesValidationPipeline",
        oligoCountsFile = "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n3"
    
    output:
        'validationPipeline/data/reconstructed_W.txt',
        'validationPipeline/data/reconstructed_H.txt',
        'NULL' # NO plot
        
    threads: 16
    
    shell: 
        '''
        module load R/3.4.0
        Rscript R/mutationSignatures/mutationSignaturesNMF.r {params} {input} {output}
        '''

rule plot_Mutationsignatures:
    '''
    Plots the signatures obtaind from callMutationSignatures
    '''
    input:
        'validationPipeline/data/reconstructed_W.txt'
    output:
        'validationPipeline/data/reconstructed_MutationSignatures.pdf'
    shell:
        '''
        module load R/3.4.0
        Rscript R/mutationSignaturesProbPlot.r {input} {output}
        '''

rule find_unbiasedNsignatures:
    '''
    Finds unbiasly the best number of mutation signatures
    this should be equal to the number of muation signatures used to 
    create the simulated matrix
    '''
    
    input:
        'validationPipeline/data/simulatedMutationMatrix.txt'
    
    params:
        nRuns = 5000,
        nSign = len(SIGNATURES.split(",")),
        discoverNsign = "TRUE",
        contextLength = "3",
        useFreq = "FALSE",
        normilizeContext = "FALSE",
        nCores = 16,
        clusterWorkingDir = "/scratch/users/paedugar/deletemeSuperApply_mutationSignaturesValidationPipeline",
        oligoCountsFile = "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n3"
    
    output:
        'validationPipeline/data/unbiased_reconstructed_W.txt',
        'validationPipeline/data/unbiased_reconstructed_H.txt',
        'validationPipeline/data/unbiased_bestNsign.pdf'
    threads: 16

    shell: 
        '''
        module load R/3.4.0
        Rscript R/mutationSignatures/mutationSignaturesNMF.r {params} {input} {output}
        '''
    
