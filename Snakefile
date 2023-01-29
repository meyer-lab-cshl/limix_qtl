import os

from snakemake.utils import min_version

shell.prefix("set -euo pipefail;")

##### set minimum snakemake version #####
min_version("6.10.0")


##### functions #####
def extendChunk(chunk):
    chunkSplitted = chunk.split("_")
    return chunkSplitted[0]+":"+chunkSplitted[1]+"-"+chunkSplitted[2]

## Variables ##
## Files, these are populated now with the test examples and we use here plink genotypes.
genotypeFile = 'Limix_QTL/test_data/Genotypes/Geuvadis'  ##Genotype without file extension. Please update flag in the runner to reflect --plink or --bgen
kinshipFile = 'Limix_QTL/test_data/Genotypes/Geuvadis_chr1_kinship.normalized.txt'
annotationFile = 'Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot_chr1_small.txt'
phenotypeFile = 'Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_Expr.txt'
covariateFile = 'Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_covariates.txt'
randomEffFile = '' #no second random effect in this example
sampleMappingFile = 'Limix_QTL/test_data/Geuvadis_CEU_gte.txt'
outputFolder = 'OutGeuvadis'

## Settings
nGenes = 5 # recommend to set to 50; 5 used here for toy dataset
startPos = 0
endOffset = 1000000000
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
windowSize = '1000000'
hwequilibrium = '0.000001'
FDR = '0.05'
## End Variables ##

finalQTLRun = outputFolder + '/qtl_results_all.txt'
topQTL = outputFolder + '/top_qtl_results_all.txt'


rule all:
    input:
        finalQTLRun, topQTL

rule generate_chunks:
    input:
        annotation=annotationFile
    output:
        chunks=outputFolder + "/chunks.txt"
    log:
        outputFolder + "log/generate_chunks.log"
    params:
        nGenes = nGenes,
        startPos = startPos,
        endOffset = endOffset
    script:
        "Limix_QTL/scripts/generate_chunks.R"

checkpoint setup_chunk_analysis:
    input:
        chunks=outputFolder + "/chunks.txt"
    output:
        directory(outputFolder + "/chunks")
    shell:
        """
        mkdir {output}
        for chunk in $(cat {input.chunks}); do
            echo $chunk > {output}/${{chunk//[:-]/_}}_info
        done
        """

rule run_qtl_mapping:
    input:
        chunkinfo = outputFolder + "/chunks/{chunk}_info",
        af = annotationFile,
        pf = phenotypeFile,
        cf = covariateFile,
        kf = kinshipFile,
        smf = sampleMappingFile
    output:
        outputFolder + "/chunks/qtl_results_{chunk}.h5"
    params:
        chunkstr=lambda wildcards: extendChunk(wildcards.chunk),
        gen=genotypeFile,
        od = outputFolder + "/chunks",
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        hwe = hwequilibrium,
        w = windowSize,
    wildcard_constraints:
        chunk="\d{1,2}_\d*_\d*"

    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/run_QTL_analysis.py
        python Limix_QTL/run_QTL_analysis.py \
            --bgen {params.gen} \
            -af {input.af} \
            -pf {input.pf} \
            -cf {input.cf} \
            -od {params.od} \
            -rf {input.kf} \
            --sample_mapping_file {input.smf} \
            -gr {params.chunkstr} \
            -np {params.np} \
            -maf {params.maf} \
            -hwe {params.hwe} \
            -w {params.w} \
            -c -gm gaussnorm -bs 500 -rs 0.95
        """

def collect_qtl_result_files(wildcards):
    checkpoint_output = checkpoints.setup_chunk_analysis.get(**wildcards).output[0]
    wc=glob_wildcards(os.path.join(checkpoint_output, "{chunk}_info"))
    return expand(checkpoint_output + "/qtl_results_{chunk}.h5",
            chunk=wc.chunk)

rule aggregate_qtl_results:
    input:
        collect_qtl_result_files,
    output:
        finalQTLRun
    params:
        IF = outputFolder + "/chunks",
        OF = outputFolder + "/"
    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
        python Limix_QTL/post_processing/minimal_postprocess.py \
            -id {params.IF} \
            -od {params.OF} \
            -sfo
        """

rule top_feature:
    input:
        finalFile = finalQTLRun
    output:
        topQTL
    params:
        IF = outputFolder + "/chunks",
        OF = outputFolder + "/"
    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
        python Limix_QTL/post_processing/minimal_postprocess.py \
            -id {params.IF} \
            -od {params.OF} \
            -tfb \
            -sfo
        """

