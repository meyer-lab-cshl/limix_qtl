import glob
import os
from subprocess import run
import pandas as pd
import re
from os.path import join

shell.prefix("set -euo pipefail;")

def _multi_arg_start(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)

def _multi_arg_end(flag, files):
    flag = " "+flag
    return " ".join(f + flag for f in files)

def _multi_arg_both_ends(flag1, flag2, files):
    flag1 += " "
    flag2 = " "+flag2
    return " ".join(flag1 + f + flag2 for f in files)

def flatenChunk(chunk):
    return chunk.replace(":", "_").replace("-", "_")

def extendChunk(chunk):
    relChunk = chunk.pop()
    chunkSplitted = relChunk.split("_")
    return chunkSplitted[0]+":"+chunkSplitted[1]+"-"+chunkSplitted[2]

## Variables ##
## Files, these are populated now with the test examples and we use here plink genotypes.
chunkFile =  "Limix_QTL/test_data/Expression/chunks_Geuvadis_CEU_Annot.txt"
genotypeFile = 'Limix_QTL/test_data/Genotypes/Geuvadis'  ##Genotype without file extension. Please update flag in the runner to reflect --plink or --bgen
kinshipFile = 'Limix_QTL/test_data/Genotypes/Geuvadis_chr1_kinship.normalized.txt'
annotationFile = 'Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot.txt'
phenotypeFile = 'Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_Expr.txt'
covariateFile = 'Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_covariates.txt'
randomEffFile = '' #no second random effect in this example
sampleMappingFile = 'Limix_QTL/test_data/Geuvadis_CEU_gte.txt'
outputFolder = 'OutGeuvadis/'

## Settings
nGenes = 50
startPos = 0
endOffset = 1000000000
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
windowSize = '1000000'
hwequilibrium = '0.000001'
FDR = '0.05'
## End Variables ##

finalQTLRun = outputFolder+'qtl_results_all.txt'
topQTL = outputFolder+'top_qtl_results_all.txt'


if False:
    with open(chunkFile,'r') as f:
        chunks = [x.strip() for x in f.readlines()]

    qtlOutput = []
    for chunk in chunks:
        #print(chunk)
        processedChunk = flatenChunk(chunk)
        #print(processedChunk)
        processedChunk=expand(outputFolder+'{chunk}.finished',chunk=processedChunk )
        #print(processedChunk)
        qtlOutput.append(processedChunk)

    ## flatten these lists
    qtlOutput = [filename for elem in qtlOutput for filename in elem]

rule all:
    input:
        finalQTLRun, topQTL


rule generate_chunks:
    input:
        annotation=annotationFile
    output:
        chunks=chunkFile
    params:
        nGenes = nGenes,
        startPos = startPos,
        endOffset = endOffset
    script:
        "Limix_QTL/scripts/generate_chunks.R"


checkpoint run_qtl_mapping:
    input:
        af = annotationFile,
        pf = phenotypeFile,
        cf = covariateFile,
        kf = kinshipFile,
        smf = sampleMappingFile
    output:
        directory(outputFolder/{chunk})
    params:
        gen=genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        hwe = hwequilibrium,
        w = windowSize,
    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
        python ./Limix_QTL/run_QTL_analysis.py
            --bgen {params.gen} \
            -af {input.af} \
            -pf {input.pf} \
            -cf {input.cf} \
            -od {params.od} \
            -rf {input.kf} \
            --sample_mapping_file {input.smf} \
            -gr {wildcards.chunk} \
            -np {params.np} \
            -maf {params.maf} \
            -hwe {params.hwe} \
            -w {params.w} \
            -c -gm gaussnorm -bs 500 -rs 0.95
        """

rule aggregate_qtl_results:
    input:
        #finalFiles = qtlOutput
        finalFiles = expand(outputFolder + '{chunk}.finished',
                            chunk=[flatenChunk(c) for c in chunks])
    params:
        IF = outputFolder,
        OF = outputFolder,
    output:
        finalQTLRun
    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
        python ./Limix_QTL/post_processing/minimal_postprocess.py \
            -id {params.IF} \
            -od {params.OF} \
            -sfo
        """

rule top_feature:
    input:
        #IF = outputFolder,
        #OF = outputFolder,
        finalFile = finalQTLRun
    params:
        IF = outputFolder,
        OF = outputFolder,
    output:
        topQTL
    run:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
        python ./Limix_QTL/post_processing/minimal_postprocess.py \
            -id {params.IF} \
            -od {params.OF} \
            -tfb \
            -sfo
        """

