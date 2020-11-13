#!/usr/bin/env python3
import pathlib2
import os
import pandas
import peppy

#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    input_keys = ['l1r1', 'l2r1', 'l1r2', 'l2r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
tom_salmon_container = 'shub://TomHarrop/singularity-containers:salmon_0.11.1'
##try and see whether salmon works off docker
salmon_container = 'docker://combinelab/salmon:latest'

#########
# SETUP #
#########

#sample_key = pandas.read_csv(sample_key_file)

#all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/asw_salmon/{sample}_quant/quant.sf', sample=all_samples),
        expand('output/asw_mh_concat_salmon/{sample}_quant/quant.sf', sample=all_samples),
        'output/fastqc',

############################
## map to asw/mh combined ##
############################

rule asw_mh_concat_salmon_quant:
    input:
        index_output = 'output/asw_mh_concat_salmon/transcripts_index/hash.bin',
        left = 'output/bbduk_trim/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/asw_mh_concat_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/asw_mh_concat_salmon/transcripts_index',
        outdir = 'output/asw_mh_concat_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon/asw_mh_concat_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.left} '
        '-2 {input.right} '
        '-o {params.outdir} '
        '--writeUnmappedNames '
        '-p {threads} '
        '&> {log}'

rule asw_mh_concat_salmon_index:
    input:
        transcriptome_length_filtered = 'data/asw_mh_transcriptome/asw_mh_isoforms_by_length.fasta'
    output:
        'output/asw_mh_concat_salmon/transcripts_index/hash.bin'
    params:
        outdir = 'output/asw_mh_concat_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/asw_mh_concat_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

##############################
## map to asw transcriptome ##
##############################

rule asw_salmon_quant:
    input:
        index_output = 'output/asw_salmon/transcripts_index/hash.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/asw_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/asw_salmon/transcripts_index',
        outdir = 'output/asw_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon/asw_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.trimmed_r1} '
        '-2 {input.trimmed_r2} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule asw_salmon_index:
    input:
        transcriptome_length_filtered = 'data/asw-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta'
    output:
        'output/asw_salmon/transcripts_index/hash.bin'
    params:
        outdir = 'output/asw_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/asw_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule fastqc:
    input:
        expand('output/bbduk_trim/{sample}_r{n}.fq.gz',
            sample=all_samples, n=[1,2])
    output:
        directory('output/fastqc')
    shell:
        'mkdir -p {output} ; '
        'fastqc --outdir {output} {input}'

rule bbduk_trim:
    input:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/bbduk_trim/{sample}.log'
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

rule join_reads:
    input:
        unpack(get_reads)
    output:
        r1 = temp('output/joined/{sample}_r1.fq.gz'),
        r2 = temp('output/joined/{sample}_r2.fq.gz'),
    shell:
        'zcat {input.l1r1} {input.l2r1} >> {output.r1} & '
        'zcat {input.l1r2} {input.l2r2} >> {output.r2} & '
        'wait'

##OG file with sequencing for this project in different structure - two folders
##folder with small files is a small miseq run
##folder with large files are standard hiseq files to be used for analysis