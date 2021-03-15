#!/usr/bin/env python3
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
salmon_container = 'docker://combinelab/salmon:latest'
bioconductor_container = 'shub://TomHarrop/r-containers:bioconductor_3.11'

#########
# RULES #
#########

rule target:
    input:
        expand('output/asw_salmon/{sample}_quant/quant.sf', sample=all_samples),
        expand('output/asw_mh_concat_salmon/{sample}_quant/quant.sf', sample=all_samples),
        'output/fastqc',
        #'output/deseq2/asw_dual/unann/nr_blastx.outfmt3',
        'output/deseq2/asw_dual/asw_dual_dds.rds',
        'output/deseq2/mh_dual/mh_dual_dds.rds',
        expand('output/joined/{sample}_r1.fq.gz', sample=all_samples)

#####################
## RNAseq analysis ##
#####################

rule unann_degs_blastx:
    input:
        unann_deg_transcripts = 'output/deseq2/asw_dual/unann/unann_deg_transcripts.fasta'
    output:
        blastx_res = 'output/deseq2/asw_dual/unann/nr_blastx.outfmt3'
    params:
        blast_db = 'bin/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/unann_degs_blastx.log'
    shell:
        'blastx '
        '-query {input.unann_deg_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_unann_deg_transcripts:
    input:
        dual_transcriptome = 'data/asw_mh_transcriptome/asw_mh_isoforms_by_length.fasta',
        transcript_hit_ids = 'output/deseq2/asw_dual/unann/unann_degs_list.txt'
    output:
        unann_deg_transcripts = 'output/deseq2/asw_dual/unann/unann_deg_transcripts.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_unann_deg_transcripts.log'
    shell:
        'filterbyname.sh '
        'in={input.dual_transcriptome} '
        'include=t '
        'names={input.transcript_hit_ids} '
        'substring=name '
        'out={output.unann_deg_transcripts} '
        '&> {log}'

rule ID_unann_DEGs_dual:
    input:
        asw_dds = 'output/deseq2/asw_dual/asw_dual_dds.rds',
        loc_ex_int_degs = 'output/deseq2/asw_dual/location_exposure_int/sig_w_annots.csv',
        loc_degs = 'output/deseq2/asw_dual/location_pairwise/sig_w_annots.csv'
    output:
        unann_degs_list = 'output/deseq2/asw_dual/unann/unann_degs_list.txt'
    singularity:
        bioconductor_container
    log:
        'output/logs/ID_unann_DEGs_dual.log'
    script:
        'src/dual_species/asw/filter_unann_degs.R'

########################################
## map to asw-mh concat transcriptome ##
########################################

rule mh_dual_dds:
    input:
        mh_gene_trans_map = 'data/asw-mh-combined-transcriptome/output/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map',
        quant_files = expand('output/asw_mh_concat_salmon/{sample}_quant/quant.sf', sample=all_samples)
    output:
        mh_dds = 'output/deseq2/mh_dual/mh_dual_dds.rds'
    singularity:
        bioconductor_container
    log:
        'output/logs/mh_dual_dds.log'
    script:
        'src/dual_species/make_mh_dds.R'

rule asw_dual_dds:
    input:
        asw_gene_trans_map = 'data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/Trinity.fasta.gene_trans_map',
        quant_files = expand('output/asw_mh_concat_salmon/{sample}_quant/quant.sf', sample=all_samples)
    output:
        asw_dds = 'output/deseq2/asw_dual/asw_dual_dds.rds'
    singularity:
        bioconductor_container
    log:
        'output/logs/asw_dual_dds.log'
    script:
        'src/dual_species/make_asw_dds.R'

rule asw_mh_concat_salmon_quant:
    input:
        index_output = 'output/asw_mh_concat_salmon/transcripts_index/refseq.bin',
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
        transcriptome_length_filtered = 'data/asw-mh-combined-transcriptome/output/asw_mh_transcriptome/asw_mh_isoforms_by_length.fasta'
    output:
        'output/asw_mh_concat_salmon/transcripts_index/refseq.bin'
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
        index_output = 'output/asw_salmon/transcripts_index/refseq.bin',
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
        'output/asw_salmon/transcripts_index/refseq.bin'
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

######################
## prep for mapping ##
######################

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

##zcat for compressed files
rule join_reads:
    input:
        unpack(get_reads)
    output:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz',
    shell:
        'cat {input.l1r1} {input.l2r1} > {output.r1} & '
        'cat {input.l1r2} {input.l2r2} > {output.r2} & '
        'wait'

##OG file with sequencing for this project in different structure - two folders
##folder with small files is a small miseq run
##folder with large files are standard hiseq files to be used for analysis