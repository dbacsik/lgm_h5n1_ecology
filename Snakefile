# LGM_H5N1
"""
This Snakefile organizes the project. It uses Nextstrain for many steps 
and it is based on the one in the 'quickstart' repository from the public avian-flu 
build: https://github.com/nextstrain/avian-flu/blob/master/quickstart-build/Snakefile
"""

# VARIABLES

# RULES

rule all:
    input:
        subsampled_fasta = 'results/context.fasta',
        subsampled_metadata = 'results/context_metadata.tsv'

"""This rule loads data from NCBI into nextstrain format"""
rule load_context:
    message: "Loading metadata"
    input:
        fasta = "data/NCBI_H5_HA.fasta"
    output:
        metadata = 'results/metadata.tsv',
        fasta = 'results/sequences.fasta'
    params:
        fields = 'name title country host date'
    shell:
        """
        augur parse \
            --sequences {input.fasta} \
            --output-sequences {output.fasta} \
            --output-metadata {output.metadata} \
            --fields {params.fields} \
            --separator '|' \
            --fix-dates monthfirst
        """

"""Subsample"""
rule subsample_context:
    message: "Subsampling sequences"
    input:
        sequences = rules.load_context.output.fasta,
        metadata = rules.load_context.output.metadata,
        include = 'data/SAmerica_Accessions.txt'
    output:
        subsampled_fasta = 'results/context.fasta',
        subsampled_metadata = 'results/context_metadata.tsv'
    params:
        groups = 'country year',
        num_sequences = 5000
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --group-by {params.groups} \
            --subsample-max-sequences {params.num_sequences} \
            --include {input.include} \
            --output-sequences {output.subsampled_fasta} \
            --output-metadata {output.subsampled_metadata} \
        """