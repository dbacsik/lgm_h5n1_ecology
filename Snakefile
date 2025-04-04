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
        merged = 'results/sequences_merged.fasta',
        metadata = 'results/metadata_merged.tsv'

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
        subsampled_fasta = 'results/context_subsampled.fasta',
        subsampled_metadata = 'results/context_subsampled.tsv'
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

"""This rule loads data from local sequencing into nextstrain format"""
rule load_local:
    message: "Loading metadata"
    input:
        fasta = "data/Peru_H5_HA.fasta"
    output:
        metadata = 'results/peru.tsv',
        fasta = 'results/peru.fasta'
    params:
        fields = 'name country'
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

"""This rule concatenates the local sequences with the reference dataset."""
rule cat_fastas:
    message:
        """
        Concatenating the lineage reference FASTA and the local FASTA.\n
        Reference FASTA: {input.context}\n
        Local FASTA: {input.local}
        """
    input:
        context = rules.load_context.output.fasta,
        local = rules.load_local.output.fasta
    output:
        merged = 'results/sequences_merged.fasta'
    shell: 
        """
        cat {input.context} {input.local} > {output.merged}
        """

"""This rule merges metadata for local and context samples."""
rule merge_metadata:
    message:
        """
        Merging the metadata for the context and the local files.\n
        Context FASTA: {input.context}\n
        Local FASTA: {input.local}
        """
    input:
        context = rules.load_context.output.metadata,
        local = rules.load_local.output.metadata
    output:
        metadata = 'results/metadata_merged.tsv'
    shell:
        """
        augur merge \
            --metadata \
                REFERENCE={input.context} \
                LOCAL={input.local} \
            --output-metadata {output.metadata}
        """