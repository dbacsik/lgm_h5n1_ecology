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
        tree = 'output/tree/tree.treefile'

"""This rule loads data from NCBI into nextstrain format"""
rule load_context:
    message: "Loading metadata"
    input:
        fasta = "input/data/NCBI_H5_HA.fasta"
    output:
        metadata = 'output/data_cleaning/context.tsv',
        fasta = 'output/data_cleaning/context.fasta'
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
rule sample_context:
    message: "Subsampling sequences"
    input:
        sequences = rules.load_context.output.fasta,
        metadata = rules.load_context.output.metadata,
        include = 'input/data/SAmerica_Accessions.txt',
        exclude = 'input/data/exclude.txt'
    output:
        fasta = 'output/data_cleaning/context_subsampled.fasta',
        metadata = 'output/data_cleaning/context_subsampled.tsv'
    params:
        groups = 'country year',
        num_sequences = 5000,
        min_length = 1500
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --group-by {params.groups} \
            --subsample-max-sequences {params.num_sequences} \
            --include {input.include} \
            --exclude {input.exclude} \
            --output-sequences {output.fasta} \
            --output-metadata {output.metadata} \
        """

"""This rule loads data from local sequencing into nextstrain format"""
rule load_local:
    message: "Loading metadata"
    input:
        fasta = "input/data/Peru_H5_HA.fasta"
    output:
        metadata = 'output/data_cleaning/peru.tsv',
        fasta = 'output/data_cleaning/peru.fasta'
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
        context = rules.sample_context.output.fasta,
        local = rules.load_local.output.fasta
    output:
        merged = 'output/data_cleaning/merged.fasta'
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
        context = rules.sample_context.output.metadata,
        local = rules.load_local.output.metadata
    output:
        metadata = 'output/data_cleaning/metadata_merged.tsv'
    shell:
        """
        augur merge \
            --metadata \
                REFERENCE={input.context} \
                LOCAL={input.local} \
            --output-metadata {output.metadata}
        """

## ALIGN
rule align:
    message: "Aligning lineage reference sequences with local sequences"
    input:
        fasta = rules.cat_fastas.output.merged,
        reference_file = 'input/reference/LC730539.fasta'
    output:
        alignment = 'output/tree/aligned.fasta'
    shell:
        """
        augur align \
            --sequences {input.fasta} \
            --reference-sequence {input.reference_file} \
            --output {output.alignment} \
            --nthreads 12
        """

### RAW TREE
"""This rule builds the initial tree."""
rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = 'output/tree/tree.treefile'
    shell:
        """
        iqtree -s {input.alignment} \
            -m GTR+G \
            -nt 12 \
            -pre output/tree/tree
        """