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
        auspice_json = 'auspice/global.json',
        focal_json = 'auspice/focal.json',

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
        exclude = 'input/data/exclude.txt'
    output:
        fasta = 'output/data_cleaning/context_subsampled.fasta',
        metadata = 'output/data_cleaning/context_subsampled.tsv'
    params:
        groups = 'year',
        num_sequences = 200,
        min_length = 1700,
        min_date = 1995
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --min-date {params.min_date} \
            --exclude {input.exclude} \
            --subsample-seed 111 \
            --subsample-max-sequences {params.num_sequences} \
            --group-by {params.groups} \
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
        fields = 'name country date'
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
                Reference={input.context} \
                Novel={input.local} \
            --output-metadata {output.metadata} \
            --source-columns "source_{{NAME}}"
        """

### Manually parse selected sequences using Nextclade to get 2.3.4.4 clade assignments.
"""This rule adds clades assignments for 2.3.4.4 lineage seqs"""
rule merge_clades:
    message:
        """
        Merging the metadata for the context and the local files.\n
        Metadata: {input.meta}\n
        Clades: {input.clades}
        """
    input:
        meta = rules.merge_metadata.output.metadata,
        clades = 'input/data/clades.tsv'
    output:
        metadata = 'output/data_cleaning/metadata_labeled.tsv'
    shell:
        """
        augur merge \
            --metadata \
                Local={input.meta} \
                Clades={input.clades} \
            --output-metadata {output.metadata} \
        """

## ALIGN
rule align:
    message: "Aligning reference sequences with local sequences"
    input:
        fasta = rules.cat_fastas.output.merged,
        reference_file = 'input/reference/reference_h5n1_ha.fasta'
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

## BUILD TREE

### RAW TREE
"""This rule builds the initial tree."""
rule tree:
    message: "Building lineage tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = 'output/tree/tree.nwk'
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 12
        """

### Refine
"""This rule refines the tree. In this case, we basically
just use it to assign internal node names, since we aren't
going to assign time to the tree."""
rule refine:
    message:
        """
        Refining tree and assigning internal node names.
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.merge_clades.output.metadata
    output:
        tree = 'output/tree/tree_refined.nwk',
        node_data = 'output/tree/node_data.json'
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-filter-iqd 3 \
            --clock-rate 0.00513
        """

### NODES AND TRAITS
"""This rule reconstructs the ancestral node sequences.
It then annotates nucleotide mutations at each branch."""
rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
        reference = 'input/reference/reference_h5n1_ha.fasta'
    output:
        nt_muts = 'output/tree/nt_muts.json',
        sequences = 'output/tree/ancestral_sequences.fasta'
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.nt_muts} \
            --output-sequences {output.sequences} \
            --root-sequence {input.reference}
        """

"""This rule simply translates the sequence of each gene at each node,
including inferred ancestral nodes."""
rule translate:
    message: "Translating amino acid sequences and identifying mutations"
    input:
        tree = rules.refine.output.tree,
        ancestral_json = rules.ancestral.output.nt_muts,
        reference = 'input/reference/reference_h5n1_ha.gb'
    output:
        aa_muts = 'output/tree/aa_muts.json',
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.ancestral_json} \
            --reference-sequence {input.reference} \
            --output-node-data {output.aa_muts}
        """

## VISUALIZATION
"""This rule exports the results of the pipeline into JSON
for visualization in auspice."""
rule export:
    message: "Exporting JSON files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.merge_clades.output.metadata,
        node_data = [rules.refine.output.node_data,
                     rules.ancestral.output.nt_muts,
                     rules.translate.output.aa_muts],
        auspice_config = 'input/config/auspice_config.json',
        colors = 'input/config/colors.tsv'
    output:
        auspice_json = 'auspice/global.json'
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --output {output.auspice_json} \
            --include-root-sequence-inline
        """

##################
# Focal tree

"""Subsample"""
rule sample_focal:
    message: "Subsampling focal sequences"
    input:
        sequences = rules.load_context.output.fasta,
        metadata = rules.load_context.output.metadata,
        include = 'input/data/SAmerica_Accessions.txt',
        exclude = 'input/data/exclude.txt'
    output:
        fasta = 'output/data_cleaning/context_subsampled.fasta',
        metadata = 'output/data_cleaning/context_subsampled.tsv'
    params:
        groups = 'year',
        num_sequences = 200,
        min_length = 1500,
        min_date = 2010
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --min-date {params.min_date} \
            --exclude {input.exclude} \
            --include {input.include} \
            --subsample-seed 123 \
            --subsample-max-sequences {params.num_sequences} \
            --group-by {params.groups} \
            --output-sequences {output.fasta} \
            --output-metadata {output.metadata} \
        """

"""This rule concatenates the local sequences with the reference dataset."""
rule cat_focal:
    message:
        """
        Concatenating the lineage reference FASTA and the local FASTA.\n
        Focal FASTA: {input.focal}\n
        Local FASTA: {input.local}
        """
    input:
        focal = rules.sample_focal.output.fasta,
        local = rules.load_local.output.fasta
    output:
        focal_merged = 'output/data_cleaning/focal_merged.fasta'
    shell: 
        """
        cat {input.focal} {input.local} > {output.focal_merged}
        """


"""This rule merges metadata for local and context samples."""
rule merge_focal:
    message:
        """
        Merging the metadata for the context and the local files.\n
        Focal Metadata: {input.focal}\n
        Local Metadata: {input.local}
        """
    input:
        focal = rules.sample_focal.output.metadata,
        local = rules.load_local.output.metadata
    output:
        metadata = 'output/data_cleaning/focal_merged.tsv'
    shell:
        """
        augur merge \
            --metadata \
                Focal={input.focal} \
                Novel={input.local} \
            --output-metadata {output.metadata}
        """

### Here I take the output of merge_focal and manually annotate novel genomes.
"""This rule merges metadata for local and context samples."""
rule merge_countries:
    message:
        """
        Merging the metadata for the context and the local files.\n
        Focal Metadata: {input.focal}\n
        Labeled Metadata: {input.labeled}
        """
    input:
        focal = rules.merge_focal.output.metadata,
        labeled = 'input/data/focal_labels.tsv'
    output:
        metadata = 'output/data_cleaning/focal_labeled.tsv'
    shell:
        """
        augur merge \
            --metadata \
                Focal={input.focal} \
                Labeled={input.labeled} \
            --output-metadata {output.metadata}
        """

## ALIGN
rule align_focal:
    message: "Aligning focal sequences with local sequences"
    input:
        fasta = rules.cat_focal.output.focal_merged,
        reference_file = 'input/reference/reference_h5n1_ha.fasta'
    output:
        alignment = 'output/tree/focal_aligned.fasta'
    shell:
        """
        augur align \
            --sequences {input.fasta} \
            --reference-sequence {input.reference_file} \
            --output {output.alignment} \
            --nthreads 12
        """

## BUILD TREE

### RAW TREE
"""This rule builds the initial tree."""
rule tree_focal:
    message: "Building focal tree"
    input:
        alignment = rules.align_focal.output.alignment
    output:
        tree = 'output/tree/focal_tree.nwk'
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 12
        """

### Refine
"""This rule refines the tree. In this case, we basically
just use it to assign internal node names, since we aren't
going to assign time to the tree."""
rule refine_focal:
    message:
        """
        Refining focal tree and assigning internal node names.
        """
    input:
        tree = rules.tree_focal.output.tree,
        alignment = rules.align_focal.output.alignment,
        metadata = rules.merge_countries.output.metadata
    output:
        tree = 'output/tree/focal_tree_refined.nwk',
        node_data = 'output/tree/focal_node_data.json'
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate .00513
        """

### NODES AND TRAITS
"""This rule reconstructs the ancestral node sequences.
It then annotates nucleotide mutations at each branch."""
rule ancestral_focal:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine_focal.output.tree,
        alignment = rules.align_focal.output.alignment,
        reference = 'input/reference/reference_h5n1_ha.fasta'
    output:
        nt_muts = 'output/tree/focal_nt_muts.json',
        sequences = 'output/tree/focal_ancestral_sequences.fasta'
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.nt_muts} \
            --output-sequences {output.sequences} \
            --root-sequence {input.reference}
        """

"""This rule simply translates the sequence of each gene at each node,
including inferred ancestral nodes."""
rule translate_focal:
    message: "Translating amino acid sequences and identifying mutations"
    input:
        tree = rules.refine_focal.output.tree,
        ancestral_json = rules.ancestral_focal.output.nt_muts,
        reference = 'input/reference/reference_h5n1_ha.gb'
    output:
        aa_muts = 'output/tree/focal_aa_muts.json',
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.ancestral_json} \
            --reference-sequence {input.reference} \
            --output-node-data {output.aa_muts}
        """

## VISUALIZATION
"""This rule exports the results of the pipeline into JSON
for visualization in auspice."""
rule export_focal:
    message: "Exporting JSON files for for auspice"
    input:
        tree = rules.refine_focal.output.tree,
        metadata = rules.merge_countries.output.metadata,
        node_data = [rules.refine_focal.output.node_data,
                     rules.ancestral_focal.output.nt_muts,
                     rules.translate_focal.output.aa_muts],
        auspice_config = 'input/config/auspice_config.json',
        colors = 'input/config/colors.tsv'
    output:
        auspice_json = 'auspice/focal.json'
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --output {output.auspice_json} \
            --include-root-sequence-inline
        """