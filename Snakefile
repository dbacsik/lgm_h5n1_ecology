# LGM_H5N1
"""
This Snakefile organizes the project. It uses Nextstrain for many steps 
and it is based on the one in the 'quickstart' repository from the public avian-flu 
build: https://github.com/nextstrain/avian-flu/blob/master/quickstart-build/Snakefile
"""

# VARIABLES

# RULES

"""This rule uses a custom R script to clean our input metadata and prepare
it for input into the augur pipeline."""
rule clean_metadata:
    message: "Cleaning metadata"
    input:
        metadata = ['data/raw/gisaid/gisaid_H5N1_None_to_20210101.xls',
                    'data/raw/gisaid/gisaid_H5N1_20220101_to_None.xls']
    output:
        clean_metadata = 'data/raw/gisaid/metadata.tsv'
    shell:
        """
        Rscript data/raw/gisaid/clean_data_gisaid.R
        """
