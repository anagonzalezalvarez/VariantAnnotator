# STEPS FOR STRUCTURAL CHARACTERIZATION OF PATHOGENIC AND BENIGN VARIANTS

## `00_DataPreparation.ipynb`

1. Read output files from vep  
    clinvar = ../results/pipeline_results/clinvar_variants.pkl
    pavs = ../results/pipeline_results/pavs_variants.pkl
    gnomad = "../results/pipeline_results/gnomad_variants.pkl"
2. Concatenate 3 df 
3. Normalize accesibility column. RSA = Area divided by residue size 
4. Save ../results/1_merged.csv

## `01_AssignmentWrapperPhysicochem.py`

1. Assign 
    -Energy
    -Contacts
    -Conservation
    -Order
    -Exposure
    -SS
    -Physico change

2. Save as ../results/2_features_all.pkl
3. Merge with panther protein class `../data/panther/PTHR17.0_human_clean.csv'`
4. Save as ../results/3_features_byprotienclass.csv
5. Filter for only genes in both.
6. Save as '../results/4_features_geneboth.csv'

## `O2_StatisticalTest`

1. Load 2_features_all,3_features_byprotienclass,4_features_geneboth
2. Get equal number of common and pathogenic variants for each protein class
3. ODS individual
4. HEATMAP ODS by protein class
