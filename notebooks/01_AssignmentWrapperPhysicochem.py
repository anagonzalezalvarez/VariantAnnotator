

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import os

try:
    os.chdir(os.path.join(os.getcwd(), './notebooks')) # '.' if the path is to current folder
except:
    pass

def getting_ref_aa_properties(path_variants_input,path_results, AA_GIVEN = False):
    """
    Function to assign the properties treshold to the variants input file

    parameters
    ------------------
    path_variants_input: str
        path to the variants input file
    path_results: str
        path to the results file
    AA_GIVEN: bool
        if the aminoacid ref is in column REFAA and aa alt in column ALTAA is given in the input file as one letter

    
    """
    variants = pd.read_pickle(path_variants_input)

    if AA_GIVEN == False:
        variants[['REFAA', 'pos','ALTAA']] = variants['HGVS_p'].str.extract(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})')
        aa_dict = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}

        # Load data
        variants['ALTAA'] = variants['ALTAA'].map(aa_dict)
        variants['REFAA'] = variants['REFAA'].map(aa_dict)

    ## ENERGY
    variants['Stabilizing Energy Change'] = ((variants['total energy']).apply(lambda x: 'Stabilizing Energy Change' if x<0 else 'not_Stabilizing Energy Change'))
    variants['Destabilizing Energy Change (<1 kcal/mol)'] = ((variants['total energy']).apply(lambda x: 'Destabilizing Energy Change (<1 kcal/mol)' if x > 0 and x < 1 else 'not_Destabilizing Energy Change (<1 kcal/mol)'))
    variants['Destabilizing Energy Change (1-3 kcal/mol)'] = ((variants['total energy']).apply(lambda x: 'Destabilizing Energy Change (1-3 kcal/mol)' if x > 1 and x < 3 else 'not_Destabilizing Energy Change (1-3 kcal/mol)'))
    variants['Destabilizing Energy Change (>3 kcal/mol)'] = ((variants['total energy']).apply(lambda x: 'Destabilizing Energy Change (>3 kcal/mol)' if x>3 else 'not_Destabilizing Energy Change (>3 kcal/mol)'))

    variants['TotalEnergy'] = ((variants['total energy']).apply(lambda x: 'TotalEnergy' if x>3 else 'not_TotalEnergy'))
    variants['LessTotalEnergy'] = ((variants['total energy']).apply(lambda x: 'LessTotalEnergy' if x<0 else 'not_LessTotalEnergy'))
    variants['VanDerWaalsClashes'] = ((variants['Van der Waals clashes']).apply(lambda x: 'VanDerWaalsClashes' if x>0.9 or x < -0.003 else 'not_VanDerWaalsClashes'))
    variants['Disulfide'] = ((variants['disulfide']).apply(lambda x: 'Disulfide' if x>0.1 or x < 0 else 'not_Disulfide'))
   
    ## CONTACTS
    variants['Contacts'] = ((variants['intra_contacts']).apply(lambda x: 'Contacts' if len(x) != 0 else 'not_Contacts'))
    variants['No Contacts'] = ((variants['intra_contacts']).apply(lambda x: 'No Contacts' if len(x) == 0 else 'not_No Contacts'))


    ## CONSERVATION
    variants['Conserved'] = ((variants['Conservation']).apply(lambda x: 'Conserved' if x > 2 else 'not_Conserved'))

    ## CATALYTIC/DOMAIN NOT USED
    variants['Catalytic'] = ((variants['is_catalytic']).apply(lambda x: 'Catalytic' if x == True else 'not_Catalytic'))
    variants['Domain'] = ((variants['DOMAINS']).apply(lambda x: 'Catalytic' if x == True else 'not_Domain'))

    ## ORDER
    variants['DisorderpLDDT'] = ((variants['pLDDT']).apply(lambda x: 'DisorderpLDDT' if x < 50 else 'not_DisorderpLDDT'))
    variants['OrderpLDDT'] = ((variants['pLDDT']).apply(lambda x: 'OrderpLDDT' if x > 50 else 'not_OrderpLDDT'))

    ## EXPOSURE
    variants['Core (<5%)'] = ((variants['ACCESIBILITY_NORMALIZED']).apply(lambda x: 'Core (<5%)' if x < 0.05 else 'not_Core (<5%)'))
    variants['Buried (5-25%)'] = ((variants['ACCESIBILITY_NORMALIZED']).apply(lambda x: 'Buried (5-25%)' if x > 0.05 and x < 0.25 else 'not_Buried (5-25%)'))
    variants['Medium-buried (25-50%)'] = ((variants['ACCESIBILITY_NORMALIZED']).apply(lambda x: 'Medium-buried (25-50%)' if x > 0.25 and x < 0.5 else 'not_Medium-buried (25-50%)'))
    variants['Medium-exposed (50-75%)'] = ((variants['ACCESIBILITY_NORMALIZED']).apply(lambda x: 'Medium-exposed (50-75%)' if x > 0.5 and x < 0.75 else 'not_Medium-exposed (50-75%)'))
    variants['Exposed (>75%)'] = ((variants['ACCESIBILITY_NORMALIZED']).apply(lambda x: 'Exposed (>75%)' if x > 0.75 else 'not_Exposed (>75%)'))

    ## SECONDARY STRUCTURE
    helices = ['H','G','I']
    βetas = ['B','E']
    coils = ['L','S','T']  

    variants['helices'] = variants.apply(lambda row: 'helices' if (row['secondary_structure'] in helices) else 'not_helices', axis=1)
    variants['β-sheet/starnd'] = variants.apply(lambda row: 'β-sheet/starnd' if (row['secondary_structure'] in βetas) else 'not_β-sheet/starnd', axis=1)
    variants['coils'] = variants.apply(lambda row: 'coils' if (row['secondary_structure'] in coils) else 'not_coils', axis=1)
    variants['α-helix'] = ((variants['secondary_structure']).apply(lambda x: 'α-helix' if x == "H" else 'not_α-helix'))
    variants['βbridge'] = ((variants['secondary_structure']).apply(lambda x: 'βbridge' if x == "B" else 'not_βbridge'))
    variants['strandβladder'] = ((variants['secondary_structure']).apply(lambda x: 'strandβladder' if x == "E" else 'not_strandβladder'))
    variants['310helix'] = ((variants['secondary_structure']).apply(lambda x: '310helix' if x == "G" else 'not_310helix'))
    variants['π-helix'] = ((variants['secondary_structure']).apply(lambda x: 'π-helix' if x == "I" else 'not_π-helix'))
    variants['HBondTurn'] = ((variants['secondary_structure']).apply(lambda x: 'HBondTurn' if x == "T" else 'not_HBondTurn'))
    variants['Bend'] = ((variants['secondary_structure']).apply(lambda x: 'Bend' if x == "S" else 'not_Bend'))
    variants['Loop'] = ((variants['secondary_structure']).apply(lambda x: 'Loop' if x == "L" else 'not_Loop'))


    ## PHYSICOCHEMICAL CHANGE
    small = ['G','A','S']
    big = ['F','W','Y','R','K','L','I','M']
    nonpolar = ['G','A','V','C','P','L','I','M','W','F']
    positive = ['H','K','R']
    negative = ['D','E']
    charged = positive + negative
    polar_uncharged = ['S','T','Y','N','Q']
    polar = positive + negative + polar_uncharged
    aromatic = ['F','W','Y']
    uncharged = nonpolar + polar_uncharged
    hydrophobic = ['A', 'C', 'I', 'L', 'M', 'F', 'W', 'V']
    hydrophilic = ['R', 'N', 'D', 'Q', 'E', 'K']
    neutral = ['G', 'H', 'P', 'S', 'T', 'Y']

    # create pysicochemical classifications REF aminoacid and ALT aminoacid
    variants['Small to Big'] = variants.apply(lambda row: 'Small to Big' if row['REFAA'] in small and row['ALTAA'] in big else 'not_Small to Big', axis=1)
    variants['Big to Small'] = variants.apply(lambda row: 'Big to Small' if row['REFAA'] in big and row['ALTAA'] in small else 'not_Big to Small', axis=1)
    variants['Polar to NonPolar'] = variants.apply(lambda row: 'Polar to NonPolar' if row['REFAA'] in polar and row['ALTAA'] in nonpolar else 'not_Polar to NonPolar', axis=1)
    variants['NonPolar to Polar'] = variants.apply(lambda row: 'NonPolar to Polar' if row['REFAA'] in nonpolar and row['ALTAA'] in polar else 'not_NonPolar to Polar', axis=1)
    variants['Hydrophilic introduced'] = variants.apply(lambda row: 'Hydrophilic introduced' if (row['REFAA'] in hydrophobic or row['REFAA'] in neutral) and (row['ALTAA'] in hydrophilic) else 'not_Hydrophilic introduced', axis=1)
    variants['Hydrophobic introduced'] = variants.apply(lambda row: 'Hydrophobic introduced' if (row['REFAA'] in hydrophilic or row['REFAA'] in neutral) and (row['ALTAA'] in hydrophobic or row['ALTAA'] in neutral) else 'not_Hydrophobic introduced', axis=1)
    variants['Charge switch'] = variants.apply(lambda row: 'Charge switch' if (row['REFAA'] in positive and row['ALTAA'] in negative) or (row['REFAA'] in negative and row['ALTAA'] in positive) else 'not_Charge switch', axis=1)
    variants['Charge lost'] = variants.apply(lambda row: 'Charge lost' if row['REFAA'] in charged and row['ALTAA'] in uncharged else 'not_Charge lost', axis=1)
    variants['Charge gain'] = variants.apply(lambda row: 'Charge gain' if row['REFAA'] in uncharged and row['ALTAA'] in charged else 'not_Charge gain', axis=1)
    variants['Aromatic to NonAromatic'] = variants.apply(lambda row: 'Aromatic to NonAromatic' if row['REFAA'] in aromatic and row['ALTAA'] not in aromatic else 'not_Aromatic to NonAromatic', axis=1)
    variants['Aromatic to polar'] = variants.apply(lambda row: 'Aromatic to polar' if row['REFAA'] in aromatic and row['ALTAA'] not in polar else 'not_Aromatic to polar', axis=1)
    
    
    # save as pickle
    variants.to_pickle(path_results)




path_variants_input = "../results/1_merged.pkl"
path_results = '../results/2_features_all.pkl'
AA_GIVEN = True
getting_ref_aa_properties(path_variants_input,path_results,AA_GIVEN)

path_restults_protein = '../results/3_features_byprotienclass.csv'
panther_1_nas = pd.read_csv('../data/panther/PTHR17.0_human_clean.csv' )
concatenated_df_physico = pd.read_pickle(path_results)
clinvar_gnomad_class = pd.merge(concatenated_df_physico, panther_1_nas, how='left', left_on='SYMBOL', right_on='Gene')
clinvar_gnomad_class_filt = clinvar_gnomad_class[clinvar_gnomad_class.GeneralProteinClass.notna()]
clinvar_gnomad_class_filt.to_csv(path_restults_protein, index=False)


### subset for only genes in both
pathogenic_rows = clinvar_gnomad_class_filt[clinvar_gnomad_class_filt['Source'] == 'Pathogenic']
common_rows = clinvar_gnomad_class_filt[clinvar_gnomad_class_filt['Source'] == 'Common']

# Find the values in 'SYMBOL' column that exist in both filtered DataFrames
symbols_with_both = set(pathogenic_rows['SYMBOL']).intersection(common_rows['SYMBOL'])
clinvar_gnomad_class_filt_both = clinvar_gnomad_class_filt[clinvar_gnomad_class_filt['SYMBOL'].isin(symbols_with_both)]
path_restults_filter = '../results/4_features_geneboth.csv'
clinvar_gnomad_class_filt_both.to_csv(path_restults_filter)
