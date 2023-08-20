
# ========== Packages ==========
from re import X
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from scipy.stats import levene
import scipy.stats as stats
from scipy.stats import fisher_exact
import numpy as np
import matplotlib.patches as mpatches
from scipy.stats import ttest_ind
from scipy.stats import levene
import scipy.stats as stats
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import matplotlib.font_manager
import seaborn as sns
import math

import matplotlib.path as mpath
import matplotlib.patches as mpatches






# ========== Functions ODS + Plotting==========



def calculating_ODS_individual(feature,df):

    """
    Obtains fisher_exact test for individual feature. 
        1. Split the name of the column to get the pos_feature value (feature present) and neg_feature value (feature not present)
        2. Does a crosstab to get the amount of variants for Pathogenic and Benign that have that feature present or not
        3. Does fisher_exact -> we obtain oddsratio and pvalue 
    Parameters
    ------------------ 
    feature: string
        name of the feature we want to test. In the following format "x_assignment"     
    
    df: pandas.DataFrame
        df that contains a column with the name of the feature in the following format "x_assignment" 
        and the possible values in this column are: "X" or "not_X" 
        features get assigned as x_assignment -> x or not_x. (Ex: disorder_assignment, disorder,not_disorder )
    """

    #pos_feature = feature.split('_')[0]  ## Takes feature (ex: 'disorder_assignment'splits it by '_' and takes first part)
    pos_feature = feature  ## Takes feature (ex: 'disorder_assignment'splits it by '_' and takes first part)
    neg_feature = f'not_{pos_feature}' ## Adds a 'not_' to the pos_feature

    try:
        frequencies_table = pd.crosstab(df[feature], df['Source']).loc[[pos_feature, neg_feature],['Pathogenic','Common']] ## Just for our convinience to be able to visualize table
        frequencies = pd.crosstab(df[feature], df['Source']).loc[[pos_feature, neg_feature],['Pathogenic','Common']].to_numpy()
        #create sum_pat with the value of pos_feature + neg_feature for Pathogenic
        pat_feat = frequencies[0][0] 
        pat_nofeat = frequencies[1][0]
        ben_feat = frequencies[0][1] 
        ben_nofeat = frequencies[1][1]

        oddsratio, pvalue = fisher_exact(frequencies)
    
    except:
        frequencies =  float("nan") ##asign them as nan if it was not possible to do the test (When function calculate_ODS_multiple_grouped_by calls this part it can
                                    ## be the case that the gene doesn't follow the criteria for all the features so just add nan to the features that do not meet criteria )
        oddsratio = float("nan")
        pvalue =  float("nan")
        pat_feat = float("nan")
        pat_nofeat = float("nan")
        ben_feat = float("nan")
        ben_nofeat = float("nan")

    return frequencies,oddsratio,pvalue,pat_feat,pat_nofeat,ben_feat,ben_nofeat

def calculate_ODS_mult(feature_list, df):
    """
    Create a dataframe "sum_results" that contains all the oddsratio(OR), p-value, q-value, lowerconfidence interval(CIL),
    upper confidence interval(CIU),significant and ORs_clasification for all our features of interest.

    Parameters
    ------------------ 
    feature_list: list
        list that contains the name of all the features we want to test. In the following format "X_assignment"     
    
    df: pandas.DataFrame
        df that contains all the features as columns with possible results being: "X" or "not_X" 
        features get assigned as x_assignment -> x or not_x. (Ex: disorder_assignment, disorder,not_disorder )
    """

    sum_results = pd.DataFrame(columns=['OR', 'p-value','CIL','CIU','pat_feat','pat_nofeat','ben_feat','ben_nofeat']) ## create df where we will store our results

    for feat in feature_list:
        wkg_freq, wkg_OR, wkg_pval,pat_feat,pat_nofeat,ben_feat,ben_nofeat = calculating_ODS_individual(feat, df) ## for each of our features we call the function to calculate individual ODS and pvalue calculating_ODS_individual
    
        #if np.logical_not(wkg_pval).any(): ##check if is not na, then we can continue 
        #if not math.isnan(wkg_OR): ##check if is not na, then we can continue
        try:
            wkg_ci_lower = np.exp(np.log(wkg_OR) - 1.96*np.sqrt(np.sum(1 / np.array(wkg_freq).flatten())))

            wkg_ci_upper = np.exp(np.log(wkg_OR) + 1.96*np.sqrt(np.sum(1 / np.array(wkg_freq).flatten())))

            
            sum_results.loc[feat] = [wkg_OR,wkg_pval,wkg_ci_lower,wkg_ci_upper, pat_feat,pat_nofeat,ben_feat,ben_nofeat] ## store the results in our df

        except:
            wkg_ci_lower = float("nan")
            wkg_ci_upper = float("nan")
            sum_results.loc[feat] = [wkg_OR,wkg_pval,wkg_ci_lower,wkg_ci_upper, pat_feat,pat_nofeat,ben_feat,ben_nofeat] ## store the results in our df  

    sum_results.insert(2, 'q-value', sum_results['p-value']*len(feature_list)) ##calculate qvalues
    sum_results['significant?'] = sum_results['q-value'].apply(
                                                    lambda x: 'significant' if x < 0.05 else 'not significant' ) ## check if its significant
    sum_results['ORs_clasification'] = sum_results['OR'].apply(
                                                    lambda x: 'pathogenic' if x > 1 else ('benign' if x < 1 else 'not_sig')) ## check if its associated with ben or pat


    return sum_results  


def calculate_ODS_multiple(feature_list, df, groupby):
    """
    Create a dataframe "sum_results" that contains all the oddsratio(OR), p-value, q-value, lowerconfidence interval(CIL),
    upper confidence interval(CIU),significant and ORs_clasification for all our features of interest.

    Parameters
    ------------------ 
    feature_list: list
        list that contains the name of all the features we want to test. In the following format "X_assignment"     
    
    df: pandas.DataFrame
        df that contains all the features as columns with possible results being: "X" or "not_X" 
        features get assigned as x_assignment -> x or not_x. (Ex: disorder_assignment, disorder,not_disorder )
    groupby: string
        in case it is given we add a column with the string 
    """

    sum_results = pd.DataFrame(columns=['OR', 'p-value','CIL','CIU','pat_feat','pat_nofeat','ben_feat','ben_nofeat']) ## create df where we will store our results

    for feat in feature_list:
        wkg_freq, wkg_OR, wkg_pval,pat_feat,pat_nofeat,ben_feat,ben_nofeat = calculating_ODS_individual(feat, df) ## for each of our features we call the function to calculate individual ODS and pvalue calculating_ODS_individual
    
        #if np.logical_not(wkg_pval).any(): ##check if is not na, then we can continue 
        #if not math.isnan(wkg_OR): ##check if is not na, then we can continue
        try:
            wkg_ci_lower = np.exp(np.log(wkg_OR) - 1.96*np.sqrt(np.sum(1 / np.array(wkg_freq).flatten())))

            wkg_ci_upper = np.exp(np.log(wkg_OR) + 1.96*np.sqrt(np.sum(1 / np.array(wkg_freq).flatten())))

            
            sum_results.loc[feat] = [wkg_OR,wkg_pval,wkg_ci_lower,wkg_ci_upper, pat_feat,pat_nofeat,ben_feat,ben_nofeat] ## store the results in our df

        except:
            wkg_ci_lower = float("nan")
            wkg_ci_upper = float("nan")
            sum_results.loc[feat] = [wkg_OR,wkg_pval,wkg_ci_lower,wkg_ci_upper, pat_feat,pat_nofeat,ben_feat,ben_nofeat] ## store the results in our df  

    sum_results.insert(2, 'q-value', sum_results['p-value']*len(feature_list)) ##calculate qvalues
    sum_results['significant?'] = sum_results['q-value'].apply(
                                                    lambda x: 'significant' if x < 0.05 else 'not significant' ) ## check if its significant
    sum_results['ORs_clasification'] = sum_results['OR'].apply(
                                                    lambda x: 'pathogenic' if x > 1 else ('benign' if x < 1 else 'not_sig')) ## check if its associated with ben or pat


    sum_results.insert(loc=0, column='groupby', value=groupby)  
    return sum_results



def calculate_ODS_multiple_grouped_by(common_elements, groupbyName, feature_list, df):
    """
    Create a dataframe "full_results_df" that contains all the oddsratio(OR), p-value, q-value, lowerconfidence interval(CIL),
    upper confidence interval(CIU),significant and ORs_clasification for all our features of interest group by another column, example by protein class.

    IMPORTANT: to get common elemnts first we should apply get_groups_per_feature and get_common_elemnts function manually, 
    to be able to check if we want to remove any feature, maybe medium_assignment if we didn't get enough results

    Parameters
    ------------------ 
    common_elements: list
        list that contains from the groupbyName column that we are gouping by. The ones that meet criteria 10 variants for 4 combinations YES_FEAT BEN/PAT, NO_FEAT: BEN/PAT
    feature_list: list
        list that contains the name of all the features we want to test. In the following format "X_assignment"     
    
    df: pandas.DataFrame
        df that contains all the features as columns with possible results being: "X" or "not_X" 
        features get assigned as x_assignment -> x or not_x. (Ex: disorder_assignment, disorder,not_disorder )
        also contains columns that can be used as gouped
    
    groupbyName: string
        name of the column according to what i will group by results: CLNDN (disease), GeneSymbol, GeneralProteinClass, etc
    """
    full_results_df = pd.DataFrame({}) ## empty df
    
    full_results_df = calculate_ODS_mult(feature_list, df)
    full_results_df.insert(loc=0, column='groupby', value="All")  


    for c in common_elements:
        wkg_variants_df = df[df[groupbyName]==c]
        wkg_results_df = calculate_ODS_multiple(feature_list, wkg_variants_df, c)

        full_results_df = pd.concat([full_results_df,wkg_results_df])
        #full_results_df = full_results_df.append(wkg_results_df)

  
    return full_results_df








def add_label_band(ax, top, bottom, label, *, spine_pos=-0.45, tip_pos=-0.02,fontsizeinput=8):
    """
    Helper function to add bracket around y-tick labels.

    Parameters
    ----------
    ax : matplotlib.Axes
        The axes to add the bracket to

    top, bottom : floats
        The positions in *data* space to bracket on the y-axis

    label : str
        The label to add to the bracket

    spine_pos, tip_pos : float, optional
        The position in *axes fraction* of the spine and tips of the bracket.
        These will typically be negative

    Returns
    -------
    bracket : matplotlib.patches.PathPatch
        The "bracket" Aritst.  Modify this Artist to change the color etc of
        the bracket from the defaults.

    txt : matplotlib.text.Text
        The label Artist.  Modify this to change the color etc of the label
        from the defaults.

    """
    # grab the yaxis blended transform
    transform = ax.get_yaxis_transform()

    # add the bracket
    bracket = mpatches.PathPatch(
        mpath.Path(
            [
                [tip_pos, top],
                [spine_pos, top],
                [spine_pos, bottom],
                [tip_pos, bottom],
            ]
        ),
        transform=transform,
        clip_on=False,
        facecolor="none",
        edgecolor="k",
        linewidth=2,
    )
    ax.add_artist(bracket)

    # add the label
    txt = ax.text(
        spine_pos,
        (top + bottom) / 2,
        label,
        ha="right",
        va="center",
        #rotation="vertical",
        clip_on=False,
        transform=transform,
        fontsize=fontsizeinput,
        fontweight='bold'
    )

    return bracket, txt




def plot_ODS(test_features):


    """
    Parameters
    ------------------ 
    test_features: pandas.dataframe
        df that contains the oddsratio(OR), p-value, q-value, lowerconfidence interval(CIL),
        upper confidence interval(CIU),significant and ORs_clasification for all our features of interest.
    """
    sns.set_theme(style="whitegrid")
    
    #fig, ax = plt.subplots(figsize=(23, 32))  # Sample figsize in inches

    fig, ax = plt.subplots(nrows=1, sharex=True, sharey=True, dpi=150,figsize=(10, 8))
    for idx, row in test_features.iloc[::-1].iterrows():
        ci = [[row['OR'] - row[::-1]['CIL']], [row['CIU'] - row['OR']]] ## calculate confidence interaval
        
       
        if row['significant?'] == 'significant': ## if its significant we check if ben/pat to decide color
            if row['ORs_clasification'] == 'pathogenic': ## color red if its pat
                plt.errorbar(x=[row['OR']], y=[row.name], xerr=ci,
                    ecolor='tab:gray', capsize=3, linestyle='None', linewidth=1, marker="o", 
                            markersize=5, mfc="tab:red", mec="tab:red")
            elif row['ORs_clasification'] == 'benign': ## color blue if its ben
                plt.errorbar(x=[row['OR']], y=[row.name], xerr=ci,
                    ecolor='tab:gray', capsize=3, linestyle='None', linewidth=1, marker="o", 
                            markersize=5, mfc="tab:blue", mec="tab:blue")
            else:
                plt.errorbar(x=[row['OR']], y=[row.name], xerr=ci,
                    ecolor='tab:gray', capsize=3, linestyle='None', linewidth=1, marker="o", 
                            markersize=5, mfc="tab:gray", mec="tab:gray")
        elif row['significant?'] == 'not significant': ## if it is not significant we color gray
            plt.errorbar(x=[row['OR']], y=[row.name], xerr=ci,
                ecolor='tab:gray', capsize=3, linestyle='None', linewidth=1, marker="o", 
                        markersize=5, mfc="tab:gray", mec="tab:gray")

         # Add odds ratio value as text annotation AND BOLD
        plt.text(row['OR'], row.name, f"{row['OR']:.2f}", ha='center', va='bottom', fontsize=9, fontweight='bold')


    plt.axvline(x=1, linewidth=1.5, linestyle='-', color='black')
    plt.axhline(y=10.7, linewidth=1.2, linestyle='--', color='grey')  
    plt.axhline(y=18.7, linewidth=1.2, linestyle='--', color='grey')  
    plt.axhline(y=21.7, linewidth=1.2, linestyle='--', color='grey')  
    plt.axhline(y=23.7, linewidth=1.2, linestyle='--', color='grey')  
    plt.axhline(y=28.7, linewidth=1.2, linestyle='--', color='grey')  
    plt.axhline(y=31.7, linewidth=1.2, linestyle='--', color='grey')  

    #plt.axhline(y=12, linewidth=1.5, linestyle='--', color='black')  


    ## ADD ODS RATIO VALUE TO DOTS

    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlabel('Odds Ratio and 95% Confidence Interval', fontsize=12,fontweight='bold')
    plt.tight_layout()
    ## SET MAXIMUM AND MINIMUM VALUES FOR X-AXIS
    plt.xlim(-1, 20)

    # add_label_band(ax, 1.5, 7.4, "Physicochem Prop AA Change ")
    # add_label_band(ax, 7.65, 15.4, "Physicochem Prop REF AA  ")
    # add_label_band(ax, 15.65, 18.4, "Residue Exposure   ")
    # add_label_band(ax, 18.65, 22.4, "Secondary Structure  ")
    ## ADD XVALUE TO PLOT
    ## remove vertical gridlines
 
    plt.show()

    return 




import matplotlib.pyplot
from matplotlib.patches import Rectangle

def plot_heatmap(test_features_groupby, list_with_formated_names, list_order):
    """
    Function to obtain a heatmap 
     Parameters
    ------------------ 
    test_features_groupby: pandas.dataframe
        df that contains the groupby, oddsratio(OR), p-value, q-value, lowerconfidence interval(CIL),
        upper confidence interval(CIU), significant and ORs_clasification for all our features of interest.
        and as index has the features named as 'feat_assignment'

    list_with_formated_names: list
        list that contains the names of the features in a nice format to use it in AXIS, ojo need to be in same order as received
    
    list_order: list
        list that contains the order of the x axis, Ex protein classes in order so they can be plotted class1, class2, class3 ...
    """

    # First we create a pivot table with features -> `index`, `groupby``, and `OR`
    test_features_groupby.reset_index(inplace=True) # features are in the index, so we need to make it a column to use it for pivot. Reset index -> we get a column named 'index' with all the features
    test_features_groupby.loc[test_features_groupby['significant?'] == 'not significant', 'OR'] = np.nan ## REPLACE NOT SIGNIFICANT VALUES WITH NAN 

    glue = test_features_groupby.pivot(index="index", columns="groupby", values="OR")

    # Reorder according to the desired order
    glue = glue.reindex(list_with_formated_names)

    # Set 'features_rename' as the new index
    glue.insert(0, 'features_rename', list_with_formated_names)
    glue = glue.set_index('features_rename')

    # Reorder columns
    glue = glue[list_order]

    # Replace infinite values with 1
    #glue = glue.replace([np.inf, -np.inf], 1)

    # Plot the heatmap
    fig, ax = plt.subplots(figsize=(23, 32))  # Sample figsize in inches
    plot = sns.heatmap(glue, cmap="coolwarm", vmin=0.05, vmax=5, center=1, xticklabels=True,linewidth=.5,annot=True,annot_kws={"fontsize": 18})
    #plot = sns.heatmap(glue, cmap="coolwarm", vmin=0.05, vmax=5, center=1, xticklabels=True,linewidth=.5)

    ## set x axis size to 22


    #for i in range(len(glue.index)): ## for each feature
    for i in range(len(glue.index)): ## for each feature
        for j in range(len(glue.columns)): ## for each groupby
            if np.isnan(glue.iloc[i,j]): ## if value is nan
                plot.add_patch(Rectangle((j, i), 1, 1, fill=True, color='#ececf4'))

    #plot.set_facecolor('xkcd:silver')
    plot.tick_params(labelsize=36)
    plot.set_xticklabels(plot.get_xticklabels(), rotation=90, horizontalalignment='right')
    plot.set(xlabel=None, ylabel=None)

    sns.set(font_scale=1.4)
    return plot









