from plots_and_functions import plot_distribution
import pandas as pd
import os
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as mtick

from matplotlib.ticker import MaxNLocator, PercentFormatter
# import numpy as np
from score_names import score_names
import seaborn as sns
import numpy as np
import math
import re
import upsetplot as up
from supervenn import supervenn


def combined_plot_old(df, ylim_top, out_name):
    """
    Plot a hist and a bar plot for each search engine mentioned in the dataframe

    Args:
        df (dataframe): dataframe with all combined results
        ylim_top (_type_): search specific ylimits for the histogram plot
        out_name (_type_): where to save the plot (without extension)

    Returns:
        dataframe: summary of results
    """
    # combined plot: all scores in one plot
    df.sort_values(by='search_engine', inplace=True, key=lambda x: x.str.lower())
    fig, ax = plt.subplots(math.ceil(df['search_engine'].nunique()/2), 4, figsize=(12, 12),
                           gridspec_kw={'width_ratios': [1, 0.1, 1, 0.1]})
    ax = ax.flatten()
    i = 0
    summary = []
    for search_engine, group in df.groupby('search_engine', sort=False):

        xcol = score_names[search_engine] if score_names[search_engine] in group.columns else "Score"
        _, hp, colors = plot_distribution(df=group, x=xcol, bins=50, ax=ax[i], ylim_top=ylim_top.get(search_engine, None))
        ax[i].set_xlabel(f'{search_engine} ({score_names[search_engine]})')
        if i % 2 == 1:
            ax[i].set_ylabel('')
        ax[i].xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax[i].xaxis.set_minor_locator(MaxNLocator(nbins=10))
        ax[i].yaxis.set_major_locator(MaxNLocator(nbins=2))
        ax[i].spines[['right', 'top']].set_visible(False)
        if search_engine == 'MeroX':
            ax[i].axvline(50, color='black')  # in default settings, this is the minimal score cutoff

        # legend
        if i == 0:
            handles, labels = ax[i].get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper center')
        else:
            ax[i].legend().remove()
        i += 1

        total = group.shape[0]
        group = group.reset_index()
        EE = sum(group['EE'])
        EH = sum(group['EH'])
        HH = sum(group['HH'])
        TT_EE = sum(group['EE'][group['isTT'] == True])
        TT_EH = sum(group['EH'][group['isTT'] == True])
        TT_HH = sum(group['HH'][group['isTT'] == True])
        TD = sum(group['isTD'])
        TT = sum(group['isTT'])
        DD = sum(group['isDD'])

        if DD > TD:
            FDR = (TD+DD)/TT
            note = "*"
        else:
            FDR = (TD-DD)/TT
            note = ""
        ENTRAP_TOTAL = TT_EH + TT_HH
        ENTRAP_TTFP = ENTRAP_TOTAL + (TT_EH - TT_HH)
        ENTRAP_FPR = (EH - HH)/EE

        summary.append({
            "search_engine": search_engine,
            "TT": TT,
            "TD": TD,
            "DD": DD,
            "EE(TT)": TT_EE,
            "EH(TT)": TT_EH,
            "HH(TT)": TT_HH,
            "EE": EE,
            "EH": EH,
            "HH": HH,
            "FDR(Decoy)": FDR,
            "FDR(Entrapment)": ENTRAP_FPR,
            "TT_FPR(Entrapment)": ENTRAP_TTFP/TT
            })

        decoy_based_error = FDR * 100
        entr_based_error = ENTRAP_FPR * 100

        sns.barplot(x=['D'+note, 'E'], y=[decoy_based_error, entr_based_error], ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#E364B4']))
        #sns.barplot(x=['D'+note, 'E'], y=[decoy_based_error, ENTRAP_TOTAL / TT * 100], ax=ax[i],
        #             palette=sns.color_palette(['#83320C', '#C354A4']))
        ax[i].spines[['right', 'top']].set_visible(False)
        if (entr_based_error < 8) & (decoy_based_error < 8):
            ax[i].set_ylim(0, 8)
        else:
            ax[i].set_ylim(0, math.ceil(entr_based_error/2)*2)
        ax[i].axhline(2, color='black', linestyle='--')
        ax[i].yaxis.set_major_locator(MaxNLocator(nbins=4))
        # ax[i].yaxis.set_minor_locator(MaxNLocator(nbins=8))
        ax[i].set_ylabel('error [%]')
        i += 1

    plt.tight_layout()
    # plt.show()

    plt.savefig(out_name + '.svg')
    plt.savefig(out_name + '.png')
    plt.close('all')
    return pd.DataFrame(summary)


def bar_plots(df_csm, df_peppair, df_respair, df_ppi, out_name, possible_ppis=None):
    """
    Plot a hist and a bar plot for each search engine mentioned in the dataframe

    Args:
        df (dataframe): dataframe with all combined results
        ylim_top (_type_): search specific ylimits for the histogram plot
        out_name (_type_): where to save the plot (without extension)

    Returns:
        dataframe: summary of results
    """
    # mark each dataframe with the level
    for df, level, level_number in zip([df_csm, df_peppair, df_respair, df_ppi], ['CSM', 'Peptide Pair', 'Residue Pair', 'PPI'], range(4)): 
        if df is not None:
            df['level'] = level
            df['level_number'] = level_number

    # combine all dataframes
    df = pd.concat([x for x in [df_csm, df_peppair, df_respair, df_ppi] if x is not None], ignore_index=True)
    df['search_engine_lc'] = df['search_engine'].str.lower()
    df.sort_values(by=['search_engine_lc', "level_number"], inplace=True)

    all_possible_pairs_dict=possible_ppi_to_dict(possible_ppis)

    df['possible'] = df[['Protein1', 'Protein2']].apply(
        lambda x: 1 if any(p1 in all_possible_pairs_dict and 
                           any (p2 in all_possible_pairs_dict[p1]  for p2 in x.iloc[1].split(";"))  
                           for p1 in x.iloc[0].split(";")) else 0, axis=1)

    df['both_found'] = df[['Protein1', 'Protein2']].apply(
        lambda x: 1 if any([p1 in all_possible_pairs_dict.keys() for p1 in x.iloc[0].split(";")]) and 
        any([p2 in all_possible_pairs_dict.keys() for p2 in x.iloc[1].split(";")]) else 0, axis=1)
                           




    # combined plot: all scores in one plot
    fig, ax = plt.subplots(math.ceil(df['search_engine'].nunique()/2), 2, figsize=(12, 12))
    ax = ax.flatten()
    i = 0
    summary = []
    numbers = []
    for search_engine, group in df.groupby('search_engine', sort=False):
        if i == 0:
            handles, labels = ax[i].get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper center')
        else:
            ax[i].legend().remove()

        # summarize results for each level
        results_e_observed = []
        results_e_estimate = []
        results_te_estimate = []
        for level, level_group in group.groupby('level', sort=False):
            
            level_group = level_group.reset_index()
            level_group_possible = level_group[(level_group['isTT'] == True) & (level_group['EE'] == True)]
            #level_group_possible = pd.merge(level_group_possible.reset_index(), all_possible_pairs, on=['Protein1', 'Protein2'], how='left')

            EE = sum(level_group['EE'])
            EH = sum(level_group['EH'])
            HH = sum(level_group['HH'])
            isTT_mask = level_group['isTT'] == True
            isEE_mask = level_group['EE'] == True
            is_found_mask = level_group['both_found'] == 1
            is_possible_mask = level_group['possible'] == 1
            TT_EE = sum(level_group['EE'][isTT_mask])
            TT_EE_possible = sum(level_group['EE'][isTT_mask & is_possible_mask])
            TT_EE_found = sum(isEE_mask & isTT_mask & is_found_mask)
            sum((level_group_possible['both_found'] == 1) & (level_group_possible['isTT'] == True))
            TT_EE_impossible = TT_EE - TT_EE_possible
            TT_EH = sum(level_group['EH'][level_group['isTT'] == True])
            TT_HH = sum(level_group['HH'][level_group['isTT'] == True])
            DD = sum(level_group['isDD'])
            TD = sum(level_group['isTD'])
            TT = sum(level_group['isTT'])
            #TT_EE = TT-TT-EH-TT_HH

            level_number = level_group['level_number'].iloc[0]

            if DD > TD:
                FDR = (TD+DD)/TT
                fdr_name = "%decoys"
            else:
                FDR = (TD-DD)/TT
                fdr_name = "FDR"
            ENTRAP_TOTAL = TT_EH + TT_HH

            ENTRAP_FP = ENTRAP_TOTAL + (TT_EH - TT_HH)
            results_e_observed.extend([{
                "level": level,
                "order": level_number * 10 + 2,
                "what": "Observed entrapment",
                "value": min(ENTRAP_TOTAL / TT * 100, 100)
            },
            {
                "level": level,
                "order": level_number * 10 + 1,
                "what": fdr_name,
                "value": min(FDR * 100, 100)
            },
            {
                "level": level,
                "order": level_number * 10 + 3,
                "what": "Impossible PPI among E-Coli TT ",
                "value": (1- TT_EE_possible / TT_EE) *100
            },
            {
                "level": level,
                "order": level_number * 10 + 4,
                "what": "Total observable False Positives among TT",
                "value": ((TT_EE_impossible + ENTRAP_TOTAL)/ TT) *100
            }])        

            results_e_estimate.extend([{
                "level": level,
                "order": level_number * 10 + 2,
                "what": "Entrapment based error",
                "value": min(100, (EH - HH) / EE * 100)
            },
            {
                "level": level,
                "order": level_number * 10 + 1,
                "what": "_dummy1",
                "value": 0
            },
            {
                "level": level,
                "order": level_number * 10 + 3,
                "what": "_dummy2",
                "value": 0
            },
            {
                "level": level,
                "order": level_number * 10 + 4,
                "what": "_dummy3",
                "value": 0
            }])

            results_te_estimate.extend([{
                "level": level,
                "order": level_number * 10 + 2,
                "what": "Entrapment based error among TT",
                #  total estimate FPR = (observed Human + estimated FP)/total TT
                # total estimate FPR = (EH+HH + estimates FP) / TT
                # total estimate FP = (EH+HH + (EH-HH))/TT
                # total estimate FP = 2*EH/TT
                "value": min(100, (2*TT_EH)/TT*100)
            },
            {
                "level": level,
                "order": level_number*10 + 1,
                "what": "_dummy1",
                "value": 0
            },
            {
                "level": level,
                "order": level_number*10 + 3,
                "what": "_dummy2",
                "value": 0
            },
            {
                "level": level,
                "order": level_number * 10 + 4,
                "what": "_dummy3",
                "value": 0
            }])

            numbers.extend([{
                "search_engine": search_engine,
                "level": level,
                "TT": TT,
                "TD": TD,
                "DD": DD,
                "EE": EE,
                "EH": EH,
                "HH": HH,
                "TT_EE": TT_EE,
                "TT_EE_possible": TT_EE_possible,
                "TT_EE_found": TT_EE_found,
                "TT_EE_impossible": TT_EE_impossible,
                "TT_EH": TT_EH,
                "TT_HH": TT_HH,
                "TD_EE": sum(level_group['EE'][level_group['isTD']]),
                "TD_EH": sum(level_group['EH'][level_group['isTD']]),
                "TD_HH": sum(level_group['HH'][level_group['isTD']]),
                "DD_EE": sum(level_group['EE'][level_group['isDD']]),
                "DD_EH": sum(level_group['EH'][level_group['isDD']]),
                "DD_HH": sum(level_group['HH'][level_group['isDD']]),
                "Total known FP in TT" : (TT_EE_impossible + ENTRAP_TOTAL),
                "FDR": FDR,
                "%Known FP": (TT_EE_impossible + ENTRAP_TOTAL) / TT,
                "True Positives<": (TT_EE - TT_EE_impossible),
                "True Positives%<": (TT_EE - TT_EE_impossible) / TT,
            },
            ])

        # convert to dataframe
        df_results_e_observed = pd.DataFrame(results_e_observed)
        df_results_e_estimate = pd.DataFrame(results_e_estimate)
        df_results_te_estimate = pd.DataFrame(results_te_estimate)

        df_results_e_observed.sort_values(by='order', inplace=True)
        df_results_e_estimate.sort_values(by='order', inplace=True)
        df_results_te_estimate.sort_values(by='order', inplace=True)
                
        sns.barplot(data=df_results_e_estimate, x='level', y='value', hue="what", ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#EB92CA', '#0092CA', '#0052CA']))
        sns.barplot(data=df_results_e_observed, x='level', y='value', hue="what", ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#E364B4', '#0092CA', '#0052CA']))
        ax[i].spines[['right', 'top']].set_visible(False)
        max_error = max(max(df_results_e_estimate['value']), max(df_results_e_observed['value']))
        if max_error < 10:
            ax[i].set_ylim(0, 10)
        else:
            ax[i].set_ylim(0, math.ceil(max_error/2)*2)
            
        #if entr_based_error < 8:
        #    ax[i].set_ylim(0, 8)
        #else:
        #    ax[i].set_ylim(0, round(entr_based_error))
        ax[i].axhline(2, color='black', linestyle='--')
        ax[i].set_title(search_engine)
        #ax[i].yaxis.set_major_locator(MaxNLocator(nbins=4))
        #ax[i].yaxis.set_minor_locator(MaxNLocator(nbins=8))
        ax[i].set_ylabel('error [%]')
        #if i != 0:
        ax[i].legend().remove()
        i += 1

    legends = [ mpatches.Patch(color='#83320C', label='Decoy based error*'),
                mpatches.Patch(color='#E364B4', label='Observed entrapment'),
                mpatches.Patch(color='#EB92CA', label='Entrapment based error'),
                mpatches.Patch(color='#0092CA', label='Impossible PPI among E-Coli TT'),
                mpatches.Patch(color='#0052CA', label='Total observable False Positives among TT')]
    fig.legend(handles=legends, loc='lower center')
    
    #ax[i].legend([])
    df_numbers = pd.DataFrame(numbers)
    df_numbers.to_csv(out_name + '_numbers.csv', index=False)
    plt.tight_layout()
    # plt.show()

    plt.savefig(out_name + '.svg')
    plt.savefig(out_name + '.png')
    plt.close('all')
    return pd.DataFrame(summary)


def getPossibleProteinPairs(pgfile, ibaq_cutoff=0):
    """
    Based on an maxquant proteingroupfile get all possible protein pairs that have in any sample together with both protein ibaqs > ibaq_cutoff

    Args:
        pgfile (str): path to proteingroup file
        ibaq_cutoff (float, optional): ibaq cutoff. Defaults to 0.

    Returns:
        dataframe: possible protein pairs
    """
    df = pd.read_csv(pgfile, sep='\t')
    df = df[df['Reverse'] != '+']
    df = df[df['Potential contaminant'] != '+']
    df = df[df['Only identified by site'] != '+']
    
    # split protein groups into single proteins
    df['Protein IDs'] = df['Protein IDs'].apply(lambda x: x.split(';'))
    df = df.explode('Protein IDs')

    # iterate over all iBAQ columns and check if iBAQ > ibaq_cutoff
    # if yes, add protein to set of proteins for this sample
    # if no, add empty set
    sample_proteins = {}
    for col in df.columns:
        if col.startswith('iBAQ '):
            sample_proteins[col] = set(df['Protein IDs'][df[col] > ibaq_cutoff])

    # iterate over all samples and add all possible protein pairs to a set
    possible_pairs = set()
    for sample, proteins in sample_proteins.items():
        for p1 in proteins:
            for p2 in proteins:
                    possible_pairs.add((p1, p2))

    # convert set to dataframe
    df_possible_pairs = pd.DataFrame(list(possible_pairs), columns=['Protein1', 'Protein2'])
    return df_possible_pairs


def combined_plot(df, ylim_top, out_name, possible_ppis=None, barplots=True):
    """
    Plot a hist and a bar plot for each search engine mentioned in the dataframe

    Args:
        df (dataframe): dataframe with all combined results
        ylim_top (_type_): search specific ylimits for the histogram plot
        out_name (_type_): where to save the plot (without extension)

    Returns:
        dataframe: summary of results
    """
    # combined plot: all scores in one plot
    df.sort_values(by='search_engine', inplace=True, key=lambda x: x.str.lower())
    if barplots:
        fig, ax = plt.subplots(math.ceil(df['search_engine'].nunique()/2), 4, figsize=(12, 12),
                               gridspec_kw={'width_ratios': [1, 0.1, 1, 0.1]})
    else:
        fig, ax = plt.subplots(math.ceil(df['search_engine'].nunique()/2), 2, figsize=(12, 12))
        #fig, ax = plt.subplots(math.ceil(df['search_engine'].nunique()/3), 3, figsize=(12, 12))        
    ax = ax.flatten()
    i = 0
    if possible_ppis is not None:
        possible_ppi_dict = possible_ppi_to_dict(possible_ppis)
        df['impossibe'] = df[['Protein1', 'Protein2']].apply(
            lambda x: 1 if any(p1 in possible_ppi_dict and any(p2 in possible_ppi_dict[p1] for p2 in x.iloc[1].split(";")) for p1 in x.iloc[0].split(";")) else 0, axis=1)  

    

    handles = []
    labels = []
    summary = []

    for search_engine, group in df.groupby('search_engine', sort=False):

        xcol = score_names[search_engine] if score_names[search_engine] in group.columns else "Score"
        _, hp, colors = plot_distribution(df=group, x=xcol, bins=50, ax=ax[i], ylim_top=ylim_top.get(search_engine, None))
        #ax[i].set_xlabel(f'{search_engine} ({score_names[search_engine]})')
        ax[i].set_xlabel(f'{score_names[search_engine]}')
        if i % 2 == 1:
            ax[i].set_ylabel('')
        ax[i].xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax[i].xaxis.set_minor_locator(MaxNLocator(nbins=10))
        ax[i].yaxis.set_major_locator(MaxNLocator(nbins=2))
        ax[i].spines[['right', 'top']].set_visible(False)
        if search_engine == 'MeroX':
            ax[i].axvline(50, color='black')  # in default settings, this is the minimal score cutoff
        ax[i].text((ax[i].get_xlim()[1]-ax[i].get_xlim()[0])/2+ax[i].get_xlim()[0],
                   (ax[i].get_ylim()[1]-ax[i].get_ylim()[0])/2+ax[i].get_ylim()[0],
                   f'{search_engine}', horizontalalignment='center', verticalalignment='center')
        # legend
        # if i == 0:
        #     handles0 = [x for x in ax[0].get_legend().legend_handles]
        #     labels0 = ax[0].get_legend().get_texts()
        #     handles.extend(handles0)
        #     labels.extend(labels0)

        ax[i].legend().remove()
        i += 1

        total = group.shape[0]
        group = group.reset_index()
        EE = sum(group['EE'])
        EH = sum(group['EH'])
        HH = sum(group['HH'])
        TT_EE = sum(group['EE'][group['isTT'] == True])
        TT_EH = sum(group['EH'][group['isTT'] == True])
        TT_HH = sum(group['HH'][group['isTT'] == True])
        TD = sum(group['isTD'])
        TT = sum(group['isTT'])
        DD = sum(group['isDD'])

        if DD > TD:
            FDR = (TD+DD)/TT
            note = "*"
        else:
            FDR = (TD-DD)/TT
            note = ""
        ENTRAP_TOTAL = TT_EH + TT_HH
        ENTRAP_TTFP = ENTRAP_TOTAL + (TT_EH - TT_HH)
        ENTRAP_FPR = (EH - HH)/EE
        
        se_summary = {
            "search_engine": search_engine,
            "TT": TT,
            "TD": TD,
            "DD": DD,
            "EE(TT)": TT_EE,
            "EH(TT)": TT_EH,
            "HH(TT)": TT_HH,
            "EE": EE,
            "EH": EH,
            "HH": HH,
            "FDR(Decoy)": FDR,
            "FDR(Entrapment)": ENTRAP_FPR,
            "TT_FPR(Entrapment)": ENTRAP_TTFP/TT
            }
        if possible_ppis is not None:
            TT_possible = sum(group['impossibe'][group['isTT'] == True])
            se_summary["TT_possible"] = TT_possible
            se_summary["TT_impossible"] = TT - TT_possible
            se_summary["TT_EE_impossible"] = TT_EE - TT_possible
            TT_percent_impossible = (1 - TT_possible/TT) * 100
            se_summary["TT_impossible%"] = TT_percent_impossible

        summary.append(se_summary)

        decoy_based_error = FDR * 100
        entr_based_error = ENTRAP_FPR * 100
        barx = ['D'+note, 'E']
        bary = [decoy_based_error, entr_based_error]
        barcol = ['#83320C', '#E364B4']
        if possible_ppis is not None:
            barx.append('I')
            bary.append(TT_percent_impossible)
            barcol.append('#C354A4')

        if barplots:
            sns.barplot(x=barx, y=bary, ax=ax[i],
                        palette=sns.color_palette(barcol))
            #sns.barplot(x=['D'+note, 'E'], y=[decoy_based_error, ENTRAP_TOTAL / TT * 100], ax=ax[i],
            #             palette=sns.color_palette(['#83320C', '#C354A4']))
            ax[i].spines[['right', 'top']].set_visible(False)
            if entr_based_error < 8:
                ax[i].set_ylim(0, 8)
            else:
                ax[i].set_ylim(0, math.ceil(entr_based_error/2)*2)
            ax[i].axhline(2, color='black', linestyle='--')
            ax[i].yaxis.set_major_locator(MaxNLocator(nbins=4))
            # ax[i].yaxis.set_minor_locator(MaxNLocator(nbins=8))
            ax[i].set_ylabel('error [%]')

            # if i == 1:
            #     ax[1].legend()
            #     handles0 = [mpatches.Patch(color='#83320C', label='Decoy based Error'),
            #                 mpatches.Patch(color='#E364B4', label='Entrapment based Error')]
            #     handles.extend(handles0)
            ax[i].legend().remove()


            i += 1

    handles = [mpatches.Patch(edgecolor='#8CCFC0', facecolor='#8CCFC070', label='Target distribution'),
                mpatches.Patch(edgecolor='#83320C', facecolor='#83320C70', label='Decoy distribution'),
                mpatches.Patch(edgecolor='#E364B4', facecolor='#E364B470', label='Entrapment distribution'),
                mpatches.Patch(color='#83320C', label='Decoy based Error'),
                mpatches.Patch(color='#E364B4', label='Entrapment based Error')]
    labels = ['Target distribution', 'Decoy distribution', 'Entrapment distribution', 'Decoy based Error', 'Entrapment based Error']

    #fig.legend(handles , labels, loc='outside lower center')
    fig.legend(handles, labels, loc='outside lower center', ncols=3, bbox_to_anchor=(0.5, 0.01), bbox_transform=plt.gcf().transFigure)

    plt.tight_layout()
    # plt.show()


    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(out_name + '.svg')
    plt.savefig(out_name + '.png')
    plt.close('all')
    ret = pd.DataFrame(summary)
    ret.to_csv(out_name + '.csv')
    return ret


def possible_ppi_to_dict(possible_ppis):
    # make sure we have all combinations
    possible_ppis_rev = possible_ppis[["Protein2", "Protein1"]]
    possible_ppis_rev.columns = ["Protein1", "Protein2"]
    all_possible_pairs = pd.concat([possible_ppis, possible_ppis_rev])

    all_possible_pairs['possible'] = 1
    
    all_possible_pairs_dict = {}


    # convert possible pairs to dict of sets
    for p1, p2 in all_possible_pairs[["Protein1", "Protein2"]].values:
        if p1 in all_possible_pairs_dict:
            p1_pairs = all_possible_pairs_dict[p1] 
        else:
            p1_pairs = set()
            all_possible_pairs_dict[p1] = p1_pairs
        p1_pairs.add(p2)
        
        if p2 in all_possible_pairs_dict:
            p2_pairs = all_possible_pairs_dict[p2] 
        else:
            p2_pairs = set()
            all_possible_pairs_dict[p2] = p2_pairs
        p2_pairs.add(p1)
    return all_possible_pairs_dict

def overlap_plots(df, out_name, possible_ppis=None, EEonly=False, plot_impossible=False, plot_TT=True, plot_TD=False, plot_DD=False):
    """
    Plot the overlap of peptide pairs between search engines

    Args:
        df (dataframe): dataframe containing all peptide pairs
        out_name (str): where to save the plot (without extension)
        possible_ppis (dataframe, optional): dataframe containing all possible protein pairs. Defaults to None.
        EEonly (bool, optional): only consider EE peptide pairs. Defaults to False.
        plot_impossible (bool, optional): plot the overlap of impossible peptide pairs. Defaults to False.
    """
    # remove plink cleav
    df_peppair = df[df['search_engine'] != 'pLink2_cleav']
    if not plot_TT:
        df_peppair = df_peppair[~df_peppair['isTT']]
    if not plot_TD:
        df_peppair = df_peppair[~df_peppair['isTD']]
    if not plot_DD:
        df_peppair = df_peppair[~df_peppair['isDD']]
    # turn peptide sequences into basic - non modified sequences
    if "PepSeq1" in df_peppair.columns:
        df_peppair['Base1'] = df_peppair['PepSeq1'].apply(lambda x: re.sub(r'[^A-Z]', '', x))
        df_peppair['Base2'] = df_peppair['PepSeq2'].apply(lambda x: re.sub(r'[^A-Z]', '', x))
        # for n-terminal peptides starting with M remove it
        df_peppair.loc[(df_peppair['PepPos1'] == 1) & (df_peppair['Base1'].str.startswith('M')), 'Base1'] = df_peppair['Base1'].str[1:]
        df_peppair.loc[(df_peppair['PepPos2'] == 1) & (df_peppair['Base2'].str.startswith('M')), 'Base2'] = df_peppair['Base2'].str[1:]

        
    elif "Peptide1" in df_peppair.columns:
        df_peppair['Base1'] = df_peppair['Peptide1'].apply(lambda x: re.sub(r'[^A-Z]', '', x))
        df_peppair['Base2'] = df_peppair['Peptide2'].apply(lambda x: re.sub(r'[^A-Z]', '', x))

        df_peppair.loc[(df_peppair['Start1'] == 1) & (df_peppair['Base1'].str.startswith('M')), 'Base1'] = \
            df_peppair.loc[(df_peppair['Start1'] == 1) & (df_peppair['Base1'].str.startswith('M')), 'Base1'].str[1:]
        df_peppair.loc[(df_peppair['Start2'] == 1) & (df_peppair['Base2'].str.startswith('M')), 'Base2'] = \
            df_peppair.loc[(df_peppair['Start2'] == 1) & (df_peppair['Base2'].str.startswith('M')), 'Base2'].str[1:]
        
    elif "fromSite1" in df_peppair.columns: # make the pair based on protein sites
        df_peppair['Base1'] = df_peppair[['fromSite1', 'Protein1']].apply(lambda x: str(x.iloc[1]) + '_ ' + x.iloc[0], axis=1)
        df_peppair['Base2'] = df_peppair[['fromSite2', 'Protein2']].apply(lambda x: str(x.iloc[1]) + '_ ' + x.iloc[0], axis=1)
    else: # just protein pairs
        df_peppair['Base1'] = df_peppair['Protein1']
        df_peppair['Base2'] = df_peppair['Protein2']



    # make two columns for alphabeticelly first and last peptide pair
    if 'run' in df_peppair.columns and 'scan' in df_peppair.columns:
        df_peppair['pair'] = df_peppair[['Base1', 'Base2', 'run', 'scan']].apply(lambda x: ';'.join(sorted([x.iloc[0], x.iloc[1]]))+';'+ x[2]+';' + str(int(x[3])), axis=1)
    else:
        df_peppair['pair'] = df_peppair[['Base1', 'Base2']].apply(lambda x: ';'.join(sorted([x.iloc[0], x.iloc[1]])), axis=1)

    if possible_ppis is not None:
        possible_ppis_rev = possible_ppis[["Protein2" , "Protein1"]]
        possible_ppis_rev.columns = ["Protein1" , "Protein2"]
        all_possible_pairs=pd.concat([possible_ppis, possible_ppis_rev])

        all_possible_pairs['possible'] = 1
        
        all_possible_pairs_dict={}


        # convert possible pairs to dict of sets
        for p1 , p2 in all_possible_pairs[["Protein1" , "Protein2"]].values:
            p1_pairs = all_possible_pairs_dict[p1] if p1 in all_possible_pairs_dict else set()
            p1_pairs.add(p2)
            all_possible_pairs_dict[p1] = p1_pairs
            p2_pairs = all_possible_pairs_dict[p2] if p2 in all_possible_pairs_dict else set()
            p2_pairs.add(p1)
            all_possible_pairs_dict[p2] = p2_pairs

        df_peppair['possible'] = df_peppair[['Protein1', 'Protein2']].apply(
            lambda x: 1 if any(p1 in all_possible_pairs_dict and 
                            any (p2 in all_possible_pairs_dict[p1]  for p2 in x.iloc[1].split(";"))  
                            for p1 in x.iloc[0].split(";")) else 0, axis=1)
        df_peppair_possible = df_peppair[df_peppair['possible'] == 1]
        df_peppair_impossible = df_peppair[df_peppair['possible'] == 0]
        
        if plot_impossible:
            df_peppair = df_peppair_impossible
        else:
            df_peppair = df_peppair_possible

    # filter for unique peptide pairs - disregrading linkage sites
    df_peppair_unique = df_peppair #[df_peppair['isTT'] == True]
    #df_peppair_imposs_unique = df_peppair_impossible[df_peppair_impossible['isTT'] == True]

    if EEonly:
        df_peppair_unique = df_peppair_unique[df_peppair['EE']]

    # filter to unique peptide pairs in each search engine
    df_peppair_unique = df_peppair_unique.drop_duplicates(subset=['pair', 'search_engine'])
    # count how many unique peptide pairs are found by each search engine
    df_peppair_total = df_peppair_unique.groupby('search_engine').size().reset_index(name='counts')
    # rename counts to toal
    df_peppair_total.rename(columns={'counts': 'total'}, inplace=True)

    # count for each peptide-pair in how many search engines it is found
    df_peppair_unique['search_counts'] = df_peppair_unique.groupby('pair')['search_engine'].transform('count')
    # count how many peptide pairs are found by how many search engines for each search engine
    df_peppair_overlap = df_peppair_unique.groupby('search_engine')['search_counts'].value_counts().reset_index(name='counts')
    # rename columns
    df_peppair_overlap.rename(columns={'search_counts': 'overlap'}, inplace=True)
    # substract 1 from overlap
    df_peppair_overlap['overlap'] = df_peppair_overlap['overlap'] - 1
    # sum up the total number of peptide pairs found by each search engine
    df_peppair_overlap = pd.merge(df_peppair_overlap, df_peppair_total, on='search_engine', how='left')
    # calculate percentage of peptide pairs found by how many search engines
    df_peppair_overlap['percentage'] = df_peppair_overlap['counts'] / df_peppair_overlap['total'] * 100
    # add lower case search engine name
    df_peppair_overlap['search_engine_lc'] = df_peppair_overlap['search_engine'].str.lower()
    # sort by search engine
    df_peppair_overlap.sort_values(by=['search_engine_lc'], inplace=True)

    # plot a stacked bar plot of the overlap in absolute numbers
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    sns.barplot(data=df_peppair_overlap, x='search_engine', y='counts', hue='overlap', ax=ax[0], palette='magma')
    ax[0].set_ylabel('Number of unique peptide pairs')
    ax[0].set_xlabel('')
    ax[0].set_title('Number of unique peptide pairs found by how many search engines')
    ax[0].spines[['right', 'top']].set_visible(False)
    #ax[0].legend().remove()
    # turn the search engine names by 45%
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=45, horizontalalignment='right')
    # plot a stacked bar plot of the overlap in percentages
    sns.barplot(data=df_peppair_overlap, x='search_engine', y='percentage', hue='overlap', ax=ax[1], palette='magma')
    ax[1].set_ylabel('Percentage of unique peptide pairs')
    ax[1].set_xlabel('')
    ax[1].set_title('Percentage of unique peptide pairs found by how many search engines')
    ax[1].spines[['right', 'top']].set_visible(False)
    #ax[1].legend().remove()
    ax[1].set_xticklabels(ax[0].get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    plt.savefig(out_name + '_overlap.svg')
    plt.savefig(out_name + '_overlap.png')
    plt.close('all')

    pairwise_overlap = []
    # loop over each search engine and count the overlap to every other search engine 
    for search_engine1, group1 in df_peppair_unique.groupby('search_engine'):
        for search_engine2, group2 in df_peppair_unique.groupby('search_engine'):
            # find the overlap of Peptide pairs between search_engine1 and search_engine2
            df_peppair_overlap_se = pd.merge(group1, group2, on='pair', how='inner')
            pairwise_overlap.append({
                'search_engine1': search_engine1,
                'search_engine2': search_engine2,
                'overlap': df_peppair_overlap_se.shape[0]
            })

    # convert to dataframe
    df_pairwise_overlap = pd.DataFrame(pairwise_overlap)
    # sort by search engine lower case name
    df_pairwise_overlap['search_engine1_lc'] = df_pairwise_overlap['search_engine1'].str.lower()
    df_pairwise_overlap['search_engine2_lc'] = df_pairwise_overlap['search_engine2'].str.lower()
    df_pairwise_overlap.sort_values(by=['search_engine1_lc', 'search_engine2_lc'], inplace=True)
    # drop lower case columns
    df_pairwise_overlap.drop(columns=['search_engine1_lc', 'search_engine2_lc'], inplace=True)


    # pivot the dataframe to a matrix
    df_pairwise_overlap_pivot = df_pairwise_overlap.pivot(index='search_engine1', columns='search_engine2', values='overlap')
    # fill up the diagonal with search_engine uniqe counts
    for se in df_pairwise_overlap_pivot.index:
        df_pairwise_overlap_pivot[se][df_pairwise_overlap_pivot.index == se] = df_peppair_overlap['counts'][(df_peppair_overlap['search_engine'] == se) & (df_peppair_overlap['overlap'] == 0)].iloc[0]

    ## fill the diagonal with 0
    #np.fill_diagonal(df_pairwise_overlap_pivot.values, np.nan)
    # sort by lower case search_engine name
    df_pairwise_overlap_pivot['search_engine'] = df_pairwise_overlap_pivot.index
    df_pairwise_overlap_pivot['search_engine_lc'] = df_pairwise_overlap_pivot['search_engine'].str.lower()
    df_pairwise_overlap_pivot.sort_values(by='search_engine_lc', inplace=True)
    df_pairwise_overlap_pivot = df_pairwise_overlap_pivot[df_pairwise_overlap_pivot['search_engine'].tolist()]

    ## plot as heatmap
    #fig, ax = plt.subplots(figsize=(12, 12))
    #sns.heatmap(df_pairwise_overlap_pivot, annot=True, fmt='g', cmap='Blues', ax=ax)
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    #ax.set_title('Number of unique peptide pairs found by how many search engines')
    #plt.tight_layout()
    #plt.savefig(out_name + '_overlap_heatmap.svg')
    #plt.savefig(out_name + '_overlap_heatmap.png')
    #plt.close('all')

    ## cluster the heatmap
    #g = sns.clustermap(df_pairwise_overlap_pivot, annot=True, fmt='g', cmap='Blues', figsize=(12, 12))
    #g.ax_heatmap.set_xlabel('')
    #g.ax_heatmap.set_ylabel('')
    #g.ax_heatmap.set_title('Number of unique peptide pairs found by how many search engines')
    #plt.tight_layout()
    #plt.savefig(out_name + '_overlap_heatmap_clustered.svg')
    #plt.savefig(out_name + '_overlap_heatmap_clustered.png')
    #plt.close('all')

    # get df_pairwise_total in same order as df_pairwise_overlap_pivot
    df_peppair_total_sort = df_peppair_total.copy()
    df_peppair_total_sort['search_engine_lc'] = df_peppair_total_sort['search_engine'].str.lower()
    df_peppair_total_sort.sort_values(by='search_engine_lc', inplace=True)
    df_pairwise_overlap_pivot_percent_1 = df_pairwise_overlap_pivot.div(df_peppair_total_sort['total'].tolist(),axis='rows') * 100
    df_pairwise_overlap_pivot_percent_2 = df_pairwise_overlap_pivot.div(df_peppair_total_sort['total'].tolist(),axis='columns') * 100


    # same two plots but as percentage of total unique peptide pairs
    df_pairwise_overlap_pivot_percent = df_pairwise_overlap_pivot 
    
    # plot as heatmap
    #fig, ax = plt.subplots(figsize=(12, 24))
    #g = sns.clustermap(df_pairwise_overlap_pivot_percent_1, annot=True, fmt='g', cmap='Blues')
    ##g.gs.update(top=0.45, bottom=0.05)
    #g.ax_heatmap.set_xlabel('')
    #g.ax_heatmap.set_ylabel('')
    #g.ax_heatmap.set_title('Number of unique peptide pairs found by how many search engines')
    #plt.tight_layout()
    #plt.savefig(out_name + '_overlap_heatmap_percent_rows.svg')
    #plt.savefig(out_name + '_overlap_heatmap_percent_rows.png')
    #plt.close('all')

    # cluster the heatmap
    #g = sns.clustermap(df_pairwise_overlap_pivot_percent_2, annot=True, fmt='g', cmap='Blues')
    #g.ax_heatmap.set_xlabel('')
    #g.ax_heatmap.set_ylabel('')
    #g.ax_heatmap.set_title('Number of unique peptide pairs found by how many search engines')
    #plt.tight_layout()
    #plt.savefig(out_name + '_overlap_heatmap_clustered_percent_cols.svg')
    #plt.savefig(out_name + '_overlap_heatmap_clustered_percent_cols.png')
    #plt.close('all')

    # make an upsetplot of the peptidepairs found by each search engine
    search_engines = df_peppair_unique['search_engine'].unique().tolist()
    uniq_peps =  set(df_peppair_unique['pair'].unique().tolist())

    # get a list of all peptite pairs found by each search engine 
    search_engine_peps = {se: set(df_peppair_unique[df_peppair_unique['search_engine'] == se]['pair'].tolist()) for se in search_engines}
    # convert ot booleans
    search_engine_peps = [ {se: pep in search_engine_peps[se] for se in search_engines} for pep in uniq_peps]

    df_peppair_up = pd.DataFrame(search_engine_peps)
    # remove search engine unique


    # convert to boolean
    df_peppair_up = df_peppair_up.groupby(search_engines).size()
    #up.plot(df_peppair_up, orientation='horizontal', min_subset_size=0, sort_by='-degree', show_counts='%d', show_percentages=False, min_degree=1)
    #plt.savefig(out_name + '_upset.svg')
    #plt.savefig(out_name + '_upset.png')
    #plt.close('all')

    ## convert to boolean
    # df_peppair_up_nonunique = df_peppair_up[df_peppair_up.sum(axis=1) > 1]
    # df_peppair_up_nonunique = df_peppair_up_nonunique.groupby(search_engines).size()
    # up.plot(df_peppair_up, orientation='horizontal', min_subset_size=50, sort_by='-degree', show_counts='%d', show_percentages=False, min_degree=2)
    # plt.savefig(out_name + '_upset_nonunique.svg')
    # plt.savefig(out_name + '_upset_nonunique.png')
    # plt.close('all')

    # extract peptide pairs for each search engine
    #df_peppair['search_engine']
    search_engine_peps = [set(df_peppair_unique[df_peppair_unique['search_engine'] == se]['pair'].tolist()) for se in search_engines]
    # count the occurence of each peptide pair
    pep_counts = {}
    for se in search_engine_peps:
        for pep in se:
            pep_counts[pep] = pep_counts.get(pep,0) + 1

    
    # remove all peptide pairs that are not in pep_counts
    search_engine_peps_nonunique = [{pep for pep in se if pep_counts[pep] > 1} for se in search_engine_peps]
    search_engine_peps_unique = [{pep for pep in se if pep_counts[pep] == 1} for se in search_engine_peps]
    # count all and non-uniuqe peptide pairs per search engine
    search_engine_peps_unique_percent = [1-len(se[1])/len(se[0]) for se in zip(search_engine_peps, search_engine_peps_nonunique)]
    
    # count how many peptides are shared between how many search engines
    pep_shared = [0]*(len(search_engines)+1)
    for pep in pep_counts.keys():
        pep_shared[pep_counts[pep]] += 1


    # sort  by unique percent
    search_engine_peps_sort = [se for _, se in sorted(zip(search_engine_peps_unique_percent, search_engine_peps), key=lambda pair: pair[0], reverse=True)]
    search_engine_sort = [se for _, se in sorted(zip(search_engine_peps_unique_percent, search_engines), key=lambda pair: pair[0], reverse=True)]
    search_engine_peps_unique_sort = [se for _, se in sorted(zip(search_engine_peps_unique_percent, search_engine_peps_unique), key=lambda pair: pair[0], reverse=True)]
    search_engine_peps_nonunique_sort = [se for _, se in sorted(zip(search_engine_peps_unique_percent, search_engine_peps_nonunique), key=lambda pair: pair[0], reverse=True)]  
    search_engine_peps_unique_percent_sort = [0] + sorted(search_engine_peps_unique_percent, reverse=True)


    plt.figure(figsize=(24, 8))
    sv = supervenn(search_engine_peps_sort, search_engine_sort, sets_ordering=None,chunks_ordering='occurrence', min_width_for_annotation=200)
    ax= sv.axes['right_side_plot']
    #clear the plot in ax
    ax.clear()
    # plot the percentages 
    ax.barh([0] + [x+0.5 for x in range(len(search_engine_peps))], search_engine_peps_unique_percent_sort, color='gray', tick_label="")
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(out_name + '_supervenn.svg')
    plt.savefig(out_name + '_supervenn.png')
    search_eninge_summary=[]
    oberlap_summary = []
    with open(out_name + '_supervenn.txt', 'w') as f:
        #write out a tbale with number of unique and non-unique peptide pairs for each search engine
        f.write('Search engine\tPairs\tshared\nnonshared')
        for se_name, se, se_shared, se_nonshared in zip(search_engine_sort, search_engine_peps_sort, search_engine_peps_nonunique_sort, search_engine_peps_unique_sort):
            search_eninge_summary.append({
                'search_engine': se_name,
                'pairs': len(se),
                'shared': len(se_shared),
                'non_shared': len(se_nonshared)
            })
            f.write(f'{se_name}\t{len(se)}\t{len(se_shared)}\t{len(se_nonshared)} \n')
        # write out the number of peptide pairs shared between how many search engines
        f.write('\nShared\tCount\n')
        for i, count in enumerate(pep_shared):
            oberlap_summary.append({
                'overlap': i,
                'count': count
            })
            f.write(f'{i}\t{count}\n')

    plt.close('all')
    print('overlap done')
    return pd.DataFrame(search_eninge_summary), pd.DataFrame(oberlap_summary)

    


def plot_eFDR_TTfpp_over_score(df):
    """ plot the TT false positive percent and eFDR over the FDR for each search engine"""
    df = df.copy()
    # sort by search engine and score descending
    df.sort_values(by=['search_engine', 'Score'], inplace=True, ascending=[True, False])

    grp_by = df.groupby('search_engine')
    x = 2
    y = math.ceil(grp_by.ngroups / 2)

    fig, ax = plt.subplots(y, x, figsize=(12, 6*y)) #, sharex=True, sharey=True)
    all_ax = ax.flatten()
    i = 0
    # cummulative sum of EE, EH and HH for each search engine and plot
    for search_engine, group in df.groupby('search_engine'):
        thisax = all_ax[i]
        i += 1
        score = score_names[search_engine] if score_names[search_engine] in group.columns else "Score"

        scor_col = score
        if search_engine == 'Kojak':
            scor_col = "Sum_score"
        thisax.set_title(search_engine)
        if scor_col == "iProbability":
            group.sort_values(by=scor_col, inplace=True, ascending=True)
        else:
            group.sort_values(by=scor_col, inplace=True, ascending=False)
        group['EE_cs'] = group['EE'].cumsum()
        group['EH_cs'] = group['EH'].cumsum()
        group['HH_cs'] = group['HH'].cumsum()

        group['TT_cs'] = group['isTT'].cumsum()
        group['TD_cs'] = group['isTD'].cumsum()
        group['DD_cs'] = group['isDD'].cumsum()
        
        group['EE_TT_cs'] = group[group['isTT']]['EE'].cumsum()
        group['EH_TT_cs'] = group[group['isTT']]['EH'].cumsum()
        group['HH_TT_cs'] = group[group['isTT']]['HH'].cumsum()

        group['TT_Human'] = (group['EH_TT_cs'] + group['HH_TT_cs']) / group['TT_cs']
        group['eFDR'] = (group['EH_cs'] - group['HH_cs']) / group['EE_cs']
        # if highest TD_cs < DD_cs, then the FDR is calculated as (TD+DD)/TT, otherwise (TD-DD)/TT
        if group['TD_cs'].max() < group['DD_cs'].max():
            group['FDR'] = (group['TD_cs'] + group['DD_cs']) / group['TT_cs']
        else:
            group['FDR'] = (group['TD_cs'] - group['DD_cs']) / group['TT_cs']

        # make monotonious
        #if scor_col == "iProbability":
        #    group.sort_values(by=scor_col, inplace=True, ascending=True)
        #else:
        #    group.sort_values(by=scor_col, inplace=True, ascending=False)
        #group.sort_values(by=score, inplace=True)
        #group['eFDR'] = group['eFDR'].cummin()

        #group['FDR'] = group['FDR'].cummin()
        #group.reset_index(inplace=True)


        # concatenate the two dataframes
        group.to_csv(f'{search_engine}_eFDR_TTfpp.csv')
        
        p_fdr = thisax.plot(group[scor_col], group['FDR'], label='decoy FDR')
        p_efdr = thisax.plot(group[scor_col], group['eFDR'], label='entrapment FDR')
        #thisax.plot(group.loc[df['isTT']][scor_col], group.loc[df['isTT']]['FDR'], 'o')
        p_fpp = thisax.plot(group.loc[df['isTT']][scor_col], group.loc[df['isTT']]['TT_Human'], label="False Positives among TT")
        #sns.lineplot(data=df_plot, x='FDR', y='value', hue='what', ax=thisax)

        #xlim = (0, math.ceil(1000*max(0.022, group['FDR'].iloc[-1]))/1000)
        #max_yFDR = group[(group['FDR'] > xlim[0]) & (group['FDR'] < xlim[1]) & ~group['FDR'].isna()]['eFDR'].iloc[-1]
        #max_yTTfpp = group[(group['FDR'] > xlim[0]) & (group['FDR'] < xlim[1]) & ~group['FDR'].isna()]['TT_Human'].iloc[-1]
        #ylim = (0, max(max_yFDR, max_yTTfpp, 0.022))
        #thisax.set_xlim(xlim)
        #thisax.set_ylim(ylim)
        thisax.set_xlabel(scor_col)
        thisax.set_ylabel('eFDR or %Human')
        thisax.legend()#['FDR', 'entrapment FDR', 'False Positives among TT'])
        #thisax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
        thisax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    #all_ax[-1].axis('off')
    all_ax[-1].legend([p_fdr, p_efdr, p_fpp], ['FDR', 'entrapment FDR', '%Human among TT'])
    #for i in range(len(all_ax)-1):
    #    all_ax[i].legend().remove()


    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig('eFDR_TTfpp_over_FDR.svg')
    plt.savefig('eFDR_TTfpp_over_FDR.png')


    







ylim_top = {
    'Kojak': 100,
    'MangoCometXlinkProphet': 300,
    'pLink2': 100,
    'pLink2_cleav': 100,
    'pLink2_nonCleav': 100,
    'ProteinProspector': 300,
    'xiSEARCH': 50,
    'XlinkX': None,
    'XlinkX PD2.5': None,
    'XlinkX PD3.2': 150,
    'MeroX': None,
    'OpenPepXL': None
}

out_name = os.path.join('plots', 'all_scores')
out_name_ppi = os.path.join('plots', 'all_ppi_scores')
out_name_bar = os.path.join('plots', 'all_levels')

df = pd.read_csv(os.path.join('search_results', 'processed_CSM_results.csv'), low_memory=False)
df_peppair = pd.read_csv(os.path.join('search_results', 'processed_peppair_results.csv'), low_memory=False)
df_respair = pd.read_csv(os.path.join('search_results', 'processed_respair_results.csv'), low_memory=False)
df_ppi = pd.read_csv(os.path.join('search_results', 'processed_ppi_results.csv'), low_memory=False)

# rm plink2_cleav
df = df[df['search_engine'] != 'pLink2_cleav']
df_peppair = df_peppair[df_peppair['search_engine'] != 'pLink2_cleav']
df_respair = df_respair[df_respair['search_engine'] != 'pLink2_cleav']
df_ppi = df_ppi[df_ppi['search_engine'] != 'pLink2_cleav']

# rename plink2_nonCleav to plink2
df['search_engine'].loc[df['search_engine'] == 'pLink2_nonCleav'] = 'pLink2'
df_peppair['search_engine'].loc[df_peppair['search_engine'] == 'pLink2_nonCleav'] = 'pLink2'
df_respair['search_engine'].loc[df_respair['search_engine'] == 'pLink2_nonCleav'] = 'pLink2'
df_ppi['search_engine'].loc[df_ppi['search_engine'] == 'pLink2_nonCleav'] = 'pLink2'

ibaq_cutoff=1e6

possible_ppis = getPossibleProteinPairs(os.path.join('search_results', 'proteinGroups.txt'), ibaq_cutoff=ibaq_cutoff)

plot_eFDR_TTfpp_over_score(df)

df_all_scores_summary = combined_plot(df, ylim_top, out_name)
df_all_scores_summary.to_csv(out_name + ".csv")
combined_plot(df[df['search_engine'] != 'XlinkX PD3.2'], ylim_top, out_name + '_noxlx32').to_csv(out_name + "_noxlx32" + ".csv")
combined_plot(df[df['search_engine'].str.match('XlinkX PD.*')], ylim_top, out_name + '_xlx').to_csv(out_name + "_xlx" + ".csv")
combined_plot(df, ylim_top, out_name + '_dist', barplots=False).to_csv(out_name +'_dist' + ".csv")
combined_plot(df, ylim_top, out_name + "_impossible", possible_ppis=possible_ppis).to_csv(out_name + ".csv")
combined_plot_old(df_ppi, ylim_top, out_name_ppi+'_old')
combined_plot(df_ppi, ylim_top, out_name_ppi).to_csv(out_name_ppi + ".csv")
combined_plot(df_ppi, ylim_top, out_name_ppi + "_impossible", possible_ppis=possible_ppis).to_csv(out_name_ppi + ".csv")

#            "FDR(Decoy)": FDR,
#            "FDR(Entrapment)": ENTRAP_FPR,
#            "search_engine": search_engine,

# make a dataframe with fdrs type and search engine
#df_fdr_decoy=df_all_scores_summary[["search_engine", "FDR(Decoy)"]]
#df_fdr_decoy.columns = ["search_engine", "FDR"]
#df_fdr_decoy['type'] = "decoy"
#df_fdr_entrap=df_all_scores_summary[["search_engine", "FDR(Entrapment)"]]
#df_fdr_entrap.columns = ["search_engine", "FDR"]
#df_fdr_entrap['type'] = "Entrapment"
#df_fdr = pd.concat([df_fdr_decoy,df_fdr_entrap])
#sns.barplot(df_fdr, x='search_engine', y='FDR', hue='type')
#plt.rcParams['svg.fonttype'] = 'none'
#plt.savefig('Decoy_Entrap_FDR.svg')
#plt.savefig('Decoy_Entrap_FDR.png')
#plt.close('all')
                                   

# call the overlap plot
def all_overlap_plots(df, out_name,  ext, plot_TT=True, plot_TD=False, plot_DD=False):
    all_se_summary, _ = overlap_plots(df, out_name + ext, possible_ppis=None, plot_TT=plot_TT, plot_TD=plot_TD, plot_DD=plot_DD)
    #overlap_plots(df, out_name + ext + "_possilbe", possible_ppis=possible_ppis, plot_impossible=False, plot_TT=plot_TT, plot_TD=plot_TD, plot_DD=plot_DD)
    impossible_se_summary, _ = overlap_plots(df, out_name + ext + "_impossilbe", possible_ppis=possible_ppis, plot_impossible=True, plot_TT=plot_TT, plot_TD=plot_TD, plot_DD=plot_DD)
    #overlap_plots(df, out_name + ext + "_EEOnly", possible_ppis=None, EEonly=True, plot_TT=plot_TT, plot_TD=plot_TD, plot_DD=plot_DD)
    #overlap_plots(df, out_name + ext + "_EEOnly_impossible", possible_ppis=possible_ppis, EEonly=True, plot_impossible=True, plot_TT=plot_TT, plot_TD=plot_TD, plot_DD=plot_DD)
    # sort impossible in same search_engine order as all
    impossible_se_summary = impossible_se_summary.set_index('search_engine').reindex(all_se_summary['search_engine']).reset_index()
    # plot percent of impossible among the unique for each search engine
    all_se_summary['non_shared_impossible'] = impossible_se_summary['non_shared']
    all_se_summary['shared_impossible'] = impossible_se_summary['shared']
    all_se_summary['shared_impossible_percent'] = all_se_summary['shared_impossible'] / all_se_summary['shared'] * 100
    all_se_summary['non_shared_impossible_percent'] = all_se_summary['non_shared_impossible'] / all_se_summary['non_shared'] * 100
    all_se_summary.to_csv(out_name + ext + ".csv")
    plt.figure(figsize=(3.6, 4.3))
    sns.barplot(data=all_se_summary, x='search_engine', y='non_shared_impossible_percent', color='#d12c2c')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(out_name + ext + "_impossible_percent.svg")
    plt.savefig(out_name + ext + "_impossible_percent.png")
    plt.close('all')





all_overlap_plots(df, out_name, "_CSM")
all_overlap_plots(df_peppair, out_name, "_peppair")
all_overlap_plots(df_respair, out_name, "_respair")
all_overlap_plots(df_ppi, out_name, "_ppi")

all_overlap_plots(df_peppair[df_peppair['search_engine'] != 'XlinkX PD3.2'], out_name, "_peppair_noxlx32")
all_overlap_plots(df_peppair[df_peppair['search_engine'].str.match('XlinkX PD.*')], out_name, "_peppair_xlx")

ext = "_peppair"
#overlap_plots(df_peppair, 'TD' + out_name + ext + "_possilbe", possible_ppis=possible_ppis, plot_impossible=False, plot_TT=False, plot_TD=True, plot_DD=False)
#overlap_plots(df_peppair, 'DD' + out_name + ext + "_possilbe", possible_ppis=possible_ppis, plot_impossible=False, plot_TT=False, plot_TD=False, plot_DD=True)
#overlap_plots(df_peppair,  out_name.replace("/", "/TD_") + ext + "_impossilbe", possible_ppis=possible_ppis, plot_impossible=True, plot_TT=False, plot_TD=True, plot_DD=False)
#overlap_plots(df_peppair, out_name.replace("/", "/DD_") + ext + "_impossilbe", possible_ppis=possible_ppis, plot_impossible=True, plot_TT=False, plot_TD=False, plot_DD=True)
#overlap_plots(df_peppair, out_name.replace("/", "/TD_") + ext, possible_ppis=None, plot_impossible=False, plot_TT=False, plot_TD=False, plot_DD=True)
#overlap_plots(df_peppair, out_name.replace("/", "/DD_") + ext, possible_ppis=None, plot_impossible=False, plot_TT=False, plot_TD=True, plot_DD=False)



# call the overlap plot
overlap_plots(df, out_name + "_possiblePPI", possible_ppis=possible_ppis)

bar_plots(df, df_peppair, df_respair, df_ppi, out_name_bar, possible_ppis=possible_ppis)

#out_name_ppi = os.path.join('plots', 'all_ppi_scores_CSM_to_ppi')
#out_name_bar = os.path.join('plots', 'ppifdr_from_CSM')

# ToDo: no used for now
#df_ppi2csm = pd.read_csv(os.path.join('search_results', 'ppi_processed_CSM_results.csv'))
#df_ppi2peppair = pd.read_csv(os.path.join('search_results', 'ppi_processed_peppair_results.csv'))
#df_ppi2respair = pd.read_csv(os.path.join('search_results', 'ppi_processed_respair_results.csv'))
#df_csm_to_2ppi = pd.read_csv(os.path.join('search_results', 'processed_CSM_TO_2PercentPPI_results.csv'))

#combined_plot(df_ppi2csm, ylim_top, out_name).to_csv(out_name + ".csv")
#combined_plot(df_csm_to_2ppi, ylim_top, out_name_ppi).to_csv(out_name_ppi + ".csv")

#bar_plots(None, None, None, df_csm_to_2ppi, out_name_bar, possible_ppis=possible_ppis)



out_name_ppi = os.path.join('plots', 'all_ppi_scores_Highest_to_ppi')
out_name_bar = os.path.join('plots', 'ppifdr_from_highest')

df_highest_to_ppi = pd.read_csv(os.path.join('search_results', 'processed_highest_to_ppi_results.csv'))

#combined_plot(df_ppi2csm, ylim_top, out_name).to_csv(out_name + ".csv")
#combined_plot(df_csm_to_2ppi, ylim_top, out_name_ppi).to_csv(out_name_ppi + ".csv")

bar_plots(None, None, None, df_highest_to_ppi, out_name_bar, possible_ppis=possible_ppis)

# make a combined bar plot 

# read in the numbers for the csm and the highest to ppi bar plots
df_csm = pd.read_csv(os.path.join('plots', 'all_levels_numbers.csv'))
df_highest = pd.read_csv(os.path.join('plots', 'ppifdr_from_highest_numbers.csv'))

# reduce to csm and ppi level  
#df_csm = df_csm[df_csm['level'].isin(['CSM', 'PPI'])]
df_csm = df_csm[df_csm['level'].isin(['CSM'])]
df_csm.loc[df_csm['level'] == 'PPI', 'level'] = 'CSM to PPI'

df_highest.loc[df_highest['level'] == 'PPI', 'level'] = 'Highest to PPI'
df_highest.loc[df_highest['level'] == 'FDR', 'level'] = 'Highest to PPI'

# combine the two dataframes
df_csm_highest = pd.concat([df_csm, df_highest], ignore_index=True)
# remove plink_cleav
df_csm_highest = df_csm_highest[df_csm_highest['search_engine'] != 'pLink2_cleav']
df_csm_highest.loc[df_csm_highest['search_engine'] == 'pLink2_nonCleav', 'search_engine'] = 'pLink2'
df_csm_highest.to_csv(os.path.join('plots', 'combined_numbers' + str(ibaq_cutoff) + '.csv'), index=False)


# transform the dataframes to long format from FDR and %known FP
df_csm_highest_vert = pd.melt(df_csm_highest, id_vars=['search_engine', 'level'], value_vars=['FDR', '%Known FP'], var_name='what', value_name='value')

# seaborn map plot a bar plot for each search engine
matplotlib.rcParams['svg.fonttype'] = 'none'
g = sns.FacetGrid(df_csm_highest_vert, col="search_engine", col_wrap=2, sharey=True, sharex=True, height=3, aspect=1.5)
g.set_titles("{col_name}")
g.map(sns.barplot,'level', 'value', 'what', order=df_csm_highest_vert['level'].unique(), palette=sns.color_palette(['#BEAED4', '#fdc086']))

for ax in g.axes.flat:
    ax.axhline(y=0.02, color='black', linestyle=':')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    #ax.title.set_position([0.5, -10])
    #ax.title.set_pad(10)

# add a custom legend that includes the 2% FDR line
legends = [ mpatches.Patch(color='#BEAED4', label='Decoy based error*'),
            mpatches.Patch(color='#fdc086', label='Total observable False Positives among TT'), 
            mlines.Line2D([], [], color='black', linestyle=':', label='2% FDR')]
plt.legend(handles=legends, loc='upper center', bbox_to_anchor=(0, -0.1),
          fancybox=True, shadow=True, ncol=5)

# make some space between the plots so titles don't look like they are axis labels
plt.subplots_adjust(hspace=0.4, wspace=0.1)

plt.savefig(os.path.join('plots', 'combined' + str(ibaq_cutoff) + '.svg'))
plt.savefig(os.path.join('plots', 'combined' + str(ibaq_cutoff) + '.png'))
plt.close('all')

# progress bar style plots of known %TP
plt.barh(df_csm_highest['search_engine'], df_csm_highest['True Positives%<'], color='#fdc086')
plt.tight_layout()
plt.savefig(os.path.join('plots', 'combined_maxTP_percent_' + ('%.2E' % ibaq_cutoff) + '.svg'))
plt.savefig(os.path.join('plots', 'combined_maxTP_percent_' + ('%.2E' % ibaq_cutoff) + '.png'))

# break long search eninge names
df_csm_highest.loc[df_csm_highest['search_engine']== 'MangoCometXlinkProphet',  'search_engine'] = 'MangoComet\nXlinkProphet'
df_csm_highest.loc[df_csm_highest['search_engine']== 'ProteinProspector',  'search_engine'] = 'Protein\nProspector'

plt.close('all')

def progress_bar_plot(df_csm_highest, outname = 'ppi_progress_maxTP_percent'):
    # double plot of known %TP and %FP
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    # sort by 'True Positives%<'
    df_csm_highest.sort_values(by='True Positives%<', inplace=True)
    csm_level = 'CSM'
    ppi_level = 'Highest to PPI'
    for a, what in enumerate([csm_level, ppi_level]):
        ax1 = ax[a]
        ax1.set_title('Maximum True Positives\n' + what)
        df = df_csm_highest[df_csm_highest['level'] == what]
        ax1.barh(df['search_engine'], [1]*len(df['search_engine']), color='#fd0000')
        ax1.barh(df['search_engine'], df['True Positives%<'], color='#fdc086')
        #ax2 = ax[a].twinx()
        #ax2.set_ylim(ax1.get_ylim())
        # add a 98% line
        ax1.axvline(x=0.98, color='black', linestyle=':')
        #ax1.set_yticks(
        #    np.arange(len(df['search_engine'])),
        #    labels=df['search_engine'])
        # write the percent values at the end of the bars
        if what == csm_level:
            for i, v in enumerate(df['True Positives%<']):
                ax1.text(v - 0.01, i - 0.1, ('%.1f%%' % (v*100)), color='black', ha='right')
                #ax1.text(v + 0.01, i - 0.1, ('%.1f%%' % ((1-v)*100)), color='gray', ha='left')
        else:
            for i, (v, se) in enumerate(zip(df['True Positives%<'],df['search_engine'])):
                # fo xsiearchh, mang and kojak we have ppi - all others get marked with a '*'
                if se in ['xiSEARCH', 'MangoComet\nXlinkProphet', 'Kojak', 'Protein\nProspector']:
                    ax1.text(v - 0.01, i - 0.1, ('%.1f%%' % (v*100)), color='black', ha='right')
                    #ax1.text(v + 0.01, i - 0.1, ('%.1f%%' % ((1-v)*100)), color='gray', ha='left')
                else:
                    ax1.text(v - 0.01, i - 0.1, ('%.1f%%*' % (v*100)), color='black', ha='right')
                    #ax1.text(v + 0.01, i - 0.1, ('%.1f%%*' % ((1-v)*100)), color='gray', ha='left')
            # overwrite also the yticks with the search engine names + '*'
            ax1.set_yticks(
                range(len(df['search_engine'])),
                labels=[se + '*' if se not in ['xiSEARCH', 'MangoComet\nXlinkProphet', 'Kojak'] else se for se in df['search_engine']])
            
        #ax2.set_yticks(ax1.get_yticks(),
        #    labels=['%.1f%%' % (x*100) for x in df['True Positives%<']])
    ax[0].text(-0.2, 1.1, "A", transform=ax1.transAxes, size=15, weight='bold')
    ax[1].text(-0.2, 1.1, "B", transform=ax1.transAxes, size=15, weight='bold')

    fig.tight_layout()
    plt.savefig(os.path.join('plots', outname + '_' + ('%.2E' % ibaq_cutoff) + '.svg'))
    plt.savefig(os.path.join('plots', outname + '_' + ('%.2E' % ibaq_cutoff) + '.png'))

progress_bar_plot(df_csm_highest)
progress_bar_plot(df_csm_highest[df_csm_highest['search_engine'].str.match("XlinkX.*")] , "ppi_progress_maxTP_percent_xlinkx")


df_csm_highest_vert = pd.melt(df_csm_highest, id_vars=['search_engine', 'level'], value_vars=['True Positives%<', '%Known FP'], var_name='what', value_name='value')
g = sns.FacetGrid(df_csm_highest_vert, col="level", col_wrap=2, sharey=True, sharex=True, height=3, aspect=1.5)
g.set_titles("{col_name}")
g.map(sns.barplot,'value', 'search_engine', 'what', order=df_csm_highest_vert['search_engine'].unique())
plt.savefig(os.path.join('plots', 'combined_maxTP_percent_' + ('%.2E' % ibaq_cutoff) + '_2.svg'))
plt.savefig(os.path.join('plots', 'combined_maxTP_percent_' + ('%.2E' % ibaq_cutoff) + '_2.png'))



plt.close('all')

df_sorted = df_csm.sort_values(by='%Known FP')
df_sorted.loc[df_sorted['search_engine']== 'MangoCometXlinkProphet',  'search_engine'] = 'MangoComet\nXlinkProphet'
df_sorted.loc[df_sorted['search_engine']== 'ProteinProspector',  'search_engine'] = 'Protein\nProspector'
df_sorted['factor'] = df_sorted['%Known FP'] / 0.02
fig, ax = plt.subplots(1, 1, figsize=(6, 6))

plt.title('Observeable Known False Positives CSMs\nat 2%FDR Decoy based')
plt.bar(df_sorted['search_engine'], df_sorted['%Known FP'], color='#fdc086')
locs, labels = plt.xticks()
plt.setp(labels, rotation=80)
# add a custom legend that includes the 2% FDR line
legends = [mpatches.Patch(color='#fdc086', label='Total observable False Positives among TT'), 
            mlines.Line2D([], [], color='black', linestyle=':', label='2% FDR (Target)')]
plt.legend(handles=legends)
ax.axhline(y=0.02, color='black', linestyle=':')
ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))

fig.tight_layout()

plt.savefig(os.path.join('plots', 'KnownFP_CSM_percent_' + ('%.2E' % ibaq_cutoff) + '_2.svg'))
plt.savefig(os.path.join('plots', 'KnownFP_CSM_percent_' + ('%.2E' % ibaq_cutoff) + '_2.png'))

plt.close('all')
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
plt.title('Observeable Known False Positives CSMs\nat 2%FDR Decoy based')
plt.bar(df_sorted['search_engine'], df_sorted['factor'], color='#fdc086')
locs, labels = plt.xticks()
plt.setp(labels, rotation=80)
# add a custom legend that includes the 2% FDR line
legends = [mpatches.Patch(color='#fdc086', label='Total observable False Positives among TT'), 
            mlines.Line2D([], [], color='black', linestyle=':', label='2% FDR (Target)')]
plt.legend(handles=legends)
ax.axhline(y=1, color='black', linestyle=':')

fig.tight_layout()

plt.savefig(os.path.join('plots', 'KnownFP_CSM_factor_' + ('%.2E' % ibaq_cutoff) + '_2.svg'))
plt.savefig(os.path.join('plots', 'KnownFP_CSM_factor_' + ('%.2E' % ibaq_cutoff) + '_2.png'))




print("done")
