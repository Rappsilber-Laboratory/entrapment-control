from plots_and_functions import plot_distribution
import pandas as pd
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
# import numpy as np
from score_names import score_names
import seaborn as sns





def combined_plot(df, ylim_top, out_name):
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
    fig, ax = plt.subplots(int(df['search_engine'].nunique()/2), 4, figsize=(12, 12),
                        gridspec_kw={'width_ratios': [1, 0.1, 1, 0.1]})
    ax = ax.flatten()
    i = 0
    summary = []
    for search_engine, group in df.groupby('search_engine', sort=False):

        xcol = score_names[search_engine] if score_names[search_engine] in group.columns else "Score"
        plot_distribution(df=group, x=xcol, bins=50, ax=ax[i], ylim_top=ylim_top[search_engine])
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
        TT_EH = sum(group['EH'][group['isTT'] == True])
        TT_HH = sum(group['HH'][group['isTT'] == True])
        # kojak has no decoys for the entrapment db so need to adjust numbers here
        if search_engine.lower() == 'kojak':
            # remove all entrapment from TD and TT
            TD = sum(group['isTD'][~(group['EH'] | group['HH'])])
            TT = sum(group['isTT'][~(group['EH'] | group['HH'])])
        else: 
            TD = sum(group['isTD'])
            TT = sum(group['isTT'])

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
        ENTRAP_FP = ENTRAP_TOTAL + (TT_EH - TT_HH)

        summary.append({
            "search_engine": search_engine,
            "TT": TT,
            "TD": TD,
            "DD": DD,
            "EH(TT)": TT_EH,
            "HH(TT)": TT_HH,
            "FDR(Decoy)": FDR,
            "FDR(Entrapment)": ENTRAP_FP/TT
            })

        decoy_based_error = FDR * 100
        entr_based_error = ENTRAP_FP / TT * 100

        sns.barplot(x=['D'+note, 'E'], y=[decoy_based_error, entr_based_error], ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#EB92CA']))
        sns.barplot(x=['D'+note, 'E'], y=[decoy_based_error, ENTRAP_TOTAL / TT * 100], ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#E364B4']))
        ax[i].spines[['right', 'top']].set_visible(False)
        if entr_based_error < 8:
            ax[i].set_ylim(0, 8)
        else:
            ax[i].set_ylim(0, round(entr_based_error))
        ax[i].axhline(2, color='black', linestyle='--')
        ax[i].yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax[i].yaxis.set_minor_locator(MaxNLocator(nbins=8))
        ax[i].set_ylabel('error [%]')
        i += 1

    plt.tight_layout()
    # plt.show()

    plt.savefig(out_name + '.svg')
    plt.savefig(out_name + '.png')
    plt.close('all')
    return pd.DataFrame(summary)


def bar_plots(df_csm, df_peppair, df_respair, df_ppi, out_name):
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
    df_csm['level'] = 'CSM'
    df_peppair['level'] = 'Peptide Pair'
    df_respair['level'] = 'Residue Pair'
    df_ppi['level'] = 'PPI'
    df_csm['level_number'] = 0
    df_peppair['level_number'] = 1
    df_respair['level_number'] = 2
    df_ppi['level_number'] = 3

    # combine all dataframes
    df = pd.concat([df_csm, df_peppair, df_respair, df_ppi], ignore_index=True)
    df['search_engine_lc'] = df['search_engine'].str.lower()
    df.sort_values(by=['search_engine_lc', "level_number"], inplace=True)



    # combined plot: all scores in one plot
    fig, ax = plt.subplots(int(df['search_engine'].nunique()/2), 2, figsize=(12, 12))
    ax = ax.flatten()
    i = 0
    summary = []
    for search_engine, group in df.groupby('search_engine', sort=False):
        if i == 0:
            handles, labels = ax[i].get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper center')
        else:
            ax[i].legend().remove()

        # summarize results for each level
        results_e_observed = []
        results_e_estimate = []
        for level, level_group in group.groupby('level', sort=False):
            level_group = level_group.reset_index()
            TT_EH = sum(level_group['EH'][level_group['isTT'] == True])
            TT_HH = sum(level_group['HH'][level_group['isTT'] == True])
            DD = sum(level_group['isDD'])
            if search_engine.lower() == 'kojak':
                # remove all entrapment from TD and TT
                TD = sum(level_group['isTD'][~(level_group['EH'] | level_group['HH'])])
                TT = sum(level_group['isTT'][~(level_group['EH'] | level_group['HH'])])
            else: 
                TD = sum(level_group['isTD'])
                TT = sum(level_group['isTT'])
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
                "order": level_number*10+1,
                "what": "Observed entrapment",
                "value": ENTRAP_TOTAL/TT*100
            },
            {
                "level": level,
                "order": level_number*10,
                "what": fdr_name,
                "value": FDR*100
            }])        
            results_e_estimate.extend([{
                "level": level,
                "order": level_number*10+1,
                "what": "Entrapment based error",
                "value": ENTRAP_FP/TT*100
            },
            {
                "level": level,
                "order": level_number*10,
                "what": "",
                "value": 0
            }            
            ])

        # convert to dataframe
        df_results_e_observed = pd.DataFrame(results_e_observed)
        df_results_e_estimate = pd.DataFrame(results_e_estimate)

        df_results_e_observed.sort_values(by='order', inplace=True)
        df_results_e_estimate.sort_values(by='order', inplace=True)
                
        sns.barplot(data=df_results_e_estimate, x='level', y='value', hue="what", ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#EB92CA']))
        sns.barplot(data=df_results_e_observed, x='level', y='value', hue="what", ax=ax[i],
                    palette=sns.color_palette(['#83320C', '#E364B4']))
        ax[i].spines[['right', 'top']].set_visible(False)
        #if entr_based_error < 8:
        #    ax[i].set_ylim(0, 8)
        #else:
        #    ax[i].set_ylim(0, round(entr_based_error))
        ax[i].axhline(2, color='black', linestyle='--')
        ax[i].set_xlabel(search_engine)
        #ax[i].yaxis.set_major_locator(MaxNLocator(nbins=4))
        #ax[i].yaxis.set_minor_locator(MaxNLocator(nbins=8))
        ax[i].set_ylabel('error [%]')
        i += 1

    plt.tight_layout()
    # plt.show()

    plt.savefig(out_name + '.svg')
    plt.savefig(out_name + '.png')
    plt.close('all')
    return pd.DataFrame(summary)


ylim_top = {
    'Kojak': 100,
    'MangoCometXlinkProphet': 300,
    'pLink2': 120,
    'ProteinProspector': 300,
    'xiSEARCH': 50,
    'XlinkX': None,
    'MeroX': None,
    'OpenPepXL': None
}

out_name = os.path.join('plots', 'all_scores')
out_name_ppi = os.path.join('plots', 'all_ppi_scores')
out_name_bar = os.path.join('plots', 'all_levels')

df = pd.read_csv(os.path.join('search_results', 'processed_CSM_results.csv'))
df_peppair = pd.read_csv(os.path.join('search_results', 'processed_peppair_results.csv'))
df_respair = pd.read_csv(os.path.join('search_results', 'processed_respair_results.csv'))
df_ppi = pd.read_csv(os.path.join('search_results', 'processed_ppi_results.csv'))

combined_plot(df, ylim_top, out_name).to_csv(out_name + ".csv")
combined_plot(df_ppi, ylim_top, out_name_ppi).to_csv(out_name_ppi + ".csv")

bar_plots(df, df_peppair, df_respair, df_ppi, out_name_bar)

out_name = os.path.join('plots', 'all_scores_from_ppi')
out_name_ppi = os.path.join('plots', 'all_ppi_scores_from_ppi')
out_name_bar = os.path.join('plots', 'all_levels_from_ppi')

df_ppi2csm = pd.read_csv(os.path.join('search_results', 'ppi_processed_CSM_results.csv'))
df_ppi2peppair = pd.read_csv(os.path.join('search_results', 'processed_peppair_results.csv'))
df_ppi2respair = pd.read_csv(os.path.join('search_results', 'processed_respair_results.csv'))
df_ppi2ppi = pd.read_csv(os.path.join('search_results', 'ppi_processed_ppi_results.csv'))

combined_plot(df_ppi2csm, ylim_top, out_name).to_csv(out_name + ".csv")
combined_plot(df_ppi2ppi, ylim_top, out_name_ppi).to_csv(out_name_ppi + ".csv")

bar_plots(df_ppi2csm, df_ppi2peppair, df_ppi2respair, df_ppi2ppi, out_name_bar)
