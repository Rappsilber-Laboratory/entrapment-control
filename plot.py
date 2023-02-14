from plots_and_functions import plot_distribution
import pandas as pd
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
# import numpy as np
from score_names import score_names
import seaborn as sns


# read in the results
df = pd.read_csv(os.path.join('search_results', 'processed_results.csv'))

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


# # single plots
# # Kojak
# # score = 'E_value'
# score = 'K_score'
# # score = 'iProbability'
# kojak_df = df[df['search_engine'] == 'Kojak']
# plot_distribution(df=kojak_df, x=score, bins=50)
# plt.show()
# plot_distribution(df=kojak_df, x=score, bins=50, ylim_top=100)
# plt.show()
#
# # Mango Comet XlinkProphet
# mango_df = df[df['search_engine'] == 'MangoCometXlinkProphet']
# plot_distribution(df=mango_df, x='score', bins=50)
# plt.show()
# plot_distribution(df=mango_df, x='score', bins=50, ylim_top=300)
# plt.show()
#
# # pLink2
# plink2_df = df[df['search_engine'] == 'pLink2']
# bins = np.arange(plink2_df['score'].min(), plink2_df['score'].max(), 0.025)
# plot_distribution(df=plink2_df, x='score', bins=bins)
# plt.show()
# plot_distribution(df=plink2_df, x='score', bins=bins, ylim_top=120)
# plt.show()
#
# # proteinProspector
# pp_df = df[df['search_engine'] == 'ProteinProspector']
# plot_distribution(df=pp_df, x='score', bins=50)
# plt.show()
# plot_distribution(df=pp_df, x='score', bins=50, ylim_top=300)
# plt.show()
#
# # xiSEARCH
# xisearch_df = df[df['search_engine'] == 'xiSEARCH']
# plot_distribution(df=xisearch_df, x='score', bins=50)
# plt.show()
# plot_distribution(df=xisearch_df, x='score', bins=50, ylim_top=50)
# plt.show()
#
# # XlinkX
# xlinkx_df = df[df['search_engine'] == 'XlinkX']
# plot_distribution(df=xlinkx_df, x='score', bins=50)
# plt.show()
#
# # MeroX
# merox_df = df[df['search_engine'] == 'MeroX']
# ax = plot_distribution(df=merox_df, x='score', bins=50)
# ax.axvline(50, color='black')  # in default settings, this is the minimal score cutoff
# plt.tight_layout()
# plt.show()
#
# # OpenPepXL
# openpepxl_df = df[df['search_engine'] == 'OpenPepXL']
# plot_distribution(df=openpepxl_df, x='score', bins=50)
# plt.tight_layout()
# plt.show()

# combined plot: all scores in one plot
df.sort_values(by='search_engine', inplace=True, key=lambda x: x.str.lower())
fig, ax = plt.subplots(int(df['search_engine'].nunique()/2), 4, figsize=(12, 12),
                       gridspec_kw={'width_ratios': [1, 0.1, 1, 0.1]})
ax = ax.flatten()
i = 0
for search_engine, group in df.groupby('search_engine', sort=False):
    plot_distribution(df=group, x=score_names[search_engine], bins=50, ax=ax[i], ylim_top=ylim_top[search_engine])
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

    # summary table
    total = group.shape[0]
    summary_table = group.reset_index()['entr_group'].value_counts()

    decoy_based_error = summary_table['decoy'] * 100 / total
    entr_based_error = summary_table['entrapment'] * 100 / total
    sns.barplot(x=['D', 'E'], y=[decoy_based_error, entr_based_error], ax=ax[i],
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

out_name = os.path.join('plots', 'all_scores')
plt.savefig(out_name + '.svg')
plt.savefig(out_name + '.png')
plt.close('all')
