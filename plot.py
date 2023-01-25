from plots_and_functions import plot_distribution
import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np

# read in the results
df = pd.read_csv(os.path.join('search_results', 'processed_results.csv'))

# single plots
# Kojak
# score = 'E_value'
score = 'K_score'
# score = 'iProbability'
kojak_df = df[df['search_engine'] == 'Kojak']
plot_distribution(df=kojak_df, x=score, bins=50)
plt.show()
plot_distribution(df=kojak_df, x=score, bins=50, ylim_top=100)
plt.show()

# Mango Comet XlinkProphet
mango_df = df[df['search_engine'] == 'MangoCometXlinkProphet']
plot_distribution(df=mango_df, x='score', bins=50)
plt.show()
plot_distribution(df=mango_df, x='score', bins=50, ylim_top=300)
plt.show()

# pLink2
plink2_df = df[df['search_engine'] == 'pLink2']
bins = np.arange(plink2_df['score'].min(), plink2_df['score'].max(), 0.025)
plot_distribution(df=plink2_df, x='score', bins=bins)
plt.show()
plot_distribution(df=plink2_df, x='score', bins=bins, ylim_top=120)
plt.show()

# proteinProspector
pp_df = df[df['search_engine'] == 'ProteinProspector']
plot_distribution(df=pp_df, x='score', bins=50)
plt.show()
plot_distribution(df=pp_df, x='score', bins=50, ylim_top=300)
plt.show()

# xiSEARCH
xisearch_df = df[df['search_engine'] == 'xiSEARCH']
plot_distribution(df=xisearch_df, x='score', bins=50)
plt.show()
plot_distribution(df=xisearch_df, x='score', bins=50, ylim_top=50)
plt.show()

# XlinkX
xlinkx_df = df[df['search_engine'] == 'XlinkX']
ax = plot_distribution(df=xlinkx_df, x='score', bins=50)
plt.show()

# MeroX
merox_df = df[df['search_engine'] == 'MeroX']
ax = plot_distribution(df=merox_df, x='score', bins=50)
ax.axvline(50, color='black')  # in default settings, this is the minimal score cutoff
plt.tight_layout()
plt.show()
