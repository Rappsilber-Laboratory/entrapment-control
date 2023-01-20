import numpy as np
from plots_and_functions import *


fasta_1 = fasta_to_dict('database/uniprot-k12-filtered-proteome_UP000000625.fasta')
all_proteins = list(fasta_1.keys())

df = pd.read_csv('/home/swantje/Dropbox/Arbeit/other_decoys/plink2/DSSO_HumanDB_plinkFDR_xifdr.csv')
df['1 - Score'] = df['Score'].apply(lambda x: 1. - x)
df = df[(df['Group'] == 'between')]

df_with_dec = pd.read_csv('/home/swantje/Dropbox/Arbeit/other_decoys/plink2/DSSO_HumanDB_unfiltered_xifdr.csv')
df_with_dec['1 - Score'] = df_with_dec['Score'].apply(lambda x: 1. - x)
df_with_dec = df_with_dec[(df_with_dec['Group'] == 'between')]

df = pd.concat([df, df_with_dec[(~df_with_dec['isTT']) & (df_with_dec['1 - Score'] >= df['1 - Score'].min())]])

df['E1'] = df['Protein1'].apply(check_amb, all_proteins=all_proteins)
df['E2'] = df['Protein2'].apply(check_amb, all_proteins=all_proteins)

df['EE'] = df['E1'] & df['E2']
df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
df['HH'] = (~df['E1']) & (~df['E2'])
print(sum(df.HH | df.EH)/len(df))

df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})
df.loc[((df['isTD'])| (df['isDD'])), 'entr_group'] = 'decoy'

ax = plot_distribution(df=df, x='1 - Score', bins=np.arange(df['1 - Score'].min(), df['1 - Score'].max(), 0.025))
plt.show()

ax = plot_distribution(df=df, x='1 - Score', bins=np.arange(df['1 - Score'].min(), df['1 - Score'].max(), 0.025),
                       ylim_top=120)
plt.show()

summary_table = df.reset_index()['entr_group'].value_counts()
summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
summary_table.to_csv('plink2.csv')
