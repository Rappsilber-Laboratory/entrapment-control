from plots_and_functions import *


fasta_1 = fasta_to_dict('database/uniprot-k12-filtered-proteome_UP000000625.fasta')
# read in uniprot IDs of Ecoli as list
all_proteins = list(fasta_1.keys())

df = pd.read_csv('/home/swantje/Dropbox/Arbeit/other_decoys/xisearch/'
                 'CMS_FDR_proteome+human3_linMod_filtered_0_PSM_xiFDR1.3.36.csv')
# subset to heteromeric
df = df[df.fdrGroup.str.contains('etween')]

# are protein 1 or protein 2 from ecoli --> gives true / false in separate column
# if ambigous match contains an ecoli protein return True
df['E1'] = df['Protein1'].apply(check_amb, all_proteins=all_proteins)
df['E2'] = df['Protein2'].apply(check_amb, all_proteins=all_proteins)

# add columns differentiating Ecoli-Ecoli, Ecoli-human, human-human
# is this used later or leftover?
df['EE'] = df['E1'] & df['E2']
df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
df['HH'] = (~df['E1']) & (~df['E2'])
print(sum(df.HH | df.EH)/len(df))

# add group column for plotting purposes
df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

# mark TD and DD as decoys
df.loc[((df['isTD'])| (df['isDD'])), 'entr_group'] = 'decoy'

ax = plot_distribution(df=df, x='Score', bins=50)
plt.show()

ax = plot_distribution(df=df, x='Score', bins=50, ylim_top=50)
plt.show()

summary_table = df.reset_index()['entr_group'].value_counts()
summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
# summary_table.to_csv('xisearch.csv')
