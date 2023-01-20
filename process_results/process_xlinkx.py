from plots_and_functions import *

fasta_1 = fasta_to_dict('database/uniprot-k12-filtered-proteome_UP000000625.fasta')
all_proteins = list(fasta_1.keys())

df_targets = pd.read_csv('../xlinkx/2p_all_files_plus3_CSMs.txt', sep='\t')
df_decoys = pd.read_csv('../xlinkx/2p_all_files_plus3_DecoyCSMs.txt', sep='\t')
df_targets = df_targets[df_targets['Crosslink Type'].str.contains('Inter')]
df_decoys = df_decoys[df_decoys['Crosslink Type'].str.contains('Inter')]
df_targets['E1'] = df_targets['Protein Accession A'].apply(check_amb, all_proteins=all_proteins)
df_targets['E2'] = df_targets['Protein Accession B'].apply(check_amb, all_proteins=all_proteins)
df_targets['EE'] = df_targets['E1'] & df_targets['E2']
df_targets['EH'] = (df_targets['E1'] & (~df_targets['E2']) | (~df_targets['E1']) & df_targets['E2'])
df_targets['HH'] = (~df_targets['E1']) & (~df_targets['E2'])
print(sum(df_targets.HH | df_targets.EH)/len(df_targets))

df_targets['entr_group'] = df_targets['E1'].astype(int) + df_targets['E2'].astype(int)
df_targets['entr_group'] = df_targets['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})
df_decoys['entr_group'] = 'decoy'
df_plot = pd.concat([df_targets, df_decoys])

ax = plot_distribution(df=df_plot.reset_index(), x='XlinkX Score', bins=50)
plt.show()

summary_table = df_plot.reset_index()['entr_group'].value_counts()
summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
summary_table.to_csv('xlinkx_numbers.csv')
