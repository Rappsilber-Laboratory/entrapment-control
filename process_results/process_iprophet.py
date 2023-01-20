from plots_and_functions import *


fasta_1 = fasta_to_dict(
    '/home/swantje/Dropbox/Arbeit/other_decoys/uniprot-k12-filtered-proteome_UP000000625.fasta')
all_proteins = list(fasta_1.keys())

df = pd.read_excel('/home/swantje/Dropbox/Arbeit/other_decoys/iprophet/iprophet-withdecoys-xl-filtered_2pct.xlsx')
df = df[df.intra == 0]
df['E1'] = df['protein1'].apply(check_amb, all_proteins=all_proteins)
df['E2'] = df['protein2'].apply(check_amb, all_proteins=all_proteins)

df['EE'] = df['E1'] & df['E2']
df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
df['HH'] = (~df['E1']) & (~df['E2'])
print(sum(df.HH | df.EH)/len(df))

df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})
df.loc[((df['decoy'] == 1)| (df['decoy-decoy'] == 1)), 'entr_group'] = 'decoy'

ax = plot_distribution(df=df, x='probability', bins=50)
plt.show()

ax = plot_distribution(df=df, x='probability', bins=50, ylim_top=300)
plt.show()

summary_table = df.reset_index()['entr_group'].value_counts()
summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
summary_table.to_csv('iprophet.csv')
