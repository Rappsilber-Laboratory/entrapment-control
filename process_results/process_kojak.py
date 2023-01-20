from plots_and_functions import *


def check_amb_fdr_group(x):
    prot1, prot2 = x['alpha_Proteins'], x['beta_Proteins']
    if ';' in prot1:
        proteins_1 = prot1.split(';')
        for prot1_i in proteins_1:
            if ';' in prot2:
                proteins_2 = prot2.split(';')
                if prot1_i in proteins_2:
                    return 'self'
                else:
                    return 'between'
            else:
                if prot1_i == prot2:
                    return 'self'
                else:
                    return 'between'
    elif ';' in prot2:
        proteins_2 = prot2.split(';')
        for prot2_i in proteins_2:
            if ';' in prot1:
                proteins_1 = prot1.split(';')
                if prot2_i in proteins_1:
                    return 'self'
                else:
                    return 'between'
            else:
                if prot2_i == prot1:
                    return 'self'
                else:
                    return 'between'
    else:
        if prot1 == prot2:
            return 'self'
        else:
            return 'between'


fasta_1 = fasta_to_dict(
    '/home/swantje/Dropbox/Arbeit/other_decoys/uniprot-k12-filtered-proteome_UP000000625.fasta')
all_proteins = list(fasta_1.keys())

df = pd.read_csv('/home/swantje/Dropbox/Arbeit/other_decoys/kojak/Kojak_CSMs_NoLoopLinks.tsv', sep='\t')
df['decoy1'] = df['alpha_Proteins'].str.contains('DECOY')
df['decoy2'] = df['beta_Proteins'].str.contains('DECOY') # if ambigous will count as decoy

# TODO subset to between!
df['fdrGroup'] = df.apply(check_amb_fdr_group, axis=1)
df = df[df['fdrGroup'] == 'between']
df['E1'] = df['alpha_Proteins'].apply(check_amb, all_proteins=all_proteins) # if ambigous will count as ecoli
df['E2'] = df['beta_Proteins'].apply(check_amb, all_proteins=all_proteins)

df['EE'] = df['E1'] & df['E2']
df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
df['HH'] = (~df['E1']) & (~df['E2'])
print(sum(df.HH | df.EH)/len(df))

df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})
df.loc[((df['decoy1'] == True)| (df['decoy2'] == True)), 'entr_group'] = 'decoy'

score = 'E_value'#'iProbability' #'K_score'
ax = plot_distribution(df=df, x=score, bins=50)
plt.show()

ax = plot_distribution(df=df, x=score, bins=50, ylim_top=100)
plt.show()

summary_table = df.reset_index()['entr_group'].value_counts()
summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
# summary_table.to_csv('kojak.csv')
