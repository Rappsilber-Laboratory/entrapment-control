import pandas as pd
from plots_and_functions import *


def check_amb_fdr_group(x):
    """
    ToDo: write docstring

    :param x:
    :return:
    """
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


def process_kojak(result_file, proteins):
    """
    Process the results of Kojak.

    :param result_file: path to the tsv file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the csv file
    df = pd.read_csv(result_file, sep='\t')

    df['decoy1'] = df['alpha_Proteins'].str.contains('DECOY')
    df['decoy2'] = df['beta_Proteins'].str.contains('DECOY')  # if ambiguous will count as decoy

    # subset to heteromeric
    df['fdrGroup'] = df.apply(check_amb_fdr_group, axis=1)
    df = df[df['fdrGroup'] == 'between']

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    df['E1'] = df['alpha_Proteins'].apply(find_protein_amb, all_proteins=proteins)
    df['E2'] = df['beta_Proteins'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    df['EE'] = df['E1'] & df['E2']
    df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
    df['HH'] = (~df['E1']) & (~df['E2'])
    # print(sum(df.HH | df.EH) / len(df))

    # add group column for plotting purposes
    df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
    df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark TD and DD as decoys
    df.loc[(df['decoy1'] | df['decoy2']), 'entr_group'] = 'decoy'

    # add score column
    df['score'] = df['E_value']

    # add search engine column
    df['search_engine'] = 'Kojak'

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment'] / summary_table['decoy']
    # # summary_table.to_csv('kojak.csv')

    return df
