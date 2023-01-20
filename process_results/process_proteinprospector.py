import pandas as pd
from plots_and_functions import *


def process_proteinprospector(result_file, proteins):
    """
    Process the results of ProteinProspector.

    :param result_file: path to the csv file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the csv file
    df = pd.read_csv(result_file)
    # '/proteinprospector/ec1_dsso_CSM_2p_inter.csv')

    # subset to heteromeric - already done in ProteinProspector ('inter' file)

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    df['E1'] = df['Acc.1'].apply(find_protein_amb, all_proteins=proteins)
    df['E2'] = df['Acc.2'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    # ToDo: is this used later or leftover?
    df['EE'] = df['E1'] & df['E2']
    df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
    df['HH'] = (~df['E1']) & (~df['E2'])
    # print(sum(df.HH | df.EH) / len(df))

    # add group column for plotting purposes
    df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
    df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark TD and DD as decoys
    df.loc[df['Decoy'] != 'Target', 'entr_group'] = 'decoy'

    # add score column
    df['score'] = df['SVM.score']

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment'] / summary_table['decoy']
    # summary_table.to_csv('proteinprospector.csv')

    return df
