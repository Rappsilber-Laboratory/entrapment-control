import pandas as pd
from plots_and_functions import *


def process_plink2(result_file, unfiltered_result_file, proteins):
    """
    Process the results of pLink2.

    :param result_file: path to the csv file with the search results
    :param unfiltered_result_file: path to the csv file with the unfiltered search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the csv file
    df = pd.read_csv(result_file)
    # '/plink2/DSSO_HumanDB_plinkFDR_xifdr.csv')
    # subset to heteromeric
    df = df[(df['Group'] == 'between')]
    # reverse the score
    df['1 - Score'] = df['Score'].apply(lambda x: 1. - x)

    # read the unfiltered csv file
    df_with_dec = pd.read_csv(unfiltered_result_file)
    # 'plink2/DSSO_HumanDB_unfiltered_xifdr.csv')
    # subset to heteromeric
    df_with_dec = df_with_dec[(df_with_dec['Group'] == 'between')]
    # reverse the score
    df_with_dec['1 - Score'] = df_with_dec['Score'].apply(lambda x: 1. - x)

    # concat the two dataframes
    df = pd.concat([df, df_with_dec[
        (~df_with_dec['isTT']) & (df_with_dec['1 - Score'] >= df['1 - Score'].min())]])

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    df['E1'] = df['Protein1'].apply(find_protein_amb, all_proteins=proteins)
    df['E2'] = df['Protein2'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    df['EE'] = df['E1'] & df['E2']
    df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
    df['HH'] = (~df['E1']) & (~df['E2'])
    # print(sum(df.HH | df.EH) / len(df))

    # add group column for plotting purposes
    df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
    df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark TD and DD as decoys
    df.loc[((df['isTD']) | (df['isDD'])), 'entr_group'] = 'decoy'

    # add score column
    df['score'] = df['1 - Score']

    # add search engine column
    df['search_engine'] = 'pLink2'

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment'] / summary_table['decoy']
    # summary_table.to_csv('plink2.csv')

    return df
