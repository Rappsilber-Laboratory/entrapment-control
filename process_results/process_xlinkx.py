import pandas as pd
from plots_and_functions import *


def process_xlinkx(result_file, decoy_file, proteins):
    """
    Process the results of XlinkX.

    :param result_file: path to the txt/tsv file with the search results
    :param decoy_file:  path to the txt/tsv file with the decoy results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the files
    target_df = pd.read_csv(result_file, sep='\t')
    decoy_df = pd.read_csv(decoy_file, sep='\t')
    # subset to heteromeric
    target_df = target_df[target_df['Crosslink Type'].str.contains('Inter')]
    decoy_df = decoy_df[decoy_df['Crosslink Type'].str.contains('Inter')]

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    target_df['E1'] = target_df['Protein Accession A'].apply(find_protein_amb, all_proteins=proteins)
    target_df['E2'] = target_df['Protein Accession B'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    target_df['EE'] = target_df['E1'] & target_df['E2']
    target_df['EH'] = (target_df['E1'] & (~target_df['E2']) | (~target_df['E1']) & target_df['E2'])
    target_df['HH'] = (~target_df['E1']) & (~target_df['E2'])
    # print(sum(target_df.HH | target_df.EH)/len(target_df))

    # add group column for plotting purposes
    target_df['entr_group'] = target_df['E1'].astype(int) + target_df['E2'].astype(int)
    target_df['entr_group'] = target_df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark decoys
    decoy_df['entr_group'] = 'decoy'

    # concat target and decoy
    df = pd.concat([target_df, decoy_df])

    # add score column
    df['score'] = df['XlinkX Score']

    # add search engine column
    df['search_engine'] = 'XlinkX'

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
    # summary_table.to_csv('xlinkx_numbers.csv')

    return df
