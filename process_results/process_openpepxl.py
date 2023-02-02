import pandas as pd
from plots_and_functions import fasta_to_dict, find_protein_amb


def process_openpepxl(result_file, proteins):
    """
    Process the results of OpenPepXL.

    :param result_file: path to the csv file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the csv file
    df = pd.read_csv(result_file, thousands=',')

    # remove mono-links
    df = df[df['xl_type'] == 'cross-link']

    # subset to heteromeric
    df = df[df['XFDR:is_interprotein'] == True]

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    df['E1'] = df['accessions'].apply(find_protein_amb, all_proteins=proteins)
    df['E2'] = df['accessions_beta'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    df['EE'] = df['E1'] & df['E2']
    df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
    df['HH'] = (~df['E1']) & (~df['E2'])
    # print(sum(df.HH | df.EH)/len(df))

    # add group column for plotting purposes
    df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
    df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark TD and DD as decoys
    # TODO they have a decoy category 'target+decoy' which are ambigous matches --> right now I keep them as targets
    df.loc[df['target_decoy'] == 'decoy', 'entr_group'] = 'decoy'

    # add search engine column
    df['search_engine'] = 'OpenPepXL'

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment'] / summary_table['decoy']
    # summary_table.to_csv('xisearch.csv')

    return df
