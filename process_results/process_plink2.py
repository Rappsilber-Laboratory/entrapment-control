import pandas as pd
import re
from plots_and_functions import *
from process_results.convert_plink_to_xiFDR import convert_df


def process_plink2(result_file, unfiltered_result_file, proteins, all_peptide_proteins,
                   all_peptide_positions):
    """
    Process the results of pLink2.

    :param result_file: path to the csv file with the search results
    :param unfiltered_result_file: path to the csv file with the unfiltered search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # df = pd.read_csv(result_file)
    df = convert_df(result_file)
    # '/plink2/DSSO_HumanDB_plinkFDR_xifdr.csv')
    # subset to heteromeric
    df = df[(df['Group'] == 'between')]
    # reverse the score
    df['1 - Score'] = df['Score'].apply(lambda x: 1. - x)

    # read the unfiltered csv file
    # df_with_dec = pd.read_csv(unfiltered_result_file)
    df_with_dec = convert_df(unfiltered_result_file)
    # 'plink2/DSSO_HumanDB_unfiltered_xifdr.csv')
    # subset to heteromeric
    df_with_dec = df_with_dec[(df_with_dec['Group'] == 'between')]
    # reverse the score
    df_with_dec['1 - Score'] = df_with_dec['Score'].apply(lambda x: 1. - x)
    df_dec = df_with_dec[(~df_with_dec['isTT']) & (df_with_dec['1 - Score'] >= df['1 - Score'].min())]
    # concat the two dataframes
    df = pd.concat([df, df_dec])
    
    df.loc[~df["PepPos1"].isna(),"isDecoy1"] = False
    df.loc[~df["PepPos2"].isna(),"isDecoy2"] = False
    df.loc[df["PepPos1"].isna(),"isDecoy1"] = \
        df.loc[df["PepPos1"].isna(),"Peptide1"].apply(
        lambda x: all([i.startswith("REV_") for i in all_peptide_proteins[x]]))
    df.loc[df["PepPos2"].isna(),"isDecoy2"] = \
        df.loc[df["PepPos2"].isna(),"Peptide2"].apply(
        lambda x: all([i.startswith("REV_") for i in all_peptide_proteins[x]]))


    df.loc[df["PepPos1"].isna(),"Protein1"] = \
        df.loc[df["PepPos1"].isna(),"Peptide1"].apply(lambda x: ";".join(all_peptide_proteins[x]))
    df.loc[df["PepPos1"].isna(),"PepPos1"] = \
        df.loc[df["PepPos1"].isna(),"Peptide1"].apply(lambda x: ";".join(all_peptide_positions[x]))
    df.loc[df["PepPos2"].isna(),"Protein2"] = \
        df.loc[df["PepPos2"].isna(),"Peptide2"].apply(lambda x: ";".join(all_peptide_proteins[x]))
    df.loc[df["PepPos2"].isna(),"PepPos2"] = \
        df.loc[df["PepPos2"].isna(),"Peptide2"].apply(lambda x: ";".join(all_peptide_positions[x]))


    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    # TODO will throw error right now because decoys are nans - will get overwritten later anyway
    df.loc[df['Protein1'].isnull(), 'Protein1'] = 'decoy'
    df.loc[df['Protein2'].isnull(), 'Protein2'] = 'decoy'
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

    df.reset_index(inplace=True)
    df["PSMID"] = ["pLink2_" + str(x) for x in df.index]

    df['isTT'] = ~(df["isDecoy1"] | df["isDecoy2"])
    df['isDD'] = (df["isDecoy1"] & df["isDecoy2"])
    df['isTD'] = ~(df["isTT"] | df["isDD"])
    # transform title into UniqueScanID
    usid_p = re.compile("([^.]*)\.([0-9]*)\.([0-9]*)?.*")
    df["UniqueScanID"] = df["Title"].apply(lambda x: usid_p.sub("\\1 \\2", x))
    df["run"] = df["Title"].apply(lambda x: usid_p.sub("\\1", x))
    df["scan"] = df["Title"].apply(lambda x: usid_p.sub("\\2", x))

    return df


def process_plink2_ppi(result_file, proteins):
    """
    Process the pLink2 PPI output file and generate a dataframe with heteromeric interactions.
    :param result_file: pLink2 filtered_sites output file
    :param proteins: list of e-coli proteins

    :return: DataFrame with the following columns:
        - Protein1: protein 1
        - Protein2: protein 2
        - GroupID: group id
        - E1: is protein 1 from E. coli
        - E2: is protein 2 from E. coli
        - EE: is protein 1 and protein 2 from E. coli
        - EH: is protein 1 and protein 2 from different species
        - HH: is protein 1 and protein 2 from Human
        summary dictionary with number of EE, EH, HH and Estimated FP (Total number of known Entrapment + Estimated FP among EE)
    """
    df = pd.read_csv(result_file,header=1)
    df_ppi_header = pd.read_csv(result_file)
    # remove residuepairs
    df_ppi = df[df['Spectrum_Order'].str.isnumeric() == False]

    # rename the columns
    columns = list(df_ppi.columns)
    for c in range(len(df_ppi_header.columns)):
        columns[c] = df_ppi_header.columns[c]
    df_ppi.columns = columns

    # make sure all the same and subset protein pairs get a common group id    
    df_ppi['GroupID'] = df_ppi['Protein_Order'].str.isnumeric().cumsum()

    # Protein1 is everything before the opening bracket
    df_ppi['Protein1'] = df_ppi["Protein"].str.replace(r'\(.*', '')  
    # Protein2 is everything after - and before the next closing bracket
    df_ppi['Protein2'] = df_ppi["Protein"].str.replace(r'.*\-', '').str.replace(r'\(.*', '')

    # swap protein1 and protein2 if protein1 is alphabetically after protein2
    df_ppi.loc[df_ppi['Protein1'] > df_ppi['Protein2'],['Protein1', 'Protein2']] = df_ppi.loc[df_ppi['Protein1'] > df_ppi['Protein2'],['Protein2', 'Protein1']]

    # find protein type homomeric pairs
    df_ppi['Homomeric'] = df_ppi['Protein1'] == df_ppi['Protein2']

    # remove all group ids that have homomeric pairs
    df_ppi_hetero = df_ppi[~df_ppi['GroupID'].isin(df_ppi[df_ppi['Homomeric'] == True]['GroupID'])]

    # mark ecoli proteins
    df_ppi_hetero['E1'] = df_ppi_hetero['Protein1'].apply(find_protein_amb, all_proteins=proteins)
    df_ppi_hetero['E2'] = df_ppi_hetero['Protein2'].apply(find_protein_amb, all_proteins=proteins)


    # join Protein1 and Protein1 for each group by ';'
    df_ppi_hetero = df_ppi_hetero.groupby('GroupID').agg({'Protein1': lambda x: ';'.join(x), 'Protein2': lambda x: ';'.join(x), 
                                                          'E1': lambda x: any(x), 'E2': lambda x: any(x)}).reset_index()

    df_ppi_hetero = df_ppi_hetero.groupby(['Protein1', 'Protein2', 'E1', 'E2']).count().reset_index()

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    df_ppi_hetero['EE'] = df_ppi_hetero['E1'] & df_ppi_hetero['E2']
    df_ppi_hetero['EH'] = (df_ppi_hetero['E1'] & (~df_ppi_hetero['E2']) | (~df_ppi_hetero['E1']) & df_ppi_hetero['E2'])
    df_ppi_hetero['HH'] = (~df_ppi_hetero['E1']) & (~df_ppi_hetero['E2'])

    # summaries the dataframe
    summary = df_ppi_hetero[['EE', 'EH', 'HH']].sum(axis=0).to_dict()
    # absolute FP = known FP + estimated FP among EE
    # absolute FP = EH + HH  + EH - HH = 2 * EH
    summary["Estimated FP"] = 2 * summary["EH"]

    return df_ppi_hetero, summary
