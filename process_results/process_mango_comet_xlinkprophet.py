import pandas as pd
from plots_and_functions import find_protein_amb
import re


def find_in_sequence(protein, peptide):
    """
    Find the (first) position of a peptide in a protein.
    If peptide can't be found directly a I/L equivalency is assumed,

    :param protein: The protein that is supposed to be the source of the peptide
    :param peptide: the peptide sequence to look for

    :returns: position of the peptide in the protein
    """
    pos = str(protein).find(peptide)
    if pos < 0:
        pos = str(protein).replace("I", "L").index(peptide.replace("I", "L"))
    return str(pos)


def process_mango_comet_xlinkprophet(result_file, proteins, td_fasta_dict):
    """
    Process the results of MangoCometXlinkProphet.

    :param result_file: path to the xlsx file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the csv file
    try:
        df = pd.read_excel(result_file)
    except Exception as e:
        print(f"Opening file {result_file} failed with error {e}." +
              " Trying to open as tab separated file")
        df = pd.read_csv(result_file, sep="\t")


    # find overlap between protein1 and protein2 columns
    df['myintra'] = df.apply(
        lambda x:
        len(set(x['protein1'].replace("rev_", "").split(','))
            .intersection(set(x['protein2']
                              .replace("rev_", "").split(',')))), axis=1)
    # transform spectrum into UniqueScanID
    usid_p = re.compile("([^.]*)\.([0-9]*)\..*")
    df["UniqueScanID"] = df["spectrum"].apply(lambda x: usid_p.sub("\\1 \\2", x))
    df["run"] = df["spectrum"].apply(lambda x: usid_p.sub("\\1", x))
    df["scan"] = df["spectrum"].apply(lambda x: int(usid_p.sub("\\2", x)))

    # subset to heteromeric
    df = df[df.myintra == 0]

    # are protein 1 or protein 2 from E. coli --> gives true / false in
    # separate column if ambiguous match contains an E. coli protein
    # return True
    df['E1'] = df['protein1'].apply(find_protein_amb, all_proteins=proteins)
    df['E2'] = df['protein2'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    df['EE'] = df['E1'] & df['E2']
    df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
    df['HH'] = (~df['E1']) & (~df['E2'])
    # print(sum(df.HH | df.EH) / len(df))

    # add group column for plotting purposes
    df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
    df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark TD and DD as decoys
    df.loc[((df['decoy'] == 1) | (df['decoy-decoy'] == 1)), 'entr_group'] = 'decoy'
    df['isTT'] = df['decoy'] == 0
    df['isTD'] = (df['decoy'] - df['decoy-decoy']) == 1
    df['isDD'] = df['decoy-decoy'] == 1
    # add score column
    df['score'] = df['probability']

    # add search engine column
    df['search_engine'] = 'MangoCometXlinkProphet'

    df.reset_index(inplace=True)
    df["PSMID"] = ["MangoCometXlinkProphet_" + str(x) for x in df.index]

    # BIG FAT TODO
    df.rename(columns={
               "peptide1_len": "peplen1",
               "peptide2_len": "peplen2",
               "beta_Pep_Link": "LinkPos2",
               "parent_charge": "precursor_charge"}, inplace=True)
    df["isDecoy1"] = df["protein1"].apply(lambda x: all([y.startswith("rev_")
                                                         for y in x.split(",")]))
    df["isDecoy2"] = df["protein2"].apply(lambda x: all([y.startswith("rev_")
                                                         for y in x.split(",")]))
    # convert into semicolon separated list (xiFDR likes that more)
    df["protein1"] = df["protein1"].str.replace(',', ';')
    df["protein2"] = df["protein2"].str.replace(',', ';')
    df["protein1"] = df["protein1"].str.replace('sp|', '', regex=False)
    df["protein2"] = df["protein2"].str.replace('sp|', '', regex=False)
    df["protein1"] = df["protein1"].str.replace(r'\|[^;]*', '', regex=True)
    df["protein2"] = df["protein2"].str.replace(r'\|[^;]*', '', regex=True)
    df["protein1"] = df["protein1"].str.replace('rev_', 'REV_')
    df["protein2"] = df["protein2"].str.replace('rev_', 'REV_')

    compid = df['composite_id'].apply(lambda x: re.sub("\[^\]]*\]", "", x))
    df["LinkPos1"] = compid.apply(lambda x: re.sub("_.*", "", x).index("x"))
    df["LinkPos2"] = compid.apply(lambda x: re.sub(",*_", "", x).index("x"))

    df["basepeptide1"] = df["peptide1"].str.replace("[^A-Z]*", "", regex=True)
    df["basepeptide2"] = df["peptide2"].str.replace("[^A-Z]*", "", regex=True)
    # find peptide positions in proteins
    df["PepPos1"] = df[["basepeptide1", "protein1"]].apply(
        lambda x: ";".join([find_in_sequence(td_fasta_dict[y].seq, x.iloc[0])
                            for y in x.iloc[1].split(";")]), axis=1)
    df["PepPos2"] = df[["basepeptide2", "protein2"]].apply(
        lambda x: ";".join([find_in_sequence(td_fasta_dict[y].seq, x.iloc[0])
                            for y in x.iloc[1].split(";")]), axis=1)

    return df


def process_mango_comet_xlinkprophet_ppi(result_file, proteins, td_fasta_dict):
    """
    Process the results of MangoCometXlinkProphet.

    :param result_file: path to the xlsx file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """
    # read the csv file
    df = pd.read_csv(result_file, sep="\t")

    # split protein_pair column by - into Protein1 and Protein2
    df[['Protein1', 'Protein2']] = df['protein_pair'].str.split('-', expand=True)

    # turn proteins into lists
    df['Protein1_nondec'] = df['Protein1'].str.replace('rev_', '').str.split(',')
    df['Protein2_nondec'] = df['Protein2'].str.replace('rev_', '').str.split(',')
    df['Protein1_list'] = df['Protein1'].str.split(',')
    df['Protein2_list'] = df['Protein2'].str.split(',')

    # find overlap between Protein1 and Protein2
    df['intra'] = df.apply(lambda x:
                           len(set(x['Protein1_nondec'])
                               .intersection(set(x['Protein2_nondec']))), axis=1)
    df['fdrGroup'] = df['intra'].apply(lambda x: 'Self' if x > 0 else 'Between')
    df['decoy1'] = df.apply(lambda x: all(p.startswith('rev_') for p in x['Protein1_list']), axis=1)
    df['decoy2'] = df.apply(lambda x: all(p.startswith('rev_') for p in x['Protein2_list']), axis=1)

    # subset to heteromeric
    df = df[df.intra == 0]

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
    df['isTT'] = ~(df['decoy1'] | df['decoy2'])
    df['isDD'] = df['decoy1'] & df['decoy2']
    df['isTD'] = df['decoy1'] ^ df['decoy2']
    df.loc[~df['isTD'], 'entr_group'] = 'decoy'
    # add score column
    df['score'] = df['probability']

    # add search engine column
    df['search_engine'] = 'MangoCometXlinkProphet'

    df.reset_index(inplace=True)
    df["PSMID"] = ["MangoCometXlinkProphet_" + str(x) for x in df.index]

    # BIG FAT TODO
    df.rename(columns={
               "peptide1_len": "peplen1",
               "peptide2_len": "peplen2",
               "beta_Pep_Link": "LinkPos2",
               "parent_charge": "precursor_charge"}, inplace=True)
    df["isDecoy1"] = df["Protein1"].apply(lambda x: all([y.startswith("rev_")
                                                         for y in x.split(",")]))
    df["isDecoy2"] = df["Protein2"].apply(lambda x: all([y.startswith("rev_")
                                                         for y in x.split(",")]))
    # convert into semicolon separated list (xiFDR likes that more)
    df["Protein1"] = df["Protein1"].str.replace(',', ';')
    df["Protein2"] = df["Protein2"].str.replace(',', ';')
    df["Protein1"] = df["Protein1"].str.replace('sp|', '', regex=False)
    df["Protein2"] = df["Protein2"].str.replace('sp|', '', regex=False)
    df["Protein1"] = df["Protein1"].str.replace(r'\|[^;]*', '', regex=True)
    df["Protein2"] = df["Protein2"].str.replace(r'\|[^;]*', '', regex=True)
    df["Protein1"] = df["Protein1"].str.replace('rev_', 'REV_')
    df["Protein2"] = df["Protein2"].str.replace('rev_', 'REV_')

    return df
