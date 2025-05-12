import pandas as pd
from plots_and_functions import find_protein_amb


def process_xlinkx(result_file, decoy_file, proteins, all_peptide_proteins, all_peptide_positions):
    """
    Process the results of XlinkX.

    :param result_file: path to the txt/tsv file with the search results
    :param decoy_file:  path to the txt/tsv file with the decoy results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """

    # read the files
    target_df = None
    # try as tsv
    try:
        target_df = pd.read_csv(result_file, sep='\t')
        decoy_df = pd.read_csv(decoy_file, sep='\t')
    except Exception:
        pass
    if not isinstance(target_df, pd.DataFrame):
        # try as csv
        try:
            target_df = pd.read_csv(result_file)
            decoy_df = pd.read_csv(decoy_file)
        except Exception:
            pass
    if not isinstance(target_df, pd.DataFrame):
        # try as xls
        try:
            target_df = pd.read_excel(result_file)
            decoy_df = pd.read_excel(decoy_file)
        except Exception:
            pass
    # subset to heteromeric
    target_df = target_df[target_df['Crosslink Type'].str.contains('Inter')]
    decoy_df = decoy_df[decoy_df['Crosslink Type'].str.contains('Inter')]
    decoy_df = decoy_df[decoy_df['XlinkX Score'] >= target_df["XlinkX Score"].min()]

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

    df.rename(columns={"Protein Accession A": "Protein1",
                       "Protein Accession B": "Protein2",
                       "Sequence A": "peptide1",
                       "Sequence B": "peptide2",
                       "MSMS.Info": "Scannumber",
                       "Charge": "precursor_charge",
                       "Crosslinker Position A": "LinkPos1",
                       "Crosslinker Position B": "LinkPos2"
                       }, inplace=True)
    df["PepPos1"] = '-1'
    df["PepPos2"] = '-1'
    df["PepPos1"][~df["Leading Protein Position A"].isna()] = \
        df["Leading Protein Position A"] - df['LinkPos1'] + 1
    df["PepPos2"][~df["Leading Protein Position B"].isna()] = \
        df["Leading Protein Position B"] - df['LinkPos2'] + 1

    df['LinkPos1'].loc[df['LinkPos1'] == 0] = 1
    df['LinkPos2'].loc[df['LinkPos2'] == 0] = 1

    # for what ever reason only one protein is reported for every possible peptide position
    # so we stick with that and get the first protein
    df['isDecoy1'] = df["peptide1"].apply(
        lambda x: all([y.startswith("REV_") for y in all_peptide_proteins[x]]))
    df['isDecoy2'] = df["peptide2"].apply(
        lambda x: all([y.startswith("REV_") for y in all_peptide_proteins[x]]))

    # find the decoys proteins
    df["Protein1"][df["Leading Protein Position A"].astype(float).isna() & df['isDecoy1']] = \
        df["peptide1"][df["Leading Protein Position A"].astype(float).isna() & df['isDecoy1']].apply(
        lambda x: all_peptide_proteins[x][0])
    df["Protein2"][df["Leading Protein Position B"].astype(float).isna() & df['isDecoy2']] = \
        df["peptide2"][df["Leading Protein Position B"].astype(float).isna() & df['isDecoy2']].apply(
        lambda x: all_peptide_proteins[x][0])

    # find the target proteins
    df["Protein1"][df["Leading Protein Position A"].astype(float).isna() & ~df['isDecoy1']] = \
        df["peptide1"][df["Leading Protein Position A"].astype(float).isna() & ~df['isDecoy1']].apply(
        lambda x: [y for y in all_peptide_proteins[x] if not y.startswith("REV_")][0])
    df["Protein2"][df["Leading Protein Position B"].astype(float).isna() & ~df['isDecoy2']] = \
        df["peptide2"][df["Leading Protein Position B"].astype(float).isna() & ~df['isDecoy2']].apply(
        lambda x: [y for y in all_peptide_proteins[x] if not y.startswith("REV_")][0])

    # find the decoys positions
    df["PepPos1"][df["Leading Protein Position A"].isna() & df['isDecoy1']] = \
        df["peptide1"][df["Leading Protein Position A"].isna() & df['isDecoy1']].apply(
        lambda x: all_peptide_positions[x][0])
    df["PepPos2"][df["Leading Protein Position B"].isna() & df['isDecoy2']] = \
        df["peptide2"][df["Leading Protein Position B"].isna() & df['isDecoy2']].apply(
        lambda x: all_peptide_positions[x][0])

    # find the target positions
    df["PepPos1"][df["Leading Protein Position A"].isna() & ~df['isDecoy1']] = \
        df["peptide1"][df["Leading Protein Position A"].isna() & ~df['isDecoy1']].apply(
        lambda x: [peppos for prot, peppos in zip(all_peptide_proteins[x], all_peptide_positions[x])
                   if not prot.startswith("REV_")][0])
    df["PepPos2"][df["Leading Protein Position B"].isna() & ~df['isDecoy2']] = \
        df["peptide2"][df["Leading Protein Position B"].isna() & ~df['isDecoy2']].apply(
        lambda x:
        [peppos for prot, peppos in zip(all_peptide_proteins[x], all_peptide_positions[x])
         if not prot.startswith("REV_")][0])

    # df["isDecoy1"] = df["Protein1"].str.startswith("REV_")
    # df["isDecoy2"] = df["Protein2"].str.startswith("REV_")

    df.reset_index(inplace=True)
    df["PSMID"] = ["XilnkX_" + str(x) for x in df.index]

    df['isTT'] = ~(df["isDecoy1"] | df["isDecoy2"])
    df['isDD'] = (df["isDecoy1"] & df["isDecoy2"])
    df['isTD'] = ~(df["isTT"] | df["isDD"])

    # create UniqueScanID from "Spectrum File" and "CSMs"
    try:
        df["UniqueScanID"] = df["Spectrum File"].str.replace(".raw", "") + " " + df["First Scan "].astype(str)
        df['run'] = df['Spectrum File'].str.replace(".raw", "")
        df['scan'] = df['CSMs ']
    except KeyError:
        df["UniqueScanID"] = df["Spectrum File"].str.replace(".raw", "") + " " + df["First Scan"].astype(str)
        df['run'] = df['Spectrum File'].str.replace(".raw", "")
        df['scan'] = df["First Scan"]

    return df


def process_xlinkx_respair_to_xiFDR(result_file, decoy_file, proteins, all_peptide_proteins, all_peptide_positions):
    """
    Process the results of XlinkX.

    :param result_file: path to the txt/tsv file with the search results
    :param decoy_file:  path to the txt/tsv file with the decoy results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """

    # read the files
    target_df = None
    # try as tsv
    try:
        target_df = pd.read_csv(result_file, sep='\t')
        #decoy_df = pd.read_csv(decoy_file, sep='\t')
    except Exception:
        pass
    if not isinstance(target_df, pd.DataFrame):
        # try as csv
        try:
            target_df = pd.read_csv(result_file)
            #decoy_df = pd.read_csv(decoy_file)
        except Exception:
            pass
    if not isinstance(target_df, pd.DataFrame):
        # try as xls
        try:
            target_df = pd.read_excel(result_file)
            #decoy_df = pd.read_excel(decoy_file)
        except Exception:
            pass
    # subset to heteromeric
    target_df = target_df[target_df['Crosslink Type'].str.contains('Inter')]
    #decoy_df = decoy_df[decoy_df['Crosslink Type'].str.contains('Inter')]
    #decoy_df = decoy_df[decoy_df['XlinkX Score'] >= target_df["XlinkX Score"].min()]

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    target_df['E1'] = target_df['Accession A'].apply(find_protein_amb, all_proteins=proteins)
    target_df['E2'] = target_df['Accession B'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    target_df['EE'] = target_df['E1'] & target_df['E2']
    target_df['EH'] = (target_df['E1'] & (~target_df['E2']) | (~target_df['E1']) & target_df['E2'])
    target_df['HH'] = (~target_df['E1']) & (~target_df['E2'])
    # print(sum(target_df.HH | target_df.EH)/len(target_df))

    # add group column for plotting purposes
    target_df['entr_group'] = target_df['E1'].astype(int) + target_df['E2'].astype(int)
    target_df['entr_group'] = target_df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark decoys
    #decoy_df['entr_group'] = 'decoy'

    # concat target and decoy
    df = target_df #pd.concat([target_df, decoy_df])

    # add score column
    df['score'] = df['Max. XlinkX Score']

    # add search engine column
    df['search_engine'] = 'XlinkX'

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment']/summary_table['decoy']
    # summary_table.to_csv('xlinkx_numbers.csv')

    # fake some values to mkae it xiFDR compatible
    df["Scannumber"] = df.index
    df["precursor_charge"] = 1
    df["peptide1"] = df["Sequence A"].str.replace("[^A-Z]","",regex=True)
    df["peptide2"] = df["Sequence B"].str.replace("[^A-Z]","",regex=True)

    df.rename(columns={"Accession A": "Protein1",
                       "Accession B": "Protein2",
                       #"Sequence A": "peptide1",
                       #"Sequence B": "peptide2",
                       #"MSMS.Info": "Scannumber",
                       #"Charge": "precursor_charge",
                       "Position A": "LinkPos1",
                       "Position B": "LinkPos2"
                       }, inplace=True)
    df["PepPos1"] = 1
    df["PepPos2"] = 1
    #df["PepPos1"][~df["Protein Position A"].isna()] = \
    #    df["Protein Position A"] - df['LinkPos1'] + 1
    #df["PepPos2"][~df["Protein Position B"].isna()] = \
    #    df["Protein Position B"] - df['LinkPos2'] + 1

    df['LinkPos1'].loc[df['LinkPos1'] == 0] = 1
    df['LinkPos2'].loc[df['LinkPos2'] == 0] = 1

    # for what ever reason only one protein is reported for every possible peptide position
    # so we stick with that and get the first protein
    df['isDecoy1'] = df["peptide1"].apply(
        lambda x: all([y.startswith("REV_") for y in all_peptide_proteins[x]]))
    df['isDecoy2'] = df["peptide2"].apply(
        lambda x: all([y.startswith("REV_") for y in all_peptide_proteins[x]]))
    
    # column crosslink sequence is messing with xifdr
    df["Crosslink Sequence"] = df["Crosslink Sequence"].str.replace('[\n\r]','', regex=True)

    # find the decoys proteins
    #df["Protein1"][df["Position A"].isna() & df['isDecoy1']] = \
    #    df["peptide1"][df["Position A"].isna() & df['isDecoy1']].apply(
    #    lambda x: all_peptide_proteins[x][0])
    #df["Protein2"][df["Position B"].isna() & df['isDecoy2']] = \
    #    df["peptide2"][df["Position B"].isna() & df['isDecoy2']].apply(
    #    lambda x: all_peptide_proteins[x][0])

    # find the target proteins
    #df["Protein1"][df["Position A"].isna() & ~df['isDecoy1']] = \
    #    df["peptide1"][df["Position A"].isna() & ~df['isDecoy1']].apply(
    #    lambda x: [y for y in all_peptide_proteins[x] if not y.startswith("REV_")][0])
    #df["Protein2"][df["Position B"].isna() & ~df['isDecoy2']] = \
    #    df["peptide2"][df["Position B"].isna() & ~df['isDecoy2']].apply(
    #    lambda x: [y for y in all_peptide_proteins[x] if not y.startswith("REV_")][0])

    # # find the decoys positions
    # df["PepPos1"][df["Position A"].isna() & df['isDecoy1']] = \
    #     df["peptide1"][df["Position A"].isna() & df['isDecoy1']].apply(
    #     lambda x: all_peptide_positions[x][0])
    # df["PepPos2"][df["Position B"].isna() & df['isDecoy2']] = \
    #     df["peptide2"][df["Position B"].isna() & df['isDecoy2']].apply(
    #     lambda x: all_peptide_positions[x][0])

    # # find the target positions
    # df["PepPos1"][df["Position A"].isna() & ~df['isDecoy1']] = \
    #     df["peptide1"][df["Position A"].isna() & ~df['isDecoy1']].apply(
    #     lambda x: [peppos for prot, peppos in zip(all_peptide_proteins[x], all_peptide_positions[x])
    #                if not prot.startswith("REV_")][0])
    # df["PepPos2"][df["Position B"].isna() & ~df['isDecoy2']] = \
    #     df["peptide2"][df["Position B"].isna() & ~df['isDecoy2']].apply(
    #     lambda x:
    #     [peppos for prot, peppos in zip(all_peptide_proteins[x], all_peptide_positions[x])
    #      if not prot.startswith("REV_")][0])

    # df["isDecoy1"] = df["Protein1"].str.startswith("REV_")
    # df["isDecoy2"] = df["Protein2"].str.startswith("REV_")

    df.reset_index(inplace=True)
    df["PSMID"] = ["XilnkX_" + str(x) for x in df.index]

    df['isTT'] = ~(df["isDecoy1"] | df["isDecoy2"])
    df['isDD'] = (df["isDecoy1"] & df["isDecoy2"])
    df['isTD'] = ~(df["isTT"] | df["isDD"])

    # create UniqueScanID from "Spectrum File" and "CSMs"
    try:
        df["UniqueScanID"] = df["Spectrum File"].str.replace(".raw", "") + " " + df["First Scan "].astype(str)
        df['run'] = df['Spectrum File'].str.replace(".raw", "")
        df['scan'] = df['CSMs ']
    except KeyError:
        df["UniqueScanID"] = df.index
        df['run'] = "xlinkxpd32"
        df['scan'] = df.index

    return df
