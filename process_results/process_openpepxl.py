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

    # subset to unambigiously heteromeric links
    df = df[(df['XFDR:is_interprotein'] == True) & (df['XFDR:is_intraprotein'] == False)]

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
    df.reset_index(inplace=True)
    df["PSMID"] = ["OpenPepXL_" + str(x) for x in df.index]


    df.rename(columns={#"accessions": "Protein1",
               #"accessions_beta": "Protein2",
               "sequence": "peptide1",
               "sequence_beta": "peptide2",
               "charge": "precursor_charge",
               #"start": "PepPos1",
               "xl_pos1": "LinkPos1",
               "xl_pos2": "LinkPos2"
               }, inplace=True)
    df["PepPos1"] = df["start"]
    df["isDecoy1"] = df["xl_target_decoy_alpha"]=="decoy"
    df["isDecoy2"] = df["xl_target_decoy_beta"]=="decoy"
    df["LinkPos1"].loc[df["LinkPos1"]==0] = 1
    df["LinkPos2"].loc[df["LinkPos2"]==0] = 1
    # why is the beta start separated by "_" if start is separated by ";"?
    df["PepPos2"] = df["BetaPepEv:start"].str.replace("_", ";", regex=False)

    df["Protein1"] = df["accessions"]
    # turn target+decoy into target    
    mask_p1_td = df["xl_target_decoy_alpha"]=="target+decoy"
    df["Protein1"].loc[mask_p1_td] = df["Protein1"][mask_p1_td].apply(
        lambda x : ";".join([pr for pr in x.split(";") if not pr.startswith("DECOY")]))
    df["PepPos1"].loc[mask_p1_td] = df[["accessions","PepPos1"]][mask_p1_td].apply(
        lambda x : ";".join([pp for pr,pp in zip(x.iloc[0].split(";"), x.iloc[1].split(";")) if not pr.startswith("DECOY")]), axis=1) 

    #df["accessions_beta_old"] = df["accessions_beta"]
    df["accessions_beta"] = df["accessions_beta"].str.replace("(?<!DECOY)_(?=DECOY|sp)",";",regex=True)
    df["Protein2"] = df["accessions_beta"]
    mask_p2_td = df["xl_target_decoy_beta"]=="target+decoy"
    df["Protein2"].loc[mask_p2_td] = df["accessions_beta"][mask_p2_td].apply(
        lambda x : ";".join([pr for pr in x.split(";") if not pr.startswith("DECOY")]))
    df["PepPos2"].loc[mask_p2_td] = df[["accessions_beta","PepPos2"]][mask_p2_td].apply(
        lambda x : ";".join([pp for pr,pp in zip(x.iloc[0].split(";"), x.iloc[1].split(";")) if not pr.startswith("DECOY")]), axis=1)
    
    df["Protein1"] = df["Protein1"].str.replace("DECOY","REV",regex=False)
    df["Protein2"] = df["Protein2"].str.replace("DECOY","REV",regex=False)

    df['isTT'] = ~(df["isDecoy1"] | df["isDecoy2"])
    df['isDD'] = (df["isDecoy1"] & df["isDecoy2"])
    df['isTD'] = ~(df["isTT"] | df["isDD"])

    # convert file_origin and spectrum_reference into a unique spectrum id
    # file_origin:B190511_02_HF_LS_IN_130_ECLP_DSSO_01_SCX17_hSAX06_rep1.idXML
    # spectrum_reference:controllerType=0 controllerNumber=1 scan=1
    # UniqueScanID: B190511_02_HF_LS_IN_130_ECLP_DSSO_01_SCX17_hSAX06_rep1 1
    df["UniqueScanID"] = df["file_origin"].str.replace(".idXML","") + " " + df["spectrum_reference"].str.replace(".*scan=","")
    df['run'] = df["file_origin"].str.replace(".idXML","")
    df['scan'] = df["spectrum_reference"].str.replace(".*scan=","", regex=True)

    return df
