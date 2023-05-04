import pandas as pd
from plots_and_functions import *


def process_proteinprospector(result_file, proteins, all_protein_dict):
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

    # add search engine column
    df['search_engine'] = 'ProteinProspector'

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment'] / summary_table['decoy']
    # summary_table.to_csv('proteinprospector.csv')
    df.reset_index(inplace=True)
    df["PSMID"] = ["ProteinProspector" + str(x) for x in df.index]

    #df["PepCounts1"] = df["DB.Peptide.1"].apply(lambda x: len(all_peptide_proteins[x]) if x in all_peptide_proteins else 0)
    #df["PepCounts2"] = df["DB.Peptide.2"].apply(lambda x: len(all_peptide_proteins[x]) if x in all_peptide_proteins else 0)
    #df["PepProts1"] = df["DB.Peptide.1"].apply(lambda x: ";".join(all_peptide_proteins[x]) if x in all_peptide_proteins else "")
    #df["PepProts2"] = df["DB.Peptide.2"].apply(lambda x: ";".join(all_peptide_proteins[x]) if x in all_peptide_proteins else "")
    #df["PepPos1"] = df["DB.Peptide.1"].apply(lambda x: ";".join(all_peptide_positions[x]) if x in all_peptide_positions else "")
    #df["PepPos2"] = df["DB.Peptide.2"].apply(lambda x: ";".join(all_peptide_positions[x]) if x in all_peptide_positions else "")
    #df["PepPos2"] = df["DB.Peptide.2"].apply(lambda x: ";".join(all_peptide_positions[x]))

    df["isDecoy1"] = False
    df["isDecoy2"] = False
    df["isDecoy1"][df["Acc.1"] == "decoy"] = True
    df["isDecoy2"][df["Acc.2"] == "decoy"] = True
    df["PepPos1"] = 1
    df["PepPos2"] = 1
    # rename some columns - so they can be read directly by xiFDR
    # some points 
    # - as decoys accesion numbers are not resolved I use the protein description as accession
    # - I have some trouble finding the link-position in the peptide and the position of the peptide in proteins 
    #   seems like peptide positions are not all reported but as a hack I use the linked protein residues as 
    #   peptide link position and put the peptide as position 1
    #   this way xiFDR can merge things up "correctly" for residue pairs
    #   "correctly" in quotes, as ambiguous peptides seem to be associated to a single position in a single protein
    df.rename(columns={"Acc.1": "Protein1",
               "Acc.2": "Protein2",
               "Protein.1": "description1",
               "Protein.2": "description2",
               "DB.Peptide.1": "peptide1",
               "DB.Peptide.2": "peptide2",
               "Len.Pep.1": "PepLen1",
               "Len.Pep.2": "PepLen2",
               "MSMS.Info": "Scannumber",
               "z": "precursor_charge",
               "XLink.AA.1": "LinkPos1",
               "XLink.AA.2": "LinkPos2"
               }, inplace=True)
    # I can take the protein description of decoys to generate the "correct" accession number
    df["Protein1"].loc[df["isDecoy1"]] = df["description1"][df["isDecoy1"]].apply(lambda x: ["REV_" + re.sub(r'sp\|(.*)\|.*', r'\1', p.name) for p in all_protein_dict.values() if p.description.find(" " + x + " OS=")>0][0])
    df["Protein2"].loc[df["isDecoy2"]] = df["description2"][df["isDecoy2"]].apply(lambda x: ["REV_" + re.sub(r'sp\|(.*)\|.*', r'\1', p.name) for p in all_protein_dict.values() if p.description.find(" " + x + " OS=")>0][0])
    df['isTT'] = ~(df["isDecoy1"] | df["isDecoy2"])
    df['isDD'] = (df["isDecoy1"] & df["isDecoy2"])
    df['isTD'] = ~(df["isTT"] | df["isDD"])
    
    return df
