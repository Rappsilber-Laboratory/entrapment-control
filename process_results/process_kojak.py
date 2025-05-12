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
        proteins_1 = prot1.replace("REV_", "").split(';')
        for prot1_i in proteins_1:
            if ';' in prot2:
                proteins_2 = prot2.replace("REV_", "").split(';')
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
        proteins_2 = prot2.replace("REV_", "").split(';')
        for prot2_i in proteins_2:
            if ';' in prot1:
                proteins_1 = prot1.replace("REV_", "").split(';')
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


def process_kojak(result_file, proteins, fasta):
    """
    Process the results of Kojak.

    :param result_file: path to the tsv file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """

    # read in the fasta file
    fasta_proteins = SeqIO.parse(fasta, "fasta")
    
    prot_array = []
    decoy_to_target = {}
    target_prot_id_to_prot = {}
    target_prot_id = 1
    # match decoy to target
    for prot_id, protein in enumerate(fasta_proteins):

        prot_array.append(protein.id)
        if 'decoy' in protein.id.lower() or 'reverse' in protein.id.lower():
            decoy_to_target[protein.id] = prot_array[prot_id - 1].split('|')[1]
        else:
            target_prot_id_to_prot[target_prot_id] = protein.id.split('|')[1]
            target_prot_id += 1



    # read the csv file
    df = pd.read_csv(result_file, sep='\t')

    # actually result file does not contain any overlap of target and decoy for a single peptide
    # so no need to check for ambiguous matches when looking for decoys
    df['decoy1'] = df['alpha_Proteins'].str.contains('DECOY')
    df['decoy2'] = df['beta_Proteins'].str.contains('DECOY')  


    # convert to matchable proteins
    df['alpha_Proteins'] = df['alpha_Proteins'].apply(lambda x: ";".join(
        ["REV_" + target_prot_id_to_prot[int(y[y.rfind('_') + 1:])] 
         if 'DECOY' in  y 
         else y.split('|')[1] 
         for y in x.split(";")]))
    df['beta_Proteins'] = df['beta_Proteins'].apply(lambda x: ";".join(
        ["REV_" + target_prot_id_to_prot[int(y[y.rfind('_') + 1:])] 
         if 'DECOY' in  y 
         else y.split('|')[1] 
         for y in x.split(";")]))

    # subset to heteromeric
    df['fdrGroup'] = df.apply(check_amb_fdr_group, axis=1)
    df = df[df['fdrGroup'] == 'between']
    # merge Sample and Scan as UniqueScanID - so we can compare overlap on scan based level
    df['UniqueScanID'] = df['Sample'] + ' ' + df['Scan'].astype(str)
    df['run'] = df['Sample']
    df['scan'] = df['Scan']
    

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
    df['Score'] = df['Sum_score']

    # add search engine column
    df['search_engine'] = 'Kojak'
    df['isTT'] = ~(df['decoy1'] | df['decoy2'])
    df['isDD'] = df['decoy1'] & df['decoy2']
    df['isTD'] = ~(df['isTT'] | df['isDD'])

    # summary_table = df.reset_index()['entr_group'].value_counts()
    # summary_table['ratio_entrapment_decoy'] = summary_table['entrapment'] / summary_table['decoy']
    # # summary_table.to_csv('kojak.csv')

    df.reset_index(inplace=True)
    df["PSMID"] = ["Kojak_" + str(x) for x in df.index]

    # convert "alpha_Pep_Link" and "alpha_Prot_Link" into peptide position
    df["PepPos1"]=df[["alpha_Pep_Link","alpha_Prot_Link"]].apply(
        lambda x : ";".join([str(int(i) -x.iloc[0]+1) for i in x.iloc[1].split(";")]), axis=1 ) 
    df["PepPos2"]=df[["beta_Pep_Link","beta_Prot_Link"]].apply(
        lambda x : ";".join([str(int(i) -x.iloc[0]+1) for i in x.iloc[1].split(";")]), axis=1 ) 
    # rename some columns - so they can be read directly by xiFDR
    df.rename(columns={"alpha_Proteins": "Protein1",
               "beta_Proteins": "Protein2",
               "alpha_Peptide": "peptide1",
               "beta_Peptide": "peptide2",
               "alpha_Pep_Link": "LinkPos1",
               "beta_Pep_Link": "LinkPos2",
               "Z": "precursor_charge"}, inplace=True)
    

    return df
