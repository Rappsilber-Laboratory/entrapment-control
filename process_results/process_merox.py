import pandas as pd
from plots_and_functions import find_protein_amb


def check_amb_fdr_group(x):
    """
    ToDo: write docstring

    :param x:
    :return:
    """
    prot1, prot2 = x['Protein 1'], x['Protein 2']
    if ';' in prot1:
        proteins_1 = prot1.split(';')
        for prot1_i in proteins_1:
            if ';' in prot2:
                proteins_2 = prot2.split(';')
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
        proteins_2 = prot2.split(';')
        for prot2_i in proteins_2:
            if ';' in prot1:
                proteins_1 = prot1.split(';')
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


def process_merox(result_file, proteins):
    """
    Process the results of MeroX.

    :param result_file: path to the csv file with the search results
    :param proteins: list of all proteins in the database
    :return: DataFrame with the results
    """

    # read the csv file
    df = pd.read_csv(result_file, sep=';')

    # subset to heteromeric
    df['fdrGroup'] = df.apply(check_amb_fdr_group, axis=1)
    df = df[df['fdrGroup'] == 'between']

    # are protein 1 or protein 2 from E. coli --> gives true / false in separate column
    # if ambiguous match contains an E. coli protein return True
    df['E1'] = df['Protein 1'].apply(find_protein_amb, all_proteins=proteins)
    df['E2'] = df['Protein 2'].apply(find_protein_amb, all_proteins=proteins)

    # add columns differentiating E. coli - E.coli (EE), E. coli - human (EH), human - human (HH)
    df['EE'] = df['E1'] & df['E2']
    df['EH'] = (df['E1'] & (~df['E2']) | (~df['E1']) & df['E2'])
    df['HH'] = (~df['E1']) & (~df['E2'])
    # print(sum(df.HH | df.EH)/len(df))

    # add group column for plotting purposes
    df['entr_group'] = df['E1'].astype(int) + df['E2'].astype(int)
    df['entr_group'] = df['entr_group'].replace({2: 'E.coli', 1: 'entrapment', 0: 'entrapment'})

    # mark TD and DD as decoys
    df['isDecoy1'] = df['Protein 1'].str.startswith('DEC_')
    df['isDecoy2'] = df['Protein 2'].str.startswith('DEC_')
    df['isDD'] = ((df['isDecoy1'] == True) & (df['isDecoy2'] == True))
    df['isTD'] = ((df['isDecoy1'] == True) | (df['isDecoy2'] == True)) \
                 & (~((df['isDecoy1'] == True) & (df['isDecoy2'] == True)))
    df['isTT'] = ((df['isDecoy1'] == False) & (df['isDecoy2'] == False))
    df.loc[((df['isTD']) | (df['isDD'])), 'entr_group'] = 'decoy'

    # add score column
    df['score'] = df['Score']

    # add search engine column
    df['search_engine'] = 'MeroX'
    df.reset_index(inplace=True)
    df["PSMID"] = ["MeroX_" + str(x) for x in df.index]

    # convert protein 1 and protein2 into a format that xiFDR can handle
    df['name1'] = df['Protein 1'].str.replace(r'(?:DEC_)?>sp\|(?:[^|]*\|([^\s]*).*)',"\\1",regex=True)
    df['name2'] = df['Protein 2'].str.replace(r'(?:DEC_)?>sp\|(?:[^|]*\|([^\s]*).*)',"\\1",regex=True)
    df['description1'] = df['Protein 1'].str.replace(r'(?:DEC_)?>sp\|(?:[^|]*\|(?:[^\s]*)\s(.*))',"\\1",regex=True)
    df['description2'] = df['Protein 2'].str.replace(r'(?:DEC_)?>sp\|(?:[^|]*\|(?:[^\s]*)\s(.*))',"\\1",regex=True)
    # get reduce the name
    df['Protein1'] = df['Protein 1'].str.replace("(?:DEC_)?>sp\|([^\|]*).*","\\1",regex=True)
    df['Protein2'] = df['Protein 2'].str.replace("(?:DEC_)?>sp\|([^\|]*).*","\\1",regex=True)
    df['Protein1'].loc[df['isDecoy1']]='REV_' + df['Protein1'][df['isDecoy1']]
    df['Protein2'].loc[df['isDecoy2']]='REV_' + df['Protein2'][df['isDecoy2']]
    df['LinkPos1'] = df['best linkage position peptide 1'].str.slice(1)
    df['LinkPos2'] = df['best linkage position peptide 2'].str.slice(1)
    df['PepPos1'] = df["From"]
    df['PepPos2'] = df["From.1"]
    df['PepPos1'].loc[df['PepPos1']==0] = 1
    df['PepPos2'].loc[df['PepPos2']==0] = 1

    return df
