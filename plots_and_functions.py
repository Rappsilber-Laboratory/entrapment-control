from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt


def fasta_to_dict(fasta_file, exclude=None):
    """
    Read in a fasta file and return a dictionary.

    :param fasta_file: path to the fasta file
    :param exclude: list of fasta headers to exclude
    :return: dictionary with fasta headers as keys and sequences as values
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    record_dict_old = record_dict.copy()

    for key, value in record_dict_old.items():
        if 'Decoy' in key:
            new = 'DECOY_' + key.split('|')[1]
        else:
            new = key.split('|')[1]
        if (exclude is not None) and (new in exclude):
            print('excluded')
            del record_dict[key]
        else:
            record_dict[new] = value
            del record_dict[key]
    return record_dict


def find_protein_amb(protein, all_proteins):
    """
    Check if a protein is in the list of all proteins (considering ambiguous matches).

    :param protein: protein to search for
    :param all_proteins: list of all proteins
    :return: True if protein is in list, False otherwise
    """
    if ';' in protein:
        proteins = protein.split(';')
        for prot_i in proteins:
            if '|' in prot_i:
                prot_i = prot_i.split('|')[1]
            if prot_i in all_proteins:
                return True
            return False
    else:
        if '|' in protein:
            protein = protein.split('|')[1]
        if protein in all_proteins:
            return True
        else:
            return False


def plot_distribution(df, x, bins, ylim_top=None):
    palette = sns.color_palette(['#8CCFC0', '#83320C', '#E364B4'])
    df['ID_type'] = df['entr_group']
    fig, ax = plt.subplots()
    sns.histplot(data=df, x=x, hue='ID_type', element="step", bins=bins, linewidth=1.7,
                 # hist_kws={"linewidth": 5, "alpha": 1},
                 hue_order=['E.coli', 'decoy', 'entrapment'], palette=palette)
    if ylim_top:
        plt.ylim(top=ylim_top)
    # plt.show()
    return ax
