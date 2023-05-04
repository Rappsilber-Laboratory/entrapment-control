from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import re


def reverse_fasta(fasta_dict, KRSwap=False):
    """
    Read in a fasta file and return a dictionary.

    :param fasta_file: path to the fasta file
    :param exclude: list of fasta headers to exclude
    :return: dictionary with fasta headers as keys and sequences as values
    """
    record_dict = {}
    # Loop through the fasta dictionary
    for key, value in fasta_dict.items():
        # Reverse complement the sequence
        reversed = value.reverse_complement()
        # Optionally swap K and R residues
        if KRSwap:
            reversed = reversed.replace('K', 'X').replace('R', 'K').replace('X', 'R')
        # Add the reversed sequence to the dictionary
        record_dict['REV_' + key] = reversed
    return record_dict

def find_peptide_positions(peptide, fasta_dict, tryptic_only=True):
    """
    Find the positions of a peptide in a fasta file.
    Terribly slow, so better use digest_proteins_tryptic look up in the returned dictionaries.

    :param fasta_dict: dictionary with fasta headers as keys and sequences as values
    :param peptide: peptide to search for
    :return: dictionary with fasta headers as keys and positions as values
    """
    pep_prot = []
    pep_pos = []
    for key, value in fasta_dict.items():
        if tryptic_only:
            positions = [str(m.start()) for m in re.finditer('((?<=[KR])|(?<=^))('+peptide+')', str(value.seq))]
        else:
            positions = [str(m.start()) for m in re.finditer('(?='+peptide+')', str(value.seq))]
        pep_pos.extend(positions)
        pep_prot.extend([key] * len(positions))    

    return ";".join(pep_prot), ";".join(pep_pos)

def find_peptide_positions_from_peptide_dictionary(peptide, peptide_proteins, peptide_positions):
    """
    Find the positions of a peptide in a fasta file.

    :param fasta_dict: dictionary with fasta headers as keys and sequences as values
    :param peptide: peptide to search for
    :return: dictionary with fasta headers as keys and positions as values
    """
    if peptide in peptide_proteins.keys():
        return ";".join(peptide_proteins[peptide]), ";".join(peptide_positions[peptide])
    else:
        return '', ''


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
        if 'decoy' in key.lower() or 'reverse' in key.lower():
            new = 'REV_' + key.split('|')[1]
        else:
            new = key.split('|')[1]
        if (exclude is not None) and (new in exclude):
            print('excluded')
            del record_dict[key]
        else:
            record_dict[new] = value
            # protein prospector stores the protein description for decoys but not the protein - hence the additional lookup here
            del record_dict[key]
    return record_dict

def digest_proteins_tryptic(protein_dict={}, missed_cleavages=2, min_length=6, max_length=50):
    """
    Digest proteins with trypsin and returns twop dictionaries with peptides as keys. 
    The first one has the proteins as values and the second one the positions of the peptides in the proteins.
    
    :param protein_dict: dictionary with fasta headers as keys and sequences as values
    :param missed_cleavages: number of missed cleavages
    :param min_length: minimum length of peptides
    :param max_length: maximum length of peptides
    :return: proteins: dictionary with peptides as keys and proteins as values
    :return: positions: dictionary with peptides as keys and positions as values
    """
    proteins = {}
    positions = {}
    for cmc in range(missed_cleavages + 1):
        # preceded by either Kor R(?<=[KR]) 
        #   or beginning of sequence (?<=^) 
        #   or M at the beginning of the sequence (?<=^M)
        # the actual peptide is looked for with a lookahead (?=) containing a capturing group
        # this ensures that also overlaping peptides are matched
        # for 0 missed cleavages: the peptide lookahead looks like this: (?=([^KR]*(?:[KR]|$)))
        # meaning anything that is not a K or R followed by either K, R or the end of the sequence
        # for each misscleavage, the peptide lookahead is extended by [^KR]*[KR] - meaning it has to contain a K or R.
        pat = re.compile('(?:(?<=[KR])|(?<=^M)|(?<=^))(?=('+"".join(['[^KR]*[KR]'] * cmc) + '[^KR]*(?:[KR]|$)))')
        for key, value in protein_dict.items():
            sequence = str(value.seq)
            for match in re.finditer(pat, sequence):
                pep = match.group(1)
                start = match.start(1)
                if len(pep) >= min_length and len(pep) <= max_length:
                    if pep in proteins.keys():
                        proteins[pep].append(key)
                        positions[pep].append(str(start))
                    else:
                        proteins[pep] = [key]
                        positions[pep] = [str(start)]
    return proteins, positions

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


def plot_distribution(df, x, bins, ylim_top=None, ax=None):
    palette = sns.color_palette(['#8CCFC0', '#83320C', '#E364B4'])
    df['ID_type'] = df['entr_group']
    if ax is None:
        fig, ax = plt.subplots()
    sns.histplot(data=df, x=x, hue='ID_type', element="step", bins=bins, linewidth=1.7,
                 # hist_kws={"linewidth": 5, "alpha": 1},
                 hue_order=['E.coli', 'decoy', 'entrapment'], palette=palette, ax=ax)
    if ylim_top:
        ax.set_ylim(top=ylim_top)
    # plt.show()
    return ax
