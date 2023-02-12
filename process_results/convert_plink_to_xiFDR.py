# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 08:06:55 2020

@author: sgiese
"""

import glob
import pandas as pd
from tqdm import tqdm
from pyteomics import mgf, parser, fasta
import re
import numpy as np


def get_fasta_df(fastaf):
    """
    Get a fasta dataframe (peptide-protein).

    Parameters:
        fastaf: str, location of fasta.
    """
    print('Cleaving the proteins with trypsin...')

    peptides = []
    proteins = []
    # standard regex
    # r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    # allow cleavage before P
    tryp_regex = '([KR])|((?<=W)K(?=P))|((?<=M)R(?=P))'
    with open(fastaf, mode='rt') as gzfile:
        for description, sequence in tqdm(fasta.FASTA(gzfile)):
            # two trypsin rules
            new_peptides1 = parser.cleave(sequence, 'trypsin', missed_cleavages=2)
            new_peptides2 = parser.cleave(sequence, tryp_regex, 2)
            new_peptides = list(set(new_peptides1) | set(new_peptides2))
            peptides.extend(new_peptides)
            proteins.extend([description.split()[0]] * len(new_peptides))

            # for peptides where m is cleaved
            if sequence.startswith("M"):
                new_peptides_m3 = parser.cleave(sequence[1:], 'trypsin', missed_cleavages=2)
                new_peptides_m4 = parser.cleave(sequence[1:], tryp_regex, missed_cleavages=2)
                new_peptides_m = list(set(new_peptides_m3) | set(new_peptides_m4))
                peps_m = set(new_peptides_m) - set(new_peptides)
                peptides.extend(peps_m)
                proteins.extend([description.split()[0]] * len(peps_m))
    df_fasta = pd.DataFrame()
    df_fasta["Peptide"] = peptides
    df_fasta["Proteins"] = proteins
    df_fasta = df_fasta.set_index("Peptide")

    df_fasta = df_fasta.reset_index()
    df_fasta = df_fasta.groupby(["Peptide"], as_index=False).agg({"Proteins": ";".join})
    df_fasta = df_fasta.set_index("Peptide")
    return df_fasta


def get_protein_dic(fastaf):
    """
    Get a <protein>:<sequence> dictionary from a fasta file.

    Parameters:
        fastaf: str, location of fasta
    """
    proteins_dic = {}
    with open(fastaf, mode='rt') as gzfile:
        for description, sequence in tqdm(fasta.FASTA(gzfile)):
            proteins_dic[description.split()[0]] = sequence
            # proteins_dic["|".join(description.split()[0].split("|")[0:2])] = sequence
    return proteins_dic


# adjust for multiple proteins
def get_positions(peptide, proteins, proteins_dic):
    """
    Get a string of peptide positions.

    Parameters:
    peptide: str, peptide
    proteins: str,proteins sep by ;
    proteins_dic: dict, <protein>:<sequence>
    """
    return ";".join([str(proteins_dic[protein_i].index(peptide) + 1) for protein_i in proteins.split(";")])


def create_modseq(df_plink_cl):
    """
    Create a modified peptide sequence string

    """

    # make mod detection easier
    df_plink_cl["Modifications"] = df_plink_cl["Modifications"].fillna("")

    # psm_df = df_plink_cl[df_plink_cl["Modifications"] != ""].head(1000)
    # map common modifications
    mod_dic = {"M": "ox", "C": "cm"}
    modex = re.compile(r"\[(\w)\]\((\d+)\)")
    modsseqs1 = [""] * len(df_plink_cl)
    modsseqs2 = [""] * len(df_plink_cl)

    # compute the mod seqecne column
    idx = -1
    for ii, row in tqdm(df_plink_cl.iterrows(), total=len(df_plink_cl)):
        idx += 1
        # get modifications
        mods = row["Modifications"]

        if mods == "":
            # nod mods
            modsseqs1[idx] = row.Peptide1
            modsseqs2[idx] = row.Peptide2
            continue
        else:
            # modifications list with (aa, pos) tuples
            # go reverse because otherwise the indices dont match after the
            # first insertion
            pepseq1 = row.Peptide1
            pepseq2 = row.Peptide2
            modified1 = False
            modified2 = False
            modes_ar = re.findall(modex, mods)[::-1]
            for mod in modes_ar:
                modpos = int(mod[1])

                # mods on peptide 1
                if modpos <= row["PepLength1"]:
                    pepseq1 = pepseq1[:modpos] + mod_dic[mod[0]] + pepseq1[modpos:]
                    modified1 = True
                else:
                    # -3 for (, ), - in the peptide string?!
                    modpos_adj = modpos - row.PepLength1 - 3
                    # mods on peptide 2
                    pepseq2 = pepseq2[:modpos_adj] + mod_dic[mod[0]] + pepseq2[modpos_adj:]
                    modified2 = True

        # check if there were changes and if not, write the unmodified sequence
        if modified1:
            modsseqs1[idx] = pepseq1
        else:
            modsseqs1[idx] = row.Peptide1

        if modified2:
            modsseqs2[idx] = pepseq2
        else:
            modsseqs2[idx] = row.Peptide2

    df_plink_cl["ModSeq1"] = modsseqs1
    df_plink_cl["ModSeq2"] = modsseqs2


# get peptide1,2,linkpos
def split_peptides(df_plink_cl):
    regex = re.compile(r"(\w+)\((\d+)\)-(\w+)\((\d+)\)")
    df_info = df_plink_cl["Peptide"].str.extract(regex)
    df_info.columns = ["Peptide1", "LinkPos1", "Peptide2", "LinkPos2"]
    df_plink_cl = df_plink_cl.join(df_info)
    return df_plink_cl


def split_proteins(df_plink_cl):
    regex = re.compile(r"(\S+) \((\d+)\)-(\S+) \((\d+)\)")
    df_info = df_plink_cl["Proteins"].str.extract(regex)
    df_info.columns = ["Protein1", "PepPos1", "Protein2", "PepPos2"]
    df_plink_cl = df_plink_cl.join(df_info)
    return df_plink_cl


def decoy_type_annotation(df_plink_cl):
    # Target_Decoy: the identification is target or decoy.
    # 0 for Decoy-Decoy,
    # 1 for Target-Decoy (or Decoy-Target), and
    # 2 for Target-Target.
    TD_dic = {2: "TT", 1: "TD", 0: "DD"}
    try:
        df_plink_cl["IDType"] = df_plink_cl["Target_Decoy"].map(TD_dic)
        df_plink_cl["isTT"] = df_plink_cl["IDType"] == "TT"
        df_plink_cl["isTD"] = df_plink_cl["IDType"] == "TD"
        df_plink_cl["isDD"] = df_plink_cl["IDType"] == "DD"
    except KeyError:
        df_plink_cl["isTT"] = True
        df_plink_cl["isTD"] = False
        df_plink_cl["isDD"] = False


def group_annotation(df_plink_cl):
    # convert within / between
    # Protein_Type: the same as the Protein_Type in spectra level described above,
    # but with 0 for Regular/Common, 1 for Intra-protein, and 2 for Inter-protein.
    type_dict = {1: "within", 2: "between", 'Intra-Protein': 'within', 'Inter-Protein': 'between'}
    df_plink_cl["Group"] = [type_dict[t] for t in df_plink_cl["Protein_Type"]]


def get_protein_short(prot):
    """Get short protein identifier to work with xiFDR."""
    return (";".join([i.split("|")[0] for i in prot.replace("sp|", "").split(";")]))


def process_plink(result_file):
    output = result_file.replace(".csv", "_xifdr.csv")

    # plink file
    df_plink = pd.read_csv(result_file)
    df_plink["PSMID"] = np.arange(1, len(df_plink) + 1)
    print(df_plink.shape)

    # this is needed if the filtered data is used as input
    if "filtered_cross-linked_" in result_file:
        df_plink["Peptide_Type"] = 3
        df_plink_cl = df_plink[df_plink["Peptide_Type"] == 3]
        df_plink_cl = df_plink_cl.sort_values(by="Score", ascending=False)
    else:
        # get crosslinks, only. 0=linear, 1=monolinked, 2=loop-linked, 3=crosslinked
        df_plink_cl = df_plink[df_plink["Peptide_Type"] == 3]
        df_plink_cl = df_plink_cl.sort_values(by="Q-value", ascending=True)

    print("Reorganize Peptide/protein storage...")
    df_plink_cl = split_peptides(df_plink_cl)
    df_plink_cl = split_proteins(df_plink_cl)
    df_plink_cl["PepLength1"] = df_plink_cl["Peptide1"].apply(len)
    df_plink_cl["PepLength2"] = df_plink_cl["Peptide2"].apply(len)

    print("Decoy Annotation..")
    # isTT, isTD, columns
    decoy_type_annotation(df_plink_cl)

    print("Group Annotation..")
    # self, between annotation
    group_annotation(df_plink_cl)

    # write normal crosslink dataframe

    df_plink_cl_between = df_plink_cl[(df_plink_cl['Protein_Type'] == 'Inter-Protein') | (df_plink_cl['Protein_Type'] == 2)]

    df_plink_cl.to_csv(output)

    # # TODO all these are modifications necessary to put it into xifdr for PPIs - taken out for now
    # # FASTA
    # fastaf = "/Users/lenz/Documents/entrapment-control/database/ecoli+entrapment+decoys.fasta"
    #
    # # get sequence data peptides, proteins
    # print("Sequence analysis...")
    # df_fasta = get_fasta_df(fastaf)
    # proteins_dic = get_protein_dic(fastaf)
    #
    # # target / decoy proteins
    # df_plink_cl = df_plink_cl.merge(df_fasta, left_on="Peptide1", right_index=True, suffixes=("", "_1"))
    # df_plink_cl = df_plink_cl.merge(df_fasta, left_on="Peptide2", right_index=True, suffixes=("", "_2"))
    #
    # # some weird list formatting was introduce .. reavel to single entries
    # df_plink_cl["Protein1"] = np.ravel(df_plink_cl.Proteins_1)
    # df_plink_cl["Protein2"] = np.ravel(df_plink_cl.Proteins_2)
    #
    # df_plink_cl["Protein1_short"] = df_plink_cl["Protein1"].apply(get_protein_short)
    # df_plink_cl["Protein2_short"] = df_plink_cl["Protein2"].apply(get_protein_short)
    #
    # # TODO this does not take into account ambiguity - needs fixing
    # df_plink_cl["Decoy1"] = df_plink_cl["Protein1"].str.contains("REVERSE")
    # df_plink_cl["Decoy2"] = df_plink_cl["Protein2"].str.contains("REVERSE")
    #
    # #(((df_plink_cl['Decoy2'] + df_plink_cl['Decoy1']) == 0) == (df_plink_cl['isTT'])).all()
    #
    # # create the modification seqeuences
    # print("Creating Modified Sequences...")
    # create_modseq(df_plink_cl)
    #
    # # peptide positions
    # pep_pos1 = np.array([get_positions(pep, prot, proteins_dic) for prot, pep in
    #                      zip(df_plink_cl["Protein1"], df_plink_cl["Peptide1"])])
    #
    # pep_pos2 = np.array([get_positions(pep, prot, proteins_dic) for prot, pep in
    #                      zip(df_plink_cl["Protein2"], df_plink_cl["Peptide2"])])
    #
    # # assign to dataframe, these are peptide positions! not protein link positions
    # df_plink_cl["PepPos1"] = pep_pos1
    # df_plink_cl["PepPos2"] = pep_pos2
    # df_plink_cl["Description1"] = df_plink_cl.Proteins_1
    # df_plink_cl["Description2"] = df_plink_cl.Proteins_2
    # df_plink_cl["Fasta1"] = df_plink_cl.Proteins_1
    # df_plink_cl["Fasta2"] = df_plink_cl.Proteins_2
    # df_plink_cl["Run"] = df_plink_cl["Title"]
    # try:
    #     df_plink_cl["exp m/z"] = df_plink_cl["Precursor_MH"] / df_plink_cl["Charge"]
    # except KeyError:
    #     df_plink_cl["Precursor_MH"] = df_plink_cl["Precursor_Mass"]
    #     df_plink_cl["exp m/z"] = df_plink_cl["Precursor_MH"] / df_plink_cl["Charge"]
