from process_results import *
from process_results.process_plink2 import process_plink2_ppi
from plots_and_functions import fasta_to_dict, digest_proteins_tryptic
from xiFDRWrapper.XiFdrWrapper import XiFdrWrapper
from pathlib import Path
import os
import pandas as pd
import logging
import sys
import re
from score_names import score_names


def read_single_xiFDR_result_file(subpath, name, file, ecoli_proteins, FDR_Folder='FDR'):
    """
    Read in the xiFDR ppi results.

    :param subpath: path to the file
    :param name: name of the search engine
    :param file: regex to match the file
    :return: DataFrame with the results
    """

    basepath = os.path.join('search_results', subpath, FDR_Folder)
    files = [f for f in os.listdir(basepath) if re.match(file, f)]
    PPI_df = process_xisearch(os.path.join(basepath, files[0]), ecoli_proteins)
    PPI_df['search_engine'] = name
    return PPI_df




def read_all_xiFDR_result_files(file, ecoli_proteins, FDR_Folders={}):
    """
    Read in the xiFDR ppi results.

    :param file: regex to match the file
    :return: DataFrame with the results
    """

    # read in the xiFDR ppi results
    #kojak
    FDR_Folder = FDR_Folders.get("kojak", "FDR")
    kojak_PPI_df = read_single_xiFDR_result_file('kojak', 'Kojak', file, ecoli_proteins, FDR_Folder)

    # mango
    FDR_Folder = FDR_Folders.get("mango", "FDR")
    mango_PPI_df = read_single_xiFDR_result_file('mango_comet_xlinkprophet', 'MangoCometXlinkProphet', file, ecoli_proteins, FDR_Folder)

    # plink2 cleavable
    FDR_Folder = "FDR"
    if "plink2" in FDR_Folders:
        FDR_Folder = FDR_Folders["plink2"]
    plink2_cleav_PPI_df = read_single_xiFDR_result_file('plink2_cleav', 'pLink2_cleav', file, ecoli_proteins, FDR_Folder)

    # plink2 non-cleavable
    FDR_Folder = "FDR"
    if "plink2" in FDR_Folders:
        FDR_Folder = FDR_Folders["plink2"]
    plink2_nonCleav_PPI_df = read_single_xiFDR_result_file('plink2_nonCleav', 'pLink2_nonCleav', file, ecoli_proteins, FDR_Folder)

    # proteinprospector
    FDR_Folder = FDR_Folders.get("pp", "FDR")
    pp_PPI_df = read_single_xiFDR_result_file('proteinprospector', 'ProteinProspector', file, ecoli_proteins, FDR_Folder)

    # xisearch
    FDR_Folder = FDR_Folders.get("xisearch", "FDR")
    xisearch_PPI_df = read_single_xiFDR_result_file('xisearch', 'xiSEARCH', file, ecoli_proteins, FDR_Folder)

    # xlinkx
    FDR_Folder = FDR_Folders.get("xlinkx", "FDR")
    xlinkx_PPI_df = read_single_xiFDR_result_file('xlinkx', 'XlinkX', file, ecoli_proteins, FDR_Folder)

    # merox
    FDR_Folder = FDR_Folders.get("merox", "FDR")
    merox_PPI_df = read_single_xiFDR_result_file('merox', 'MeroX', file, ecoli_proteins, FDR_Folder)

    # openpepxl
    FDR_Folder = FDR_Folders.get("openpepxl", "FDR")
    openpepxl_PPI_df = read_single_xiFDR_result_file('openpepxl', 'OpenPepXL', file, ecoli_proteins, FDR_Folder)

    # concatenate all dataframes

    return pd.concat([kojak_PPI_df, mango_PPI_df, plink2_cleav_PPI_df, plink2_nonCleav_PPI_df,
                      pp_PPI_df, xisearch_PPI_df, xlinkx_PPI_df, merox_PPI_df, openpepxl_PPI_df],
                    ignore_index=True)




# read in E. coli fasta file
ecoli_fasta = fasta_to_dict('database/uniprot-k12-filtered-proteome_UP000000625.fasta')
# read in uniprot IDs of E. coli as list
ecoli_proteins = list(ecoli_fasta.keys())
full_fasta = fasta_to_dict('database/ecoli_human_comb_rever.fasta')

all_peptide_proteins, all_peptide_positions = digest_proteins_tryptic(
    full_fasta, missed_cleavages=7, min_length=1, max_length=80)



xiFDR_jar = os.path.join('xiFDR', 'xiFDR-2.2.betaB-jar-with-dependencies.jar')
add_xifdr_args = ["--uniquePSMs=0", "--minPeptideLength=1", "--forward=" + "|".join(score_names.values())]
#
# process the results
#
# Kojak
kojak_result_file = os.path.join('search_results', 'kojak', 'Kojak_CSMs_NoLoopLinks.tsv')
kojak_fasta_file = os.path.join('search_results', 'kojak', 'EcoliDECOY_HumanENTRAP.fasta')
kojak_df = process_kojak(kojak_result_file, ecoli_proteins, kojak_fasta_file)

kojak_xifdr_output = os.path.join('search_results', 'kojak', 'FDR')
Path.mkdir(Path(kojak_xifdr_output), exist_ok=True, parents=True)
kojak_xifdr_input = os.path.join(kojak_xifdr_output, 'kojak_xiFDR_Input.csv')
kojak_df.to_csv(kojak_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(kojak_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(kojak_xifdr_input, kojak_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# Mango Comet XlinkProphet
mango_result_file = os.path.join('search_results', 'mango_comet_xlinkprophet',
                                 'iprophet-withdecoys-xl-filtered_2pct.xlsx')
mango_df = process_mango_comet_xlinkprophet(mango_result_file, ecoli_proteins, full_fasta)

mango_xifdr_output = os.path.join('search_results', 'mango_comet_xlinkprophet', 'FDR')
Path.mkdir(Path(mango_xifdr_output), exist_ok=True, parents=True)
mango_xifdr_input = os.path.join(mango_xifdr_output, 'mango_xiFDR_Input.csv')
mango_df.to_csv(mango_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(mango_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(mango_xifdr_input, mango_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# pLink2 cleavable
plink2_cleave_res_folder = os.path.join('search_results', 'plink2_cleav')
plink2_cleav_result_file = os.path.join(
    plink2_cleave_res_folder,
    'LK_2.3.9_cleavable_ecoli_human_2023.01.02.filtered_cross-linked_spectra.csv')
plink2_cleav_unfiltered_result_file = os.path.join(
    plink2_cleave_res_folder, 'LK_2.3.9_cleavable_ecoli_human_2023.01.02.csv')
plink2_cleav_df = process_plink2(plink2_cleav_result_file, plink2_cleav_unfiltered_result_file,
                                 ecoli_proteins, all_peptide_proteins, all_peptide_positions)
plink2_cleav_df['search_engine'] = 'pLink2_cleav'

plink2_cleav_xifdr_output = os.path.join('search_results', 'plink2_cleav', 'FDR')
Path.mkdir(Path(plink2_cleav_xifdr_output), exist_ok=True, parents=True)
plink2_cleav_xifdr_input = os.path.join(plink2_cleav_xifdr_output, 'plink2_cleav_xiFDR_Input.csv')
plink2_cleav_df.to_csv(plink2_cleav_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(plink2_cleav_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(plink2_cleav_xifdr_input, plink2_cleav_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# ToDo: (not used for now) use for residue pair to PPI analysis
# plink2_cleav_residue_ppi_file = os.path.join('search_results', 'plink2_cleav',
#                               'LK_2.3.9_cleavable_ecoli_human_2023.01.02.filtered_cross-linked_sites.csv')

# plink2_cleav_ppi_df, plink2_ppi_summary = process_plink2_ppi(plink2_cleav_residue_ppi_file,
#                                                              ecoli_proteins)
# plink2_cleav_ppi_df.to_csv(os.path.join(
#     'search_results', 'plink2_cleav', '_extracted_ppi.csv'))


# pLink2 nonCleaveable
plink2_nonCleave_res_folder = os.path.join('search_results', 'plink2_nonCleav')
plink2_nonCleav_result_file = os.path.join(
    plink2_nonCleave_res_folder,
    'ecoli_human_2023.01.24.filtered_cross-linked_spectra.csv')
# unfiltered is too big for repository adjust datastore path if needed
plink2_nonCleav_unfiltered_result_folder = \
    '/data/rappstore/users/lfischer/search_comparison/results/plink2_LK/2.3.11_nonCleavable/reports/'
plink2_nonCleav_unfiltered_result_file = os.path.join(
    plink2_nonCleav_unfiltered_result_folder, 'ecoli_human_2023.01.24.csv')

plink2_nonCleav_df = process_plink2(plink2_nonCleav_result_file,
                                    plink2_nonCleav_unfiltered_result_file,
                                    ecoli_proteins, all_peptide_proteins, all_peptide_positions)
plink2_nonCleav_df['search_engine'] = 'pLink2_nonCleav'

plink2_nonCleav_xifdr_output = os.path.join('search_results', 'plink2_nonCleav', 'FDR')
Path.mkdir(Path(plink2_nonCleav_xifdr_output), exist_ok=True, parents=True)
plink2_nonCleav_xifdr_input = os.path.join(plink2_nonCleav_xifdr_output, 'plink2_nonCleav_xiFDR_Input.csv')
plink2_nonCleav_df.to_csv(plink2_nonCleav_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(plink2_nonCleav_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(plink2_nonCleav_xifdr_input, plink2_nonCleav_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )


# proteinProspector
pp_result_file = os.path.join('search_results', 'proteinprospector', 'ec1_dsso_CSM_2p_inter.csv')
pp_df = process_proteinprospector(pp_result_file, ecoli_proteins, full_fasta)

pp_xifdr_output = os.path.join('search_results', 'proteinprospector', 'FDR')
Path.mkdir(Path(pp_xifdr_output), exist_ok=True, parents=True)
pp_xifdr_input = os.path.join(pp_xifdr_output, 'proteinprospector_xiFDR_Input.csv')
pp_df.to_csv(pp_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(pp_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(pp_xifdr_input, pp_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# xiSEARCH
xisearch_result_file = os.path.join(
    'search_results', 'xisearch',
    'xiSEARCH_CSM_FDR_proteome+human3_linMod_filtered_0_PSM_xiFDR1.3.36.csv')
#xisearch_result_file = "/data/rappstore/users/lfischer/search_comparison/results/xisearch/PPI_FDR_proteome+human3_linMod_filtered_0_PSM_xiFDR1.3.36.csv"
xisearch_df = process_xisearch(xisearch_result_file, ecoli_proteins)

xisearch_xifdr_output = os.path.join('search_results', 'xisearch', 'FDR')
Path.mkdir(Path(xisearch_xifdr_output), exist_ok=True, parents=True)
xisearch_xifdr_input = os.path.join(xisearch_xifdr_output, 'xisearch_xiFDR_Input.csv')
xisearch_df.to_csv(xisearch_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(xisearch_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(xisearch_xifdr_input, xisearch_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# ToDo: not used for now
# xisearch_xifdr_output = os.path.join('search_results', 'xisearch', 'PPI_FDR')
# xifdr.xifdr_execution(xisearch_xifdr_input, xisearch_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
#                       additional_xifdr_arguments=["--ppifdr=2",], xifdr_filename=xiFDR_jar
#                       )

# XlinkX
xlinkx_result_file = os.path.join('search_results', 'xlinkx', '2p_all_files_plus3_CSMs.txt')
xlinkx_decoy_file = os.path.join('search_results', 'xlinkx', '2p_all_files_plus3_DecoyCSMs.txt')
xlinkx_df = process_xlinkx(xlinkx_result_file, xlinkx_decoy_file, ecoli_proteins, all_peptide_proteins,all_peptide_positions)

xlinkx_xifdr_output = os.path.join('search_results', 'xlinkx', 'FDR')
Path.mkdir(Path(xlinkx_xifdr_output), exist_ok=True, parents=True)
xlinkx_xifdr_input = os.path.join(xlinkx_xifdr_output, 'xlinkx_xiFDR_Input.csv')
xlinkx_df.to_csv(xlinkx_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(xlinkx_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(xlinkx_xifdr_input, xlinkx_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# MeroX
merox_result_file = os.path.join('search_results', 'merox', '2p_FDR_230113.csv')
merox_df = process_merox(merox_result_file, ecoli_proteins)

merox_xifdr_output = os.path.join('search_results', 'merox', 'FDR')
Path.mkdir(Path(merox_xifdr_output), exist_ok=True, parents=True)
merox_xifdr_input = os.path.join(merox_xifdr_output, 'merox_xiFDR_Input.csv')
merox_df.to_csv(merox_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(merox_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(merox_xifdr_input, merox_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )

# OpenPepXL
openpepxl_result_file = os.path.join('search_results', 'openpepxl', 'OpenPepXL_DSSO.csv')
openpepxl_df = process_openpepxl(openpepxl_result_file, ecoli_proteins)

openpepxl_xifdr_output = os.path.join('search_results', 'openpepxl', 'FDR')
Path.mkdir(Path(openpepxl_xifdr_output), exist_ok=True, parents=True)
openpepxl_xifdr_input = os.path.join(openpepxl_xifdr_output, 'openpepxl_xiFDR_Input.csv')
openpepxl_df.to_csv(openpepxl_xifdr_input)
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join(openpepxl_xifdr_output, 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
xifdr.xifdr_execution(openpepxl_xifdr_input, openpepxl_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
                      )


# concatenate all dataframes
df = pd.concat([kojak_df, mango_df, plink2_cleav_df, plink2_nonCleav_df,
                pp_df, xisearch_df, xlinkx_df, merox_df, openpepxl_df],
               ignore_index=True)

# save to csv
# join the xiFDR results for each level
df_csm = read_all_xiFDR_result_files(file=r'.*_CSM_.*.csv', ecoli_proteins=ecoli_proteins)
df_csm.to_csv(os.path.join('search_results', 'processed_CSM_results.csv'), index=False)

df_peppair = read_all_xiFDR_result_files(file=r'.*_PeptidePairs_.*.csv', ecoli_proteins=ecoli_proteins)
df_peppair.to_csv(os.path.join('search_results', 'processed_peppair_results.csv'), index=False)

df_respair = read_all_xiFDR_result_files(file=r'.*_Links_.*.csv', ecoli_proteins=ecoli_proteins)
df_respair.to_csv(os.path.join('search_results', 'processed_respair_results.csv'), index=False)

df_ppi = read_all_xiFDR_result_files(file=r'.*_ppi_.*.csv', ecoli_proteins=ecoli_proteins)
df_ppi.to_csv(os.path.join('search_results', 'processed_ppi_results.csv'), index=False)


df_ppi = read_all_xiFDR_result_files(file=r'.*_ppi_.*.csv', ecoli_proteins=ecoli_proteins)
df_ppi.to_csv(os.path.join('search_results', 'processed_ppi_results.csv'), index=False)

# ToDo: Not used for now
# # join the xiFDR results for each level based on highest FDR
# df_csm = read_all_xiFDR_result_files(file=r'.*_CSM_.*.csv', ecoli_proteins=ecoli_proteins, FDR_Folders={'xisearch': 'PPI_FDR'})
# df_csm.to_csv(os.path.join('search_results', 'ppi_processed_CSM_results.csv'), index=False)
#
#
# df_peppair = read_all_xiFDR_result_files(file=r'.*_PeptidePairs_.*.csv', ecoli_proteins=ecoli_proteins, FDR_Folders={'xisearch': 'PPI_FDR'})
# df_peppair.to_csv(os.path.join('search_results', 'ppi_processed_peppair_results.csv'), index=False)
#
# df_respair = read_all_xiFDR_result_files(file=r'.*_Links_.*.csv', ecoli_proteins=ecoli_proteins, FDR_Folders={'xisearch': 'PPI_FDR'})
# df_respair.to_csv(os.path.join('search_results', 'ppi_processed_respair_results.csv'), index=False)
#
# df_ppi = read_all_xiFDR_result_files(file=r'.*_ppi_.*.csv', ecoli_proteins=ecoli_proteins, FDR_Folders={'xisearch': 'PPI_FDR'})
# df_ppi.to_csv(os.path.join('search_results', 'ppi_processed_ppi_results.csv'), index=False)
#
#
# df_ppi_base = read_all_xiFDR_result_files(file=r'.*_ppi_.*.csv', ecoli_proteins=ecoli_proteins, FDR_Folders={'xisearch': 'PPI_FDR'})
# df_ppi_base.to_csv(os.path.join('search_results', 'ppi_processed_ppi_results.csv'), index=False)
