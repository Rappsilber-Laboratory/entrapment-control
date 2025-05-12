from process_results import *
from process_results.process_mango_comet_xlinkprophet import process_mango_comet_xlinkprophet, process_mango_comet_xlinkprophet_ppi
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
    files = [f for f in os.listdir(basepath) if re.match(file, f) and not re.match(".*_NAPs_.*", f)]
    PPI_df = process_xisearch(os.path.join(basepath, files[0]), ecoli_proteins)
    PPI_df['search_engine'] = name
    return PPI_df

def file_exists(path, pattern):
    if not os.path.exists(path):
        return False
    return len([f for f in os.listdir(path) if re.match(pattern, f)]) > 0



def read_all_xiFDR_result_files(file, ecoli_proteins, FDR_Folders={}, default_FDR_Folder='FDR'):
    """
    Read in the xiFDR ppi results.

    :param file: regex to match the file
    :return: DataFrame with the results
    """

    # read in the xiFDR ppi results
    #kojak
    FDR_Folder = FDR_Folders.get("kojak", default_FDR_Folder)
    kojak_PPI_df = read_single_xiFDR_result_file('kojak', 'Kojak', file, ecoli_proteins, FDR_Folder)

    # mango
    FDR_Folder = FDR_Folders.get("mango", default_FDR_Folder)
    mango_PPI_df = read_single_xiFDR_result_file('mango_comet_xlinkprophet', 'MangoCometXlinkProphet', file, ecoli_proteins, FDR_Folder)

    ## plink2 cleavable
    #FDR_Folder = default_FDR_Folder
    #if "plink2" in FDR_Folders:
    #    FDR_Folder = FDR_Folders["plink2"]
    #plink2_cleav_PPI_df = read_single_xiFDR_result_file('plink2_cleav', 'pLink2_cleav', file, ecoli_proteins, FDR_Folder)

    # plink2 non-cleavable
    FDR_Folder = default_FDR_Folder
    if "plink2" in FDR_Folders:
        FDR_Folder = FDR_Folders["plink2"]
    plink2_nonCleav_PPI_df = read_single_xiFDR_result_file('plink2_nonCleav', 'pLink2_nonCleav', file, ecoli_proteins, FDR_Folder)

    # proteinprospector
    FDR_Folder = FDR_Folders.get("pp", default_FDR_Folder)
    pp_PPI_df = read_single_xiFDR_result_file('proteinprospector', 'ProteinProspector', file, ecoli_proteins, FDR_Folder)

    # xisearch
    FDR_Folder = FDR_Folders.get("xisearch", default_FDR_Folder)
    xisearch_PPI_df = read_single_xiFDR_result_file('xisearch', 'xiSEARCH', file, ecoli_proteins, FDR_Folder)

    # xlinkx
    FDR_Folder = FDR_Folders.get("xlinkx", default_FDR_Folder)
    xlinkx_PPI_df = read_single_xiFDR_result_file('xlinkx', 'XlinkX PD2.5', file, ecoli_proteins, FDR_Folder)

    # xlinkx_pd32
    FDR_Folder = FDR_Folders.get("xlinkx_pd32", default_FDR_Folder)
    xlinkx_pre_PPI_df = read_single_xiFDR_result_file('xlinkx_pd32', 'XlinkX PD3.2', file, ecoli_proteins, FDR_Folder)

    # merox
    FDR_Folder = FDR_Folders.get("merox", default_FDR_Folder)
    merox_PPI_df = read_single_xiFDR_result_file('merox', 'MeroX', file, ecoli_proteins, FDR_Folder)

    # openpepxl
    FDR_Folder = FDR_Folders.get("openpepxl", default_FDR_Folder)
    openpepxl_PPI_df = read_single_xiFDR_result_file('openpepxl', 'OpenPepXL', file, ecoli_proteins, FDR_Folder)

    # concatenate all dataframes

    return pd.concat([kojak_PPI_df, mango_PPI_df, plink2_nonCleav_PPI_df,# plink2_cleav_PPI_df,
                      pp_PPI_df, xisearch_PPI_df, xlinkx_PPI_df, xlinkx_pre_PPI_df, merox_PPI_df, openpepxl_PPI_df],
                    ignore_index=True)




# read in E. coli fasta file
ecoli_fasta = fasta_to_dict('database/uniprot-k12-filtered-proteome_UP000000625.fasta')
# read in uniprot IDs of E. coli as list
ecoli_proteins = list(ecoli_fasta.keys())
full_fasta = fasta_to_dict('database/ecoli_human_comb_rever.fasta')

all_peptide_proteins, all_peptide_positions = digest_proteins_tryptic(
    full_fasta, missed_cleavages=7, min_length=1, max_length=80)

#df_csm = read_all_xiFDR_result_files(file=r'.*_CSM_.*.csv', ecoli_proteins=ecoli_proteins)


xiFDR_jar = os.path.join('xiFDR', 'xiFDR-2.2-jar-with-dependencies.jar')
add_xifdr_args = ["--uniquePSMs=0", "--lengthgroups=0", "--minPeptideLength=1", "--forward=iProbability|E_value|Sum_score|" + "|".join(score_names.values())] # "--ec-filter", 
xifdr = XiFdrWrapper()
xifdr.logger.addHandler(logging.FileHandler(os.path.join('search_results', 'xiFDR.log')))
xifdr.logger.addHandler(logging.StreamHandler(sys.stdout))
xifdr.logger.setLevel(logging.INFO)
add_xifdr_args_forward = ["--uniquePSMs=0", "--lengthgroups=0", "--minPeptideLength=1", "--forward=iProbability|E_value|Sum_score|" + "|".join(score_names.values())]
add_xifdr_args_ppi = add_xifdr_args_forward + ["--ppifdr=2"]


#
# process the results
#
# Kojak
#kojak_result_file = os.path.join('search_results', 'kojak', 'Kojak_CSMs_NoLoopLinks.tsv')
kojak_result_file = os.path.join('search_results', 'kojak', 'results.csm.inter.txt')
kojak_fasta_file = os.path.join('search_results', 'kojak', 'EcoliDECOY_HumanENTRAP.fasta')
kojak_df = process_kojak(kojak_result_file, ecoli_proteins, kojak_fasta_file)
kojak_xifdr_output = os.path.join('search_results', 'kojak', 'FDR')
Path.mkdir(Path(kojak_xifdr_output), exist_ok=True, parents=True)
kojak_xifdr_input = os.path.join(kojak_xifdr_output, 'kojak_xiFDR_Input.csv')
kojak_df.to_csv(kojak_xifdr_input)
if not file_exists(kojak_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(kojak_xifdr_input, kojak_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# Kojak CSM to 2%PPI FDR
kojak_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'kojak', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(kojak_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(kojak_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(kojak_xifdr_input, kojak_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                          )

# highest to PPI
kojak_ppi_result_file = os.path.join('search_results', 'kojak', 'results.ppi.inter.txt')
kojak_ppi_df = process_kojak(kojak_ppi_result_file, ecoli_proteins, kojak_fasta_file)
kojak_ppi_xifdr_output = os.path.join('search_results', 'kojak', 'FDR_PPI')
Path.mkdir(Path(kojak_ppi_xifdr_output), exist_ok=True, parents=True)
kojak_ppi_xifdr_input = os.path.join(kojak_ppi_xifdr_output, 'kojak_PPI_xiFDR_Input.csv')
kojak_ppi_df.to_csv(kojak_ppi_xifdr_input)
if not file_exists(kojak_ppi_xifdr_output, ".*_ppi_xiFDR.*"):
    xifdr.xifdr_execution(kojak_ppi_xifdr_input, kojak_ppi_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )


# Mango Comet XlinkProphet
mango_result_file = os.path.join('search_results', 'mango_comet_xlinkprophet',
                                 'iprophet-withdecoys-xl-filtered_2pct.xlsx')
mango_df = process_mango_comet_xlinkprophet(mango_result_file, ecoli_proteins, full_fasta)

mango_xifdr_output = os.path.join('search_results', 'mango_comet_xlinkprophet', 'FDR')
Path.mkdir(Path(mango_xifdr_output), exist_ok=True, parents=True)
mango_xifdr_input = os.path.join(mango_xifdr_output, 'mango_xiFDR_Input.csv')
mango_df.to_csv(mango_xifdr_input)
if not file_exists(mango_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(mango_xifdr_input, mango_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# Mango CSM to 2%PPI FDR
mango_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'mango_comet_xlinkprophet', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(mango_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(mango_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(mango_xifdr_input, mango_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                            additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                            )

# Mango Highest to PPI
mango_ppi_result_file = os.path.join('search_results', 'mango_comet_xlinkprophet',
                                     'iprophet-withdecoys-xl-protpairs-filtered_2pct.tsv')
mango_ppi_output = os.path.join('search_results', 'mango_comet_xlinkprophet', 'FDR_PPI')

mango_ppi_df = process_mango_comet_xlinkprophet_ppi(mango_ppi_result_file, ecoli_proteins, full_fasta)
Path.mkdir(Path(mango_ppi_output), exist_ok=True, parents=True)
mango_ppi_df.to_csv(os.path.join(mango_ppi_output, 'mango_ppi_FDR.csv'))

## pLink2 cleavable
#plink2_cleave_res_folder = os.path.join('search_results', 'plink2_cleav')
#plink2_cleav_result_file = os.path.join(
#    plink2_cleave_res_folder,
#    'LK_2.3.9_cleavable_ecoli_human_2023.01.02.filtered_cross-linked_spectra.csv')
#plink2_cleav_unfiltered_result_file = os.path.join(
#    plink2_cleave_res_folder, 'LK_2.3.9_cleavable_ecoli_human_2023.01.02.csv')
#plink2_cleav_df = process_plink2(plink2_cleav_result_file, plink2_cleav_unfiltered_result_file,
#                                 ecoli_proteins, all_peptide_proteins, all_peptide_positions)
#plink2_cleav_df['search_engine'] = 'pLink2_cleav'
#
#plink2_cleav_xifdr_output = os.path.join('search_results', 'plink2_cleav', 'FDR')
#Path.mkdir(Path(plink2_cleav_xifdr_output), exist_ok=True, parents=True)
#plink2_cleav_xifdr_input = os.path.join(plink2_cleav_xifdr_output, 'plink2_cleav_xiFDR_Input.csv')
#plink2_cleav_df.to_csv(plink2_cleav_xifdr_input)
#xifdr.xifdr_execution(plink2_cleav_xifdr_input, plink2_cleav_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
#                      additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
#                      )
## plink2 cleavable CSM to 2%PPI FDR
#plink2_cleav_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'plink2_cleav', 'FDR_CSM_TO_2PercentPPI')
#Path.mkdir(Path(plink2_cleav_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
#add_xifdr_args = ["--uniquePSMs=0", "--minPeptideLength=1", "--ppifdr=2", "--forward=" + "|".join(score_names.values())]
#xifdr.xifdr_execution(plink2_cleav_xifdr_input, plink2_cleav_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
#                        additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
#                        )


# pLink2 nonCleaveable
plink2_nonCleave_res_folder = os.path.join('search_results', 'plink2_nonCleav')
plink2_nonCleav_result_file = os.path.join(
    plink2_nonCleave_res_folder,
    'ecoli_human_2023.01.24.filtered_cross-linked_spectra.csv')
plink2_nonCleav_unfiltered_result_file = os.path.join(
    plink2_nonCleave_res_folder,
    'ecoli_human_2023.01.24.unfiltered_cross-linked_spectra.csv')
## unfiltered is too big for repository adjust datastore path if needed
# plink2_nonCleav_unfiltered_result_folder = \
#     '/data/rappstore/users/lfischer/search_comparison/results/plink2_LK/2.3.11_nonCleavable/reports/'
# plink2_nonCleav_unfiltered_result_file = os.path.join(
#     plink2_nonCleav_unfiltered_result_folder, 'ecoli_human_2023.01.24.csv')

plink2_nonCleav_df = process_plink2(plink2_nonCleav_result_file,
                                    plink2_nonCleav_unfiltered_result_file,
                                    ecoli_proteins, all_peptide_proteins, all_peptide_positions)
plink2_nonCleav_df['search_engine'] = 'pLink2_nonCleav'

plink2_nonCleav_xifdr_output = os.path.join('search_results', 'plink2_nonCleav', 'FDR')
Path.mkdir(Path(plink2_nonCleav_xifdr_output), exist_ok=True, parents=True)
plink2_nonCleav_xifdr_input = os.path.join(plink2_nonCleav_xifdr_output, 'plink2_nonCleav_xiFDR_Input.csv')
plink2_nonCleav_df.to_csv(plink2_nonCleav_xifdr_input)
if not file_exists(plink2_nonCleav_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(plink2_nonCleav_xifdr_input, plink2_nonCleav_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# plink2 nonCleavable CSM to 2%PPI FDR
plink2_nonCleav_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'plink2_nonCleav', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(plink2_nonCleav_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
add_xifdr_args = ["--uniquePSMs=0", "--minPeptideLength=1", "--ppifdr=2", "--forward=" + "|".join(score_names.values())]
if not file_exists(plink2_nonCleav_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(plink2_nonCleav_xifdr_input, plink2_nonCleav_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                            additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                            )


# proteinProspector
pp_result_file = os.path.join('search_results', 'proteinprospector', 'ppts_csm_inter_2.txt')# 'ec1_dsso_CSM_2p_inter.csv')
pp_df = process_proteinprospector(pp_result_file, ecoli_proteins, full_fasta)

pp_xifdr_output = os.path.join('search_results', 'proteinprospector', 'FDR')
Path.mkdir(Path(pp_xifdr_output), exist_ok=True, parents=True)
pp_xifdr_input = os.path.join(pp_xifdr_output, 'proteinprospector_xiFDR_Input.csv')
pp_df.to_csv(pp_xifdr_input)
if not file_exists(pp_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(pp_xifdr_input, pp_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# proteinProspector CSM to 2%PPI FDR
pp_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'proteinprospector', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(pp_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(pp_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(pp_xifdr_input, pp_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                            additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                            )

# ppi data
pp_ppi_result_file = os.path.join('search_results', 'proteinprospector', 'ppts_ppi_inter_2.txt')
pp_ppi_df = process_proteinprospector(pp_ppi_result_file, ecoli_proteins, full_fasta)
print (f'pp_ppi_df.shape = {pp_ppi_df.shape}')

pp_ppi_xifdr_output = os.path.join('search_results', 'proteinprospector', 'FDR_PPI')
Path.mkdir(Path(pp_ppi_xifdr_output), exist_ok=True, parents=True)
pp_ppi_xifdr_input = os.path.join(pp_ppi_xifdr_output, 'proteinprospector_xiFDR_Input.csv')
pp_ppi_df.to_csv(pp_ppi_xifdr_input)
print (f'pp_ppi_xifdr_input = {pp_ppi_xifdr_input}')
if not file_exists(pp_ppi_xifdr_output, ".*_PPI_xiFDR.*"):
    xifdr.xifdr_execution(pp_ppi_xifdr_input, pp_ppi_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )


# xiSEARCH
#prepare from basic xiSEARCH results
xisearch_result_file = os.path.join(
    'search_results', 'xisearch',
    'ID2all+human3_linMod_filtered_0.csv')
xisearch_xifdr_output = os.path.join(os.getcwd(), 'search_results', 'xisearch', 'FDR')
Path.mkdir(Path(xisearch_xifdr_output), exist_ok=True, parents=True)
xisearch_xifdr_args = ["--psmfdr=2", "--boost=PSM", "--boost-between"]
if not file_exists(xisearch_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xisearch_result_file, xisearch_xifdr_output, pepfdr="100", memory="6G", reportfactor="10000",
                          additional_xifdr_arguments=xisearch_xifdr_args, xifdr_filename=xiFDR_jar
                          )

xisearch_xifdr_output = os.path.join('search_results', 'xisearch', 'FDR_PPI')
Path.mkdir(Path(xisearch_xifdr_output), exist_ok=True, parents=True)
xisearch_xifdr_args = [ "--ppifdr=2", "--boost=ppi", "--boost-between"]
if not file_exists(xisearch_xifdr_output, ".*_ppi_xiFDR.*"):
    xifdr.xifdr_execution(xisearch_result_file, xisearch_xifdr_output, pepfdr="100", memory="6G", reportfactor="10000",
                          additional_xifdr_arguments=xisearch_xifdr_args, xifdr_filename=xiFDR_jar
                          )


xisearch_result_file = os.path.join(
    'search_results', 'xisearch', 'FDR',
    'FDR_CSM_xiFDR2.2.csv')

#xisearch_result_file = "/data/rappstore/users/lfischer/search_comparison/results/xisearch/PPI_FDR_proteome+human3_linMod_filtered_0_PSM_xiFDR1.3.36.csv"
xisearch_df = process_xisearch(xisearch_result_file, ecoli_proteins)

#xisearch_xifdr_output = os.path.join('search_results', 'xisearch', 'FDR')
#Path.mkdir(Path(xisearch_xifdr_output), exist_ok=True, parents=True)
#xisearch_xifdr_input = os.path.join(xisearch_xifdr_output, 'xisearch_xiFDR_Input.csv')
#xisearch_df.to_csv(xisearch_xifdr_input)
#xifdr.xifdr_execution(xisearch_xifdr_input, xisearch_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
#                      additional_xifdr_arguments=add_xifdr_args, xifdr_filename=xiFDR_jar
#
#                       )

# xiSEARCH CSM to 2%PPI FDR
xisearch_xifdr_input = os.path.join(
    'search_results', 'xisearch', 'FDR', 'FDR_CSM_xiFDR2.2.csv')

xisearch_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'xisearch', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(xisearch_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(xisearch_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xisearch_xifdr_input, xisearch_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                            additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                            )

# xiSEARCH Highest to PPI
xisearch_ppi_result_file = os.path.join(
    'search_results', 'xisearch', 'FDR_PPI',
    'FDR_ppi_xiFDR2.2.csv')
xisearch_ppi_df = process_xisearch(xisearch_ppi_result_file, ecoli_proteins)
xisearch_ppi_output = os.path.join('search_results', 'xisearch', 'FDR_PPI')
Path.mkdir(Path(xisearch_ppi_output), exist_ok=True, parents=True)
xisearch_ppi_df.to_csv(os.path.join(xisearch_ppi_output, 'xisearch_ppi_FDR.csv'))


# XlinkX
xlinkx_result_file = os.path.join('search_results', 'xlinkx', '2p_all_files_plus3_CSMs.txt')
xlinkx_decoy_file = os.path.join('search_results', 'xlinkx', '2p_all_files_plus3_DecoyCSMs.txt')
xlinkx_df = process_xlinkx(xlinkx_result_file, xlinkx_decoy_file, ecoli_proteins, all_peptide_proteins,all_peptide_positions)

xlinkx_xifdr_output = os.path.join('search_results', 'xlinkx', 'FDR')
Path.mkdir(Path(xlinkx_xifdr_output), exist_ok=True, parents=True)
xlinkx_xifdr_input = os.path.join(xlinkx_xifdr_output, 'xlinkx_xiFDR_Input.csv')
xlinkx_df.to_csv(xlinkx_xifdr_input)
if not file_exists(xlinkx_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xlinkx_xifdr_input, xlinkx_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# XlinkX CSM to 2%PPI FDR
xlinkx_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'xlinkx', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(xlinkx_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(xlinkx_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xlinkx_xifdr_input, xlinkx_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                          )

# XlinkX pd3.2
xlinkx_pre_result_file = os.path.join('search_results', 'xlinkx_pd32', 'build410CSMs_all_true.xlsx')
xlinkx_pre_decoy_file = os.path.join('search_results', 'xlinkx_pd32', 'build410decoyCSMs_all_true.xlsx')
xlinkx_pre_df = process_xlinkx(xlinkx_pre_result_file, xlinkx_pre_decoy_file, ecoli_proteins, all_peptide_proteins,all_peptide_positions)

xlinkx_pre_xifdr_output = os.path.join('search_results', 'xlinkx_pd32', 'FDR')
Path.mkdir(Path(xlinkx_pre_xifdr_output), exist_ok=True, parents=True)
xlinkx_pre_xifdr_input = os.path.join(xlinkx_pre_xifdr_output, 'xlinkx_pre_xiFDR_Input.csv')
xlinkx_pre_df.to_csv(xlinkx_pre_xifdr_input)
if not file_exists(xlinkx_pre_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xlinkx_pre_xifdr_input, xlinkx_pre_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# XlinkX pd3.2
xlinkx_pre_xl_result_file = os.path.join('search_results', 'xlinkx_pd32', 'build410XL_all_true.xlsx')
xlinkx_pre_xl_decoy_file = None # = os.path.join('search_results', 'xlinkx_pd32', 'build410decoyCSMs_all_true.xlsx')
xlinkx_pre_xl_df = process_xlinkx_respair_to_xiFDR(xlinkx_pre_xl_result_file, xlinkx_pre_xl_decoy_file, ecoli_proteins, all_peptide_proteins,all_peptide_positions)

xlinkx_pre_xl_xifdr_output = os.path.join('search_results', 'xlinkx_pd32', 'FDR_xl')
Path.mkdir(Path(xlinkx_pre_xl_xifdr_output), exist_ok=True, parents=True)
xlinkx_pre_xl_xifdr_input = os.path.join(xlinkx_pre_xl_xifdr_output, 'xlinkx_pre_xl_xiFDR_Input.csv')
xlinkx_pre_xl_df.to_csv(xlinkx_pre_xl_xifdr_input)
if not file_exists(xlinkx_pre_xl_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xlinkx_pre_xl_xifdr_input, xlinkx_pre_xl_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )


# XlinkX_pre CSM to 2%PPI FDR
xlinkx_pre_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'xlinkx_pd32', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(xlinkx_pre_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(xlinkx_pre_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xlinkx_pre_xifdr_input, xlinkx_pre_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                          )

# XlinkX_pre xl to PPI FDR
xlinkx_pre_xl_xifdr_XL_TO_PPI_output = os.path.join('search_results', 'xlinkx_pd32', 'FDR_XL_TO_PPI')
Path.mkdir(Path(xlinkx_pre_xl_xifdr_XL_TO_PPI_output), exist_ok=True, parents=True)
if not file_exists(xlinkx_pre_xl_xifdr_XL_TO_PPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(xlinkx_pre_xl_xifdr_input, xlinkx_pre_xl_xifdr_XL_TO_PPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )


# MeroX
merox_result_file = os.path.join('search_results', 'merox', '2p_FDR_230113.csv')
merox_df = process_merox(merox_result_file, ecoli_proteins)

merox_xifdr_output = os.path.join('search_results', 'merox', 'FDR')
Path.mkdir(Path(merox_xifdr_output), exist_ok=True, parents=True)
merox_xifdr_input = os.path.join(merox_xifdr_output, 'merox_xiFDR_Input.csv')
merox_df.to_csv(merox_xifdr_input)
if not file_exists(merox_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(merox_xifdr_input, merox_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# MeroX CSM to 2%PPI FDR
merox_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'merox', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(merox_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(merox_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(merox_xifdr_input, merox_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                            additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                            )


# OpenPepXL
openpepxl_result_file = os.path.join('search_results', 'openpepxl', 'OpenPepXL_DSSO.csv')
openpepxl_df = process_openpepxl(openpepxl_result_file, ecoli_proteins)

openpepxl_xifdr_output = os.path.join('search_results', 'openpepxl', 'FDR')
Path.mkdir(Path(openpepxl_xifdr_output), exist_ok=True, parents=True)
openpepxl_xifdr_input = os.path.join(openpepxl_xifdr_output, 'openpepxl_xiFDR_Input.csv')
openpepxl_df.to_csv(openpepxl_xifdr_input)
if not file_exists(openpepxl_xifdr_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(openpepxl_xifdr_input, openpepxl_xifdr_output, pepfdr="100", memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=add_xifdr_args_forward, xifdr_filename=xiFDR_jar
                          )

# OpenPepXL CSM to 2%PPI FDR
openpepxl_xifdr_CSM_TO_2PercentPPI_output = os.path.join('search_results', 'openpepxl', 'FDR_CSM_TO_2PercentPPI')
Path.mkdir(Path(openpepxl_xifdr_CSM_TO_2PercentPPI_output), exist_ok=True, parents=True)
if not file_exists(openpepxl_xifdr_CSM_TO_2PercentPPI_output, ".*_CSM_xiFDR.*"):
    xifdr.xifdr_execution(openpepxl_xifdr_input, openpepxl_xifdr_CSM_TO_2PercentPPI_output, pepfdr="100", memory="1G", reportfactor="10000",
                            additional_xifdr_arguments=add_xifdr_args_ppi, xifdr_filename=xiFDR_jar
                            )


# concatenate all dataframes
df = pd.concat([kojak_df, mango_df, plink2_nonCleav_df,# plink2_cleav_df, 
                pp_df, xisearch_df, xlinkx_df, xlinkx_pre_df, merox_df, openpepxl_df],
               ignore_index=True)
df.to_csv(os.path.join('search_results', 'processed.csv'), index=False)


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


df_CSM_to_2PPI = read_all_xiFDR_result_files(file=r'.*_ppi_.*.csv', ecoli_proteins=ecoli_proteins, default_FDR_Folder='FDR_CSM_TO_2PercentPPI')
df_CSM_to_2PPI.to_csv(os.path.join('search_results', 'processed_CSM_TO_2PercentPPI_results.csv'), index=False)

df_highest_to_ppi = read_all_xiFDR_result_files(file=r'.*_ppi_.*.csv', ecoli_proteins=ecoli_proteins, default_FDR_Folder='FDR', FDR_Folders={'xisearch': 'FDR_PPI', 'kojak': 'FDR_PPI', 'mango': 'FDR_PPI', 'pp': 'FDR_PPI', 'xlinkx_pre': 'FDR_XL_TO_PPI'})    
df_highest_to_ppi.to_csv(os.path.join('search_results', 'processed_highest_to_ppi_results.csv'), index=False)

