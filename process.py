from process_results import *
from plots_and_functions import fasta_to_dict

import os
import pandas as pd


# read in E. coli fasta file
ecoli_fasta = fasta_to_dict('database/uniprot-k12-filtered-proteome_UP000000625.fasta')
# read in uniprot IDs of E. coli as list
ecoli_proteins = list(ecoli_fasta.keys())

#
# process the results
#
# Kojak
kojak_result_file = os.path.join('search_results', 'kojak', 'Kojak_CSMs_NoLoopLinks.tsv')
kojak_df = process_kojak(kojak_result_file, ecoli_proteins)

# Mango Comet XlinkProphet
mango_result_file = os.path.join('search_results', 'mango_comet_xlinkprophet',
                                 'iprophet-withdecoys-xl-filtered_2pct.xlsx')
mango_df = process_mango_comet_xlinkprophet(mango_result_file, ecoli_proteins)

# pLink2
plink2_result_file = os.path.join('search_results', 'plink2',
                                  'LK_2.3.9_cleavable_ecoli_human_2023.01.02.filtered_cross-linked_spectra.csv')
pllink2_unfiltered_result_file = os.path.join('search_results', 'plink2',
                                              'LK_2.3.9_cleavable_ecoli_human_2023.01.02.csv')
plink2_df = process_plink2(plink2_result_file, pllink2_unfiltered_result_file, ecoli_proteins)

# proteinProspector
pp_result_file = os.path.join('search_results', 'proteinprospector', 'ec1_dsso_CSM_2p_inter.csv')
pp_df = process_proteinprospector(pp_result_file, ecoli_proteins)

# xiSEARCH
xisearch_result_file = os.path.join(
        'search_results', 'xisearch',
        'xiSEARCH_CSM_FDR_proteome+human3_linMod_filtered_0_PSM_xiFDR1.3.36.csv')
xisearch_df = process_xisearch(xisearch_result_file, ecoli_proteins)

# XlinkX
xlinkx_result_file = os.path.join('search_results', 'xlinkx', '2p_all_files_plus3_CSMs.txt')
xlinkx_decoy_file = os.path.join('search_results', 'xlinkx', '2p_all_files_plus3_DecoyCSMs.txt')
xlinkx_df = process_xlinkx(xlinkx_result_file, xlinkx_decoy_file, ecoli_proteins)

# MeroX
merox_result_file = os.path.join('search_results', 'merox', '2p_FDR_230113.csv')
merox_df = process_merox(merox_result_file, ecoli_proteins)

# OpenPepXL
openpepxl_result_file = os.path.join('search_results', 'openpepxl', 'OpenPepXL_DSSO.csv')
openpepxl_df = process_openpepxl(openpepxl_result_file, ecoli_proteins)

# concatenate all dataframes
df = pd.concat([kojak_df, mango_df, plink2_df, pp_df, xisearch_df, xlinkx_df, merox_df, openpepxl_df],
               ignore_index=True)

# save to csv
df.to_csv(os.path.join('search_results', 'processed_results.csv'), index=False)
