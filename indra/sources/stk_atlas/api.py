__all__ = ['process_local', 'process_from_web', 'process_from_dataframe']

import logging
import pandas as pd
from .processor import StkAtlasProcessor

logger = logging.getLogger(__name__)

kinase_table = '41586_2022_5575_MOESM3_ESM.xlsx'
kinase_sheet_name = 'Table S1 Data'
annot_table = '41586_2022_5575_MOESM5_ESM.xlsx'
annot_sheet_name = 'Supplementary Table 3'

base_url = ('https://static-content.springer.com/esm/'
            'art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/')


# The maximum rank to consider for annotations to construct Statements
DEFAULT_MAX_RANK = 15
# The minimum percentile cutoff for annotations to construct INDRA Statements
DEFAULT_MIN_PERCENTILE = 90


def process_from_web(cache_local=False, max_rank=DEFAULT_MAX_RANK,
                     min_percentile=DEFAULT_MIN_PERCENTILE):
    kinase_table_path = base_url + kinase_table
    annotation_table_path = base_url + annot_table
    # If we use caching, we override the URLs to be local paths
    if cache_local:
        import pystow
        kinase_table_path = pystow.ensure('indra', 'stk_atlas', kinase_table,
                                          url=kinase_table_path)
        annotation_table_path = pystow.ensure('indra', 'stk_atlas',
                                              annot_table,
                                              url=annotation_table_path)
    # Now load the actual tables
    kinase_df = pd.read_excel(kinase_table_path,
                              sheet_name=kinase_sheet_name)
    # FIXME: the annotations table errors when loaded with openpyxl/pandas
    annot_df = pd.read_excel(annotation_table_path,
                             sheet_name=annot_sheet_name)
    return process_from_dataframe(kinase_df, annot_df, max_rank=max_rank,
                                  min_percentile=min_percentile)


def process_local(kinase_xlsx, annot_csv, max_rank=DEFAULT_MAX_RANK,
                  min_percentile=DEFAULT_MIN_PERCENTILE):
    kinase_df = pd.read_excel(kinase_xlsx, sheet_name=kinase_sheet_name,
                              )
    annot_df = pd.read_csv(annot_csv)
    return process_from_dataframe(kinase_df, annot_df, max_rank=max_rank,
                                  min_percentile=min_percentile)


def process_from_dataframe(kinase_df, annot_df, max_rank=DEFAULT_MAX_RANK,
                           min_percentile=DEFAULT_MIN_PERCENTILE):
    sp = StkAtlasProcessor(kinase_df, annot_df, max_rank=max_rank,
                           min_percentile=min_percentile)
    sp.process_phosphorylations()
    return sp
