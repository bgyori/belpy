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


def process_from_web(cache_local=False):
    kinase_table_url = base_url + kinase_table
    annotation_table_url = base_url + annot_table
    if cache_local:
        import pystow
        kinase_table_local = pystow.ensure('indra', 'stk_atlas', kinase_table,
                                           url=kinase_table_url)
        annotation_table_local = pystow.ensure('indra', 'stk_atlas',
                                               annot_table,
                                               url=annotation_table_url)
        kinase_df = pd.read_excel(kinase_table_local,
                                  sheet_name=kinase_sheet_name)
        annot_df = pd.read_excel(annotation_table_local,
                                 sheet_name=annot_sheet_name)
    else:
        kinase_df = pd.read_excel(kinase_table_url,
                                  sheet_name=kinase_sheet_name)
        annot_df = pd.read_excel(annotation_table_url,
                                 sheet_name=annot_sheet_name)


def process_local(kinase_xlsx, annot_csv):
    kinase_df = pd.read_excel(kinase_xlsx, sheet_name=kinase_sheet_name)
    annot_df = pd.read_csv(annot_csv)
    return process_from_dataframe(kinase_df, annot_df)


def process_from_dataframe(kinase_df, annot_df):
    sp = StkAtlasProcessor(kinase_df, annot_df)
    sp.process_phosphorylations()
    return sp
