from copy import deepcopy

import tqdm

from indra.databases import hgnc_client
from indra.statements import Phosphorylation
from indra.ontology.standardize import get_standard_agent


class StkAtlasProcessor:
    def __init__(self, kinase_df, annot_df, max_rank, min_percentile):
        self.kinase_df = kinase_df
        self.annot_df = annot_df
        self.max_rank = max_rank
        self.min_percentile = min_percentile
        self.kinase_lookup = self.make_kinase_lookup()
        self.statements = []

    def make_kinase_lookup(self):
        matrix_to_agent = {}
        for matrix_name, gene_name in zip(self.kinase_df.Matrix_name,
                                          self.kinase_df.Gene):
            hgnc_id = hgnc_client.get_hgnc_id(gene_name)
            assert hgnc_id, gene_name
            matrix_to_agent[matrix_name] = get_standard_agent(gene_name,
                                                              {'HGNC': hgnc_id})
        return matrix_to_agent

    def process_phosphorylations(self):
        for _, row in tqdm.tqdm(self.annot_df.iterrows()):
            substrate_uniprot = row['Uniprot Primary Accession']
            substrate_agent = get_standard_agent(substrate_uniprot,
                                                 {'UP': substrate_uniprot})
            site_str = row['Phosphosite']
            residue, position = site_str[0], site_str[1:]
            # NOTE: Database Uniprot Accession might be relevant too
            rank_rows = row.filter(regex='^.*_rank')
            top_ranks = rank_rows[rank_rows <= self.max_rank]
            for kinase_col, rank in top_ranks.items():
                kinase_matrix_name = kinase_col.split('_')[0]
                if row[f'{kinase_matrix_name}_percentile'] < self.min_percentile:
                    continue
                kinase_agent = \
                    deepcopy(self.kinase_lookup[kinase_matrix_name])
                stmt = Phosphorylation(kinase_agent, substrate_agent,
                                       residue, position)
                self.statements.append(stmt)
