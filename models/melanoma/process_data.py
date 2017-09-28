import pandas as pd
import numpy as np
from indra.databases import uniprot_client
from indra.literature.pubmed_client import get_ids_for_gene
import collections
import json
import itertools

data_fname = 'copies_per_cell.csv'

data = read_data()
cell_lines = get_cell_lines(data)


def read_data(fname=data_fname):
    data = pd.read_csv(fname)
    return data

def get_cell_lines(data):
    cell_lines = list(set([col_name.split('_')[0] for col_name
                           in data.columns[3:] if 'Bridge' not in col_name]))
    return cell_lines

def get_gene_names(data):
    up_ids = data.iloc[:,0]
    gene_names = []
    for up_id in up_ids:
        up_id = up_id.split('-')[0]
        gene_name = uniprot_client.get_gene_name(up_id)
        if gene_name is not None:
            gene_names.append(gene_name)
    return sorted(list(set(gene_names)))

def get_pmids(genes):
    pmids = []
    for gene in genes:
        try:
            gene_pmids = get_ids_for_gene(gene)
        except ValueError:
            print('Invalid gene symbol: %s' % gene)
            continue
        pmids += gene_pmids
    return sorted(list(set(pmids)))


def get_data_fold_changes(data, intra_cl_FC_cutoff=0.3,
                          lower_limit_FC_cutoff=-2):
    """
    Processes the data from file to combine replicates by averaging and
    return fold changes for each cell_line, gene combo
    intra_cl_FC_cutoff : int or float
        Fold change cutoff to compare for any two replicates. Sets the value to
        NaN if the FC difference exceeds the cutoff.
    lower_limit_FC_cutoff : int or float
        Lower limit for FC over the whole dataset. In order to not exaggerate
        the significance of a large reduction in values (which is much less
        likely to be detected), fold changes cannot dip below this value
    """

    data_fc = pd.DataFrame(columns=(['Gene_Symbol'] + cell_lines))
    data_fc['Gene_Symbol'] = data['Gene_Symbol']
    for line in cell_lines:
        line_cols = data.columns.tolist()
        line_cols = [x for x in line_cols if line in x]
        s1 = data[line + '_Tr1'] / data[line + '_Cr1']
        s2 = data[line + '_Tr2'] / data[line + '_Cr2']
        f = abs(s1 - s2)
        s = (s1 + s2) / 2
        # this makes all intra cell line variability NaN if it is greater than FC 1
        s[np.log10(f) > intra_cl_FC_cutoff] = np.nan
        # don't want fold changes way below detection limit throwing us off
        s[s < lower_limit_FC_cutoff] = lower_limit_FC_cutoff
        s[np.isinf(s)] = np.nan
        data_fc[line] = np.log10(s)
    return data_fc


def build_extremes_table(data_fc, extreme_limit=0.4):
    """
    This function builds a csv table with headers (human browsable)
    and one without headers for all extreme perturbations
    """
    liml = -extreme_limit
    limu = extreme_limit
    extremes = []
    for cl1, cl2 in itertools.combinations(cell_lines, 2):
        filt = ( ( (data_fc[cl1] < liml) & (data_fc[cl2] > limu) ) |
                 ( (data_fc[cl1] > limu) & (data_fc[cl2] < liml) ) )
        data_fc_f = data_fc[filt]
        data_fc_f = data_fc_f[['Gene_Symbol', cl1, cl2]]
        extreme_genes = data_fc_f['Gene_Symbol'].tolist()
        difference = (data_fc_f[cl1] - data_fc_f[cl2]).tolist()
        for ex_g, d in zip(extreme_genes, difference):
            extremes.append((cl1, cl2, ex_g, d))
    extremes = sorted(extremes, key=lambda x: abs(x[3]), reverse=True)
    print(extremes)
    df_extremes = pd.DataFrame(columns=['Cell_line1', 'Cell_line2',
                                        'Gene_Symbol', 'Difference'],
                               data = sorted(extremes, key=lambda x: abs(x[3]),
                                             reverse=True))
    df_extremes.to_csv('extremes.csv', index=False)
    df_extremes.to_csv('extremes_nohead.csv', index=False, header=None)
    return df_extremes


def build_extremes_each_line(data_fc, extreme_limit):
    """
    This function builds a csv table with headers (human browsable)
    and one without headers for all extreme perturbations
    """
    liml = -extreme_limit
    limu = extreme_limit
    extremes = []
    for cl in cell_lines:
        filt = ( ( (data_fc[cl] < liml) | (data_fc[cl] > limu) ) )
        data_fc_f = data_fc[filt]
        data_fc_f = data_fc_f[['Gene_Symbol', cl]]
        extreme_genes = data_fc_f['Gene_Symbol'].tolist()
        values = data_fc_f[cl].tolist()
        for ex_g, d in zip(extreme_genes, values):
            extremes.append((cl, ex_g, d))
    extremes = sorted(extremes, key=lambda x: abs(x[2]), reverse=True)
    print(extremes)
    df_extremes = pd.DataFrame(columns=['Cell_line',
                                        'Gene_Symbol', 'Difference'],
                               data = sorted(extremes, key=lambda x: abs(x[2]),
                                             reverse=True))
    df_extremes.to_csv('extremes_each_cl.csv', index=False)
    df_extremes.to_csv('extremes_each_cl_nohead.csv', index=False, header=None)
    return df_extremes




def build_context_json_for_cyjs(data_fc, gene_list):
    """
    Make a (human readable) melted table of the fold change data and a json
    file for cyjs context setting limited to the gene_list provided
    """
    data_fc = data_fc[data_fc['Gene_Symbol'].isin(gene_names)]
    data_fc = data_fc.set_index('Gene_Symbol')
    data_fc = data_fc.transpose()
    data_fc = data_fc.reset_index(drop=False)
    data_fc = data_fc.melt(id_vars=['index'])
    data_fc.index.name = 'id'
    data_fc.reset_index(inplace=True)
    data_fc.drop(['id'], axis=1, inplace=True)
    data_fc.rename(columns={'index': 'id'}, inplace=True)
    data_fc.rename(columns={'Gene_Symbol': 'gene'}, inplace=True)
    data_fc['antibody'] = [None] * len(data_fc)
    data_fc['site'] = [None] * len(data_fc)
    recursivedict = lambda: collections.defaultdict(recursivedict)
    j = recursivedict()
    df1 = data_fc
    ids = df1['id'].unique().tolist()
    for i in ids:
        df2 = df1[df1['id'] == i]
        genes = df2['gene'].unique().tolist()
        for g in genes:
            df3 = df2[df2['gene'] == g]
            df4 = df3[~df3['site'].isin([None])]
            j[i][g]['members']['antibodies'] = df4['antibody'].tolist()
            j[i][g]['members']['sites'] = df4['site'].tolist()
            j[i][g]['members']['values'] = df4['value'].tolist()
            df5 = df3[df3['site'].isin([None])]
            j[i][g]['node']['antibodies'] = df5['antibody'].tolist()
            j[i][g]['node']['values'] = df5['value'].tolist()
    df1.to_csv('data_fc_melted.csv', index=False)
    with open('data_fc.json', 'w') as outfile:
        json.dump(j, outfile, indent=2, sort_keys=True)
