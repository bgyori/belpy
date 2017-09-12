import pandas as pd
from indra.databases import uniprot_client
from indra.literature.pubmed_client import get_ids_for_gene

data_fname = 'copies_per_cell.csv'

def read_data(fname=data_fname):
    data = pd.read_csv(fname)
    return data

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


def get_data_fold_changes(data, intra_cl_FC_cutoff, lower_limit_FC_cutoff):
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

    cell_lines = list(set([x.split('_')[0] for x in cols
                           if x not in ['Unnamed: 0',
                                        'Uniprot_Id',
                                        'Gene_Symbol']]))
    data.drop(['Unnamed: 0', 'Uniprot_Id'], axis=1, inplace=True)
    cell_lines = [x for x in cell_lines if 'Bridge' not in x]
    retain_cols = []
    for c in cols:
        for x in cell_lines:
            if x in c:
                retain_cols.append(c)
    retain_cols = ['Gene_Symbol'] + retain_cols
    data = data[retain_cols]
    data_fc = pd.DataFrame(columns=(['Gene_Symbol']+cell_lines))
    data_fc['Gene_Symbol'] = df['Gene_Symbol']
    for line in cell_lines:
        line_cols = df.columns.tolist()
        line_cols = [x for x in line_cols if line in x]
        s1 = np.log2(df[line + '_Tr1'] / df[line + '_Cr1'])
        s2 = np.log2(df[line + '_Tr2'] / df[line + '_Cr2'])
        f = abs(s1 - s2)
        s = (s1 + s2)/2
        # this makes all intra cell line variability NaN
        # if abs(log2(FC)) > 1
        s[f > intra_cl_FC_cutoff] = float('NaN')
        # this makes all things below -2 equal to -2
        # don't want fold changes way below detection to be exaggerated
        s[s < lower_limit_FC_cutoff] = lower_limit_FC_cutoff
        data_fc[line] = s
    return data_fc


def build_extremes_table(data_fc, extreme_limit):
    """
    This function builds a csv table with headers (human browsable)
    and one without headers for all extreme perturbations
    """
    liml = -extreme_limit
    limu = extreme_limit
    extremes = []
    for cl1, cl2 in itertools.combinations(cell_lines, 2):
        filt = (((data_fc[cl1] < liml) & (data_fc[cl2] > limu)) |
                ((data_fc[cl1] > limu) & (data_fc[cl2] < liml)))
        data_fc_f = data_fc[filt]
        data_fc_f = data_fc_f[['Gene_Symbol', cl1, cl2]]
        extreme_genes = data_fc_f['Gene_Symbol'].tolist()
        difference = (data_fc_f[cl1] - data_fc_f[cl2]).tolist()
        for ex_g, d in zip(extreme_genes, difference):
            extremes.append((cl1, cl2, ex_g, d))
    extremes = sorted(extremes, key=lambda x: abs(x[3]), reverse=True)
    df_extremes = pd.DataFrame(columns=['Cell_line1', 'Cell_line2',
                                        'Gene_Symbol', 'Difference'],
                               data=sorted(extremes, key=lambda x: abs(x[3]),
                                           reverse=True))
    df_extremes.to_csv('extremes.csv',
                       index=False)
    df_extremes.to_csv('extremes_nohead.csv',
                       index=False, header=None)
    return df_extremes


