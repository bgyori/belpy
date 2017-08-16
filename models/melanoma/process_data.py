import pandas
from indra.databases import uniprot_client
from indra.literature.pubmed_client import get_ids_for_gene

data_fname = 'copies_per_cell.csv'

def read_data(fname=data_fname):
    data = pandas.read_csv(fname)
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
