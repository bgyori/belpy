import os
import sys
import process_data
from random import shuffle
from indra.sources import biopax
import indra.tools.assemble_corpus as ac
from indra.assemblers import CxAssembler
from indra.literature.pubmed_client import get_ids_for_gene
from indra.statements import *
from indra.assemblers import CyJSAssembler
from tqdm import tqdm
from indra.tools.reading.submit_reading_pipeline_aws import \
    submit_run_reach, wait_for_complete
from indra.sources import biopax, bel
import pandas as pd

if sys.version_info[0] < 3:
    raise Exception('Run this only in Python 3')

def read_gene_list(path):
    df = pd.read_csv(path, sep='\t', error_bad_lines=False, header=None)
    gene_list = df[0].tolist()
    return sorted(list(set(gene_list)))

def get_biopax_stmts(gene_list):
    bp_stmts_path = 'stmts_from_dbs/bp_stmts.pkl'
    if os.path.isfile(bp_stmts_path):
        bp_stmts = ac.load_statements(bp_stmts_path)
    else:
        bp_stmts = []
        for gene in tqdm(gene_list):
            bp = biopax.process_pc_neighborhood([gene])
            bp_stmts += bp.statements
        ac.dump_statements(bp_stmts, bp_stmts_path)
    return bp_stmts

def get_ndex_stmts(gene_list):
    ndex_stmts_path = 'stmts_from_dbs/ndex_stmts.pkl'
    if os.path.isfile(ndex_stmts_path):
        ndex_stmts = ac.load_statements(ndex_stmts_path)
    else:
        ndex_stmts = []
        for gene in tqdm(gene_list):
            ndex = bel.process_ndex_neighborhood([gene])
            ndex_stmts += ndex.statements
        ac.dump_statements(ndex_stmts, ndex_stmts_path)
    return ndex_stmts

def chunked_assembly(chunk_stmts, save_file, gene_names):
    chunk_stmts = ac.map_grounding(chunk_stmts)
    chunk_stmts = ac.filter_grounded_only(chunk_stmts)
    chunk_stmts = ac.filter_human_only(chunk_stmts)
    chunk_stmts = ac.expand_families(chunk_stmts)
    chunk_stmts = ac.filter_gene_list(chunk_stmts, gene_names, 'all')
    chunk_stmts = ac.map_sequence(chunk_stmts)
    ac.dump_statements(chunk_stmts, save_file)
    return chunk_stmts

def run_reading(pmids):
    pmid_fname = 'melanoma_pmids.txt'
    with open(pmid_fname, 'wt') as fh:
        shuffle(pmids)
        for pmid in pmids:
            fh.write('%s\n' % pmid)
    # Submit reading
    job_list = submit_run_reach('melanoma', pmid_fname)
    # Wait for reading to complete
    reading_res = wait_for_complete(job_list)

def run_assembly(stmts, save_file, gene_names):
    stmts_assembled = []
    chunk_size = 25000
    stmt_chunk_paths = []
    for i in tqdm(range(0, len(stmts), chunk_size)):
        path = 'chunks/final_stmts_' + str(i) + '.pkl'
        stmt_chunk_paths.append(path)
        if os.path.isfile(path):
            print('Chunk ' + str(i) + ' already processed')
        else:
            stmts_chunk = stmts[i:i + chunk_size]
            chunked_assembly(stmts_chunk, path, gene_names)
    stmts = []
    for p in stmt_chunk_paths:
        s = ac.load_statements(p)
        stmts += s
    stmts = ac.run_preassembly(stmts, return_toplevel=False)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_enzyme_kinase(stmts)
    ac.dump_statements(stmts, 'final_stmts.pkl')
    return stmts

def assemble_cx(stmts, save_file):
    cxa = CxAssembler(stmts)
    cxa.make_model(add_indra_json=False)
    cxa.save_model(save_file)
    return cxa

def assemble_cyjs(stmts, save_file):
    cja = CyJSAssembler()
    cja.add_statements(stmts)
    cja.make_model(grouping=True)
    cja.save_json(save_file)
    return cja

if __name__ == '__main__':
    run_reading = False
    data = process_data.read_data()
    gene_names = process_data.get_gene_names(data)
    if run_reading:
        pmids = process_data.get_pmids(gene_names)
        run_reading(pmids)
    stmts = ac.load_statements('stmts.pkl')
    ras_gene_names = get_gene_names('../../data/ras_pathway_proteins.csv')
    msb2015_gene_names = get_gene_names('MohammadFS_MSB_2015_gene_list.csv')
    msb2017_gene_names = get_gene_names('MohammadFS_MSB_2017_gene_list.csv')
    gene_list = ras_gene_names + msb2015_gene_names + msb2017_gene_names
    gene_list = sorted(list(set(gene_list)))
    stmts = run_assembly(stmts, 'final_stmts.pkl', gene_list)
    assemble_cyjs(stmts, 'model')
