from random import shuffle
from indra.tools.reading.submit_reading_pipeline_aws import \
    submit_run_reach, wait_for_complete
import indra.tools.assemble_corpus as ac
import process_data


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

def run_assembly(stmts):
    stmts = ac.map_grounding(stmts)
    # Further assembly steps here
    return stmts

if __name__ == '__main__':
    data = process_data.read_data()
    gene_names = process_data.get_gene_names(data)
    pmids = process_data.get_pmids(gene_names)
    run_reading(pmids)
    # Here the pkl needs to be combined and copied into the folder
    stmts = ac.load_statements('stmts.pkl')
    stmts = run_assembly(stmts)
