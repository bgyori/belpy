import os
import pickle
from indra.statements import *
import indra.tools.assemble_corpus as ac
from indra.assemblers import PysbAssembler
from indra.explanation.model_checker import ModelChecker
from process_data import *

def make_stmts_to_check():
    # Vemurafenib treatment results in decreased amount of SPRY1, 2 and 4
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    sprys = [Agent('SPRY1', db_refs={'HGNC': '11269', 'UP': 'O43609'}),
             Agent('SPRY2', db_refs={'HGNC': '11270', 'UP': 'O43957'}),
             Agent('SPRY4', db_refs={'HGNC': '15533', 'UP': 'Q9C004'})]
    stmts = [IncreaseAmount(braf, spry) for spry in sprys]

    # Vemurafenib treatment results in decreased MITF amount
    mitf = Agent('MITF', db_refs={'HGNC': '7105', 'UP': 'O75030'})
    stmt = IncreaseAmount(braf, mitf)
    stmts.append(stmt)

    # Vemurafenib treatment results in increased CDKN1B
    cdkn1b = Agent('CDKN1B', db_refs={'HGNC': '1785', 'UP': 'P46527'})
    stmt = DecreaseAmount(braf, cdkn1b)
    stmts.append(stmt)
    return stmts

def export_paths(paths, model, stmts):
    """Export paths for pathway map in JSON-like format."""
    conc = 1
    time = 1
    cell_line = 'MMACSF'
    drug = 'Vemurafenib'
    label = '%s_%s_%s_%s' % (drug, time, conc, cell_line)
    res = {label: []}
    for path in paths:
        entry = {'meta': [], 'path': []}
        path_stmts = stmts_for_path(path, model, stmts)
        agent_names = [stmt.agent_list()[0].name for stmt in path_stmts]
        start_gene = agent_names[0]
        end_gene = path_stmts[-1].agent_list()[-1].name
        agent_names.append(end_gene)
        agent_names_str = ' > '.join(agent_names)
        length = len(path_stmts)
        score = 0.0
        uuids = [stmt.uuid for stmt in path_stmts]
        entry['path'] = uuids
        entry['meta']= [agent_names_str, start_gene, end_gene, length, score]
        res[label].append(entry)
    return res

def stmts_for_path(path, model, stmts):
    path_stmts = []
    for path_rule, sign in path:
        for rule in model.rules:
            if rule.name == path_rule:
                stmt = _stmt_from_rule(model, path_rule, stmts)
                assert stmt is not None
                path_stmts.append(stmt)
    return path_stmts

def _stmt_from_rule(model, rule_name, stmts):
    """Return the INDRA Statement corresponding to a given rule by name."""
    stmt_uuid = None
    for ann in model.annotations:
        if ann.predicate == 'from_indra_statement':
            if ann.subject == rule_name:
                stmt_uuid = ann.object
                break
    if stmt_uuid:
        for stmt in stmts:
            if stmt.uuid == stmt_uuid:
                return stmt

if __name__ == '__main__':
    based = '.' # SET BASE FOLDER HERE TO ACCESS PICKLES
    stmts_to_check = make_stmts_to_check()
    stmts = ac.load_statements(os.path.join(based, 'fallahi_eval_pysb_stmts.pkl'))
    with open(os.path.join(based, 'fallahi_eval_pysb_model.pkl'), 'rb') as fh:
        model = pickle.load(fh)

    mc = ModelChecker(model, stmts_to_check)
    mc.prune_influence_map()
    res = mc.check_model(100, 10)
    # Get the paths for CDKN1B
    cdkn1b_paths = res[4][1].paths
    path_export = export_paths(cdkn1b_paths, model, stmts)
