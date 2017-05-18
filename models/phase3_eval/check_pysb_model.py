import pickle
from indra.tools.model_checker import ModelChecker
#from manual_stmts import stmts as manual_stmts
from assemble_pysb import set_context, add_observables
import process_data
from indra.util import write_unicode_csv
from indra.assemblers import PysbAssembler
import make_stmts_for_checking as make_stmts

print("Processing data")

data = process_data.read_data(process_data.data_file)
data_genes = process_data.get_all_gene_names(data)


print('Loading data statements.')
data_stmts, data_values = make_stmts.run(dec_thresh=0.5, inc_thresh=1.5)

with open('korkut_stmts_no_ev.pkl', 'rb') as f:
    print('Loading korkut_model_pysb statements.')
    base_stmts = pickle.load(f)

# Merge the sources of statements
# stmts = manual_stmts + base_stmts
stmts = base_stmts
#stmts = manual_stmts

# Assemble model
pa = PysbAssembler()
pa.add_statements(stmts)
model = pa.make_model()

#with open('korkut_pysb.pkl', 'wb') as f:
#    pickle.dump(pa.model, f)

# Preprocess and assemble the pysb model
#model = assemble_pysb(combined_stmts, data_genes, '')

mc = ModelChecker(model)

# Iterate over each drug/ab statement subset
results = []
for drug_name, ab_dict in data_stmts.items():
    for ab, stmt_list in ab_dict.items():
        value = data_values[drug_name][ab]
        # For each subset, check statements; if any of them checks out, we're
        # good and can move on to the next group
        print("-- Checking the effect of %s on %s --" % (drug_name, ab))
        relation = 'positive' if value > 1 else 'negative'
        path_found = 0
        path = ''
        for stmt in stmt_list:
            print("Checking: %s" % stmt)
            result = mc.check_statement(stmt)
            if result:
                print("Path found, skipping rest")
                path_found = 1
                path = str(result)
                break
            else:
                print("No path found")

        results.append((drug_name, ab, relation, value, path_found, path))
write_unicode_csv('model_check_results.csv', results)