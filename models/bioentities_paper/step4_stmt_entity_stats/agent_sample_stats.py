import os
import sys
import pickle
import numpy as np
import pandas as pd
from collections import defaultdict
from matplotlib import pyplot as plt
from indra.util import write_unicode_csv, plot_formatting as pf
from indra.tools.reading.reading_results_stats import report_grounding


def get_ungrounded_stats():
    fname = '../step3_sample_training_test/bioentities_test_stmts_mapped.pkl'
    with open(fname, 'rb') as fh:
        stmts = pickle.load(fh)
    allu_test, anyu_test = report_grounding(stmts)

    fname = '../step3_sample_training_test/training_pmid_stmts.pkl'
    stmts = []
    with open(fname, 'rb') as fh:
        st = pickle.load(fh)
        for k, v in st.items():
            stmts += v
    allu_train, anyu_train = report_grounding(stmts)
    return (allu_test, anyu_test, allu_train, anyu_train)


def plot_ungrounded_stats(allu_test, anyu_test, allu_train, anyu_train):
    pf.set_fig_params()
    plt.figure(figsize=(3, 2), dpi=300)
    xticks = np.array([0, 1])
    col_width = 0.35
    btrain = plt.bar(xticks - 0.5*col_width, [allu_train, anyu_train], col_width,
                    align='center', linewidth=0.5, color='r')
    btest = plt.bar(xticks + 0.5*col_width, [allu_test, anyu_test], col_width,
                    align='center', linewidth=0.5, color='b')
    plt.xticks(xticks, ('All args ungrounded', 'Any args ungrounded'))
    plt.ylabel('Pct. Extracted Events')
    plt.ylim((0, 35))
    ax = plt.gca()
    pf.format_axis(ax)
    plt.legend((btrain, btest), ('Training corpus', 'Test corpus'),
               loc='upper left', frameon=False, fontsize=pf.fontsize)
    plt.savefig('ungrounded_stats.pdf')


def grounding_stats(data):
    cats = (['P'], ['F', 'C', 'X'], ['S'], ['B'], ['U'], ['M'])
    cat_names = ('Protein/gene', 'Family/complex', 'Small molecule',
                 'Biological process', 'Other/unknown', 'microRNA')
    rows = []
    num_agents = len(data)
    plt.figure(figsize=(2, 2), dpi=300)
    for ix, cat in enumerate(cats):
        cat_rows = data[data.EntityType.apply(lambda et: et in cat)]
        cat_number = len(cat_rows)
        cat_pct = (100 * cat_number / float(num_agents))
        cat_pct_str = '%.1f' % cat_pct
        correct_rows = cat_rows[cat_rows.Grounding == 1]
        correct_number = len(correct_rows)
        correct_pct = (100 * correct_number / float(cat_number)) if \
            cat_number > 0 else 0
        correct_pct_of_total = (100 * correct_number) / float(num_agents)
        correct_pct_str = '%.1f' % correct_pct
        rows.append((cat, cat_number, cat_pct_str, correct_number,
                     correct_pct_str))
        def stderr(k, n):
            return np.sqrt(((k/float(n)) * (1-(k/float(n)))) / float(n))
        stderr_inc = 100 * stderr(cat_number - correct_number, num_agents)
        inc_handle = plt.bar(ix, cat_pct, color='red', align='center',
                             yerr=stderr_inc, linewidth=0.5)
        stderr_corr = 100 * stderr(correct_number, num_agents)
        corr_handle = plt.bar(ix, correct_pct_of_total, color='blue',
                              align='center', yerr=stderr_corr,
                              linewidth=0.5)
    plt.xticks(range(len(cats)), cat_names, rotation=90)
    plt.ylabel('Pct. Curated Entities')
    plt.subplots_adjust(left=0.18, bottom=0.43, top=0.96)
    ax = plt.gca()
    pf.format_axis(ax)
    plt.legend((corr_handle, inc_handle), ('Correct', 'Incorrect'),
               loc='upper right', frameon=False, fontsize=pf.fontsize)
    plt.show()
    print(rows)
    write_unicode_csv('%s_agents_sample_stats.csv' % mode, rows)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()
    mode = sys.argv[1]
    if mode not in ('training', 'test'):
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()


    pf.set_fig_params()
    plt.ion()
    family_cats = ('F', 'C', 'X')

    fname = 'training_agents_sample_curated.csv' if mode == 'training' else \
            'test_agents_with_be_sample_curated.csv'
    with open(fname, 'rb') as f:
        data = pd.read_csv(f)

    num_curated = 300
    curated = data[0:num_curated]
    grounding_stats(curated)
