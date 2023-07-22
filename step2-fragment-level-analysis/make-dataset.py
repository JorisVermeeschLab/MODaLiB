import time

import numpy as np
import pandas as pd
from modalib.core import Sample


BLACKLIST = {
    'GC1002975',  # Duplicate of GC110694
    'GC1002979',  # Duplicate of GC110695
    'GC114178',  # Pregnant
    'GC114150',  # Died - no cancer info
    'GC114151',  # Died - no cancer info
    'GC114152',  # Died - no cancer info
    'GC114153',  # Died - no cancer info
    'GC114154',  # Died - no cancer info
    'GC114155',  # Died - no cancer info
    'GC114156',  # Died - no cancer info
    'GC114157',  # Died - no cancer info
    'GC114158',  # Died - no cancer info
    'GC114126',  # No consent
    'GC112910',  # No consent
    'GC112908',  # No consent
    'GC112889',  # No info
    'GC112886',  # No info
    'GC112885',  # No consent
    'GC119371',  # No info
    'GC116988',  # Died - no cancer info
    'GC114161',  # No info
    'GC114149',  # carcinoma?
}
BLACKLISTED_CATEGORIES = {
    'EXTERNAL PROJECT', 'Healthy pregnant', 'Pregnant Aneuploidies', 'grp_newsample',
    'Cancer in Pregnancy', 'GRP_abn'
}

df = pd.read_csv('../samples-details.tsv', delimiter='\t')
categories = df['category'].to_numpy()
gc_codes = df['GC CODE'].to_numpy()
sample_types = df['Sample type'].to_numpy()
cancer_types = df['Cancer Type'].to_numpy()
cat = {}


selected_gc_codes, selected_categories = [], []
for gc_code, category, sample_type, cancer_type in zip(gc_codes, categories, sample_types, cancer_types):
    if sample_type.lower() == 'cfdna':
        if isinstance(gc_code, str) and (gc_code not in BLACKLIST):
            if category not in BLACKLISTED_CATEGORIES:
                #if 'GRP' in category:
                #    continue
                selected_gc_codes.append(gc_code)
                selected_categories.append(category)
                print(gc_code, category, cancer_type)
selected_categories = np.asarray(selected_categories, dtype=object)
print(len(selected_categories))

print(selected_gc_codes)

