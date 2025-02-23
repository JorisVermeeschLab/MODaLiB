import numpy as np
import pandas as pd


def is_old_sample(gc_code: str) -> bool:
    if not isinstance(gc_code, str):
        return False
    assert gc_code.startswith('GC')
    gc_code = gc_code.split('_')[0]
    if len(gc_code) == 9:
        return True
    elif len(gc_code) == 8:
        i = int(gc_code[2:])
        return i < 121000
    else:
        print(gc_code)
        assert False


df = pd.read_excel('metadata/Metadata_EMseqsamples.xlsx')

# Get breast cancer cases (BRCA)
mask = df['category'].apply(lambda s: s in {'grp_cancer', 'GRP_cancer', 'cancer', 'Cancer'})
mask &= df['Cancer Type'].apply(lambda s: s in {'breast cancer', 'Breast Cancer', 'Breast cancer'})
mask &= df['GC CODE'].apply(lambda s: is_old_sample(s))
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
print(int(np.sum(mask)))
print(list(df['GC CODE'][mask].to_numpy()))

# Get colorectal cancer cases (COLCA)
mask = df['category'].apply(lambda s: s in {'cancer', 'Cancer', 'GRP_cancer', 'grp_cancer'})
mask &= df['Cancer Type'].apply(lambda s: s in {'Colorectal cancer', 'colorectal cancer', 'Colorectal Cancer'})
mask &= df['GC CODE'].apply(lambda s: not is_old_sample(s))
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
mask &= df['other info'].apply(lambda s: 'nivolumab' not in str(s).lower())
mask &= df['other info'].apply(lambda s: str(s).lower() not in {'before therapy'})
gc_codes = list(df['GC CODE'][mask].to_numpy())
gc_codes = list(set(gc_codes).union({"GC142079"}))
gc_codes = list(set(gc_codes).difference({'GC146664_CAATAGCC-TCTACGCA'}))
print(len(gc_codes))
print(gc_codes)

# Get controls for COLCA (LVT)
mask = df['category'].apply(lambda s: s in {'Healthy controls', 'GRP_healthy'})
mask &= df['GC CODE'].apply(lambda s: not is_old_sample(s))
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
gc_codes = list(df['GC CODE'][mask].to_numpy())
gc_codes = list(set(gc_codes).difference({'GC121334', 'GC121340', 'GC145429', 'GC145430', 'GC146665_TCATCTCC-CTTGCTAG', 'GC121345', 'GC121337', 'GC121340, GC121334', 'GC121344', 'GC121341', 'GC121338', 'GC121342', 'GC121336', 'GC121331', 'GC121348', 'GC121347', 'GC121349', 'GC121333', 'GC121343', 'GC121335', 'GC121346', 'GC121332', 'GC121339'}))
gc_codes = list(set(gc_codes).difference({'GC146769_CTTGCTAG-TCATCTCC'}))  # # GC131898, GC146769_CTTGCTAG-TCATCTCC, GC146769_GTCGATTG-AAGACCGT
print(len(gc_codes))
print(gc_codes)


print(set(df['other info']))
