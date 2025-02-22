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
print(df.columns)

# Get controls for BRCA (GRP)
#print(set(df['category']))
#mask = df['category'].apply(lambda s: s in {'GEP Healthy controls', 'Healthy controls_repl', 'GEP_Healthy controls', 'Healthy controls', 'GRP_healthy'})
##mask &= df['GC CODE'].apply(lambda s: is_old_sample(s))
#mask &= df['sex'].apply(lambda s: s == 'F')
#mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
#mask &= df['status'].apply(lambda s: s == 'sequenced')
#mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
#mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
#print(int(np.sum(mask)))
#print(list(df['GC CODE'][mask].to_numpy()))

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

mask = df['category'].apply(lambda s: s in {'cancer', 'Cancer', 'GRP_cancer', 'grp_cancer'})
mask &= df['Cancer Type'].apply(lambda s: s not in {'Colorectal cancer', 'colorectal cancer', 'Colorectal Cancer', 'breast cancer', 'Breast Cancer', 'Breast cancer'})
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
print(list(df['GC CODE'][mask].to_numpy()))

# Get colorectal cancer cases (COLCA)
mask = df['category'].apply(lambda s: s in {'cancer', 'Cancer'})
mask &= df['Cancer Type'].apply(lambda s: s in {'Colorectal cancer', 'colorectal cancer', 'Colorectal Cancer'})
mask &= df['GC CODE'].apply(lambda s: not is_old_sample(s))
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
print(int(np.sum(mask)))
print(list(df['GC CODE'][mask].to_numpy()))

# Get controls for COLCA (LVT)
mask = df['category'].apply(lambda s: s in {'Healthy controls'})
mask &= df['GC CODE'].apply(lambda s: not is_old_sample(s))
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl'})
mask &= df['notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
print(int(np.sum(mask)))
print(list(df['GC CODE'][mask].to_numpy()))

print(set(df['category']))
