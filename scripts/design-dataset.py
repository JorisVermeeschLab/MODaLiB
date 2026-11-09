from typing import List, Tuple

import numpy as np
import datetime
import collections
import pandas as pd
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedKFold
import lightgbm
from scipy.optimize import linear_sum_assignment


RUN_TO_DATETIME = {
    'Emseq_Nov21': datetime.datetime(2021, 11, 1),
    'Emseq_Jun23': datetime.datetime(2023, 6, 1),
    'Emseq_Jan23 ': datetime.datetime(2023, 1, 1),
    'EMBCseq_Jan2022': datetime.datetime(2022, 1, 1),
    'Emseq_Jan23': datetime.datetime(2023, 1, 1),
    '230928_EMSeq_September2023': datetime.datetime(2023, 9, 28),
    'Emseq_February2024_pt2': datetime.datetime(2024, 2, 1),
    'Emseq_February2024': datetime.datetime(2024, 2, 1),
    'Emseq_August2025': datetime.datetime(2025, 8, 1),
    '200812_EMseq_methyl_BCCC': datetime.datetime(2020, 8, 12),
    'EMseq_Aug2021': datetime.datetime(2021, 8, 1),
    'EMseq_May2021': datetime.datetime(2021, 5, 1),
    'EMseq_Nov2020': datetime.datetime(2020, 11, 1),
    'Emseq_Dec22': datetime.datetime(2022, 12, 1),
    'EMseq_July23': datetime.datetime(2023, 7, 1),
    'Emseq_July21': datetime.datetime(2021, 7, 1),
    'EMseqPILOT': np.nan,
    'EMseq_Nov2020_PT2': datetime.datetime(2020, 11, 1),
    'Emseq_Nov21_part2': datetime.datetime(2021, 11, 1),
    'Emseq_February2024_pt3': datetime.datetime(2024, 2, 1),
    'Emseq_March2022': datetime.datetime(2022, 3, 1),
    'Emseq_Feb2022': datetime.datetime(2022, 2, 1),
    'Emseq_August23': datetime.datetime(2023, 8, 1),
    'Emseq_Nov22': datetime.datetime(2022, 11, 1),
    'Emseq_February2024_pt4': datetime.datetime(2024, 2, 1),
    'Emseq_Oct21_spikein': datetime.datetime(2021, 10, 1),
}
for k, v in list(RUN_TO_DATETIME.items()):
    # Convert datetime to number of years since 2020-01-01 (float)
    if isinstance(v, datetime.datetime):
        RUN_TO_DATETIME[k] = v.timestamp() / (365.25 * 24 * 60 * 60) - datetime.datetime(2020, 1, 1).timestamp() / (365.25 * 24 * 60 * 60)

SEX_TO_FLOAT = {
    'M': 0.0,
    'F': 1.0,
}

AGE_TO_FLOAT = {
    '18-30': 24.0,
    '21-30': 25.5,
    '31-40': 35.5,
    '41-50': 45.5,
    '51-60': 55.5,
    '61-70': 65.5,
    '?': np.nan,
}

STAGE_TO_FLOAT = {
    0: 0, 'Stage III': 3, 'stage I': 1, 'IVC': 4, 'IVB': 4, 'Figo IIa': 2,
    'IIB': 2, 'FIGO rIB2': 1, 'III/IV': 3,
    'III': 3, 'Figo Ib': 1, 'I?': 1, 'IIIC1': 3, 'Figo I': 1, 'IVA': 4,
    'FIGO IIIA2': 3, 'IIC': 2, 'IIIB?': 3, 'II': 2, 'IIIC2': 3, 'IIIc': 3,
    'Figo IV': 4, 'IIIB': 3, 'IIIb': 3, 'Stage I': 1, 'IV (metastasis to liver)': 4,
    'IIA': 2, 'stage M5': 4, 'IC': 1, 'IB': 1, 'IV': 4, 'I': 1, 'Stage II': 2,
    'IA': 1, 'IIIa': 3, '?': np.nan, 'IIIA': 3, 'IIIC': 3, 'Figo IIb': 2
}

SUBTYPE_MSI = {
    'MSS': (np.nan, 'MSS'),
    'MSS, ascendens': ('Ascendens', 'MSS'),
    'ascendens': ('Ascendens', np.nan),
    'MSS, caecum': ('Caecum', 'MSS'),
    'MSI, hepatic flexure': ('Liver flexure', 'MSI-High'),
    'hepatic flexure': ('Liver flexure', np.nan),
    'MSS, rectum': ('Rectum', 'MSS'),
    'rectum adenocarcinoma cancer, Liver metastasis (MSS)': ('Rectum', 'MSS'),
}

STANDARDIZE_STAGE = {
    0: '0', 'Stage III': 'III', 'stage I': 'I', 'Figo IIa': 'IIA',
    'FIGO rIB2': 'IB', 'III/IV': 'III',
    'Figo Ib': 'IB', 'I?': 'I', 'IIIC1': 'IIIC', 'Figo I': 'I',
    'FIGO IIIA2': 'IIIA', 'IIIB?': 'IIIB', 'IIIC2': 'IIIC', 'IIIc': 'IIIC',
    'Figo IV': 'IV', 'IIIb': 'IIIB', 'Stage I': 'I', 'IV (metastasis to liver)': 'IV',
    'stage M5': 'IV', 'Stage II': 'II', 'IIIa': 'IIIA', '?': np.nan, 'Figo IIb': 'IIB'
}

STANDARDIZE_MSI_STATUS = {
    'MSS?': 'MSS', 'Sporadisch MSI-H': 'MSI-High', '?': np.nan
}

STANDARDIZE_LOCATION = {
    'Liver flexure': 'Hepatic flexure', 'rectum': 'Rectosigmoid', 'no ': np.nan,
    'splenic flexure':'Splenic flexure', 'Liver flexure + Transversum': 'Transversum',
    'Liver flexure + rectosigmoïd': 'Rectosigmoid', 'rectosigmoïd': 'Rectosigmoid'
}


def age_to_float(age_str: str) -> float:
    if isinstance(age_str, int):
        return float(age_str)
    if isinstance(age_str, float):
        return age_str
    return AGE_TO_FLOAT.get(age_str, np.nan)


def match(df_cases, df_controls, vars_cols: List[str], id_col: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match cases and controls. Controls are assumed to be more numerous than cases.
    """
    
    # Convert dataframes to numpy arrays for matching
    X_cases = df_cases[vars_cols].to_numpy()
    X_controls = df_controls[vars_cols].to_numpy()
    X = np.concatenate([X_cases, X_controls], axis=0)
    y = np.concatenate([np.ones(len(X_cases)), np.zeros(len(X_controls))], axis=0)

    # Get the number of cases and controls
    n = len(df_cases)
    m = len(df_controls)

    # Produce propensity scores by cross-validation of a lightgbm model
    cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=42
    )
    model = lightgbm.LGBMClassifier()
    propensity_scores = cross_val_predict(model, X, y, cv=cv, method='predict_proba')[:, 1]

    # Propensity score difference
    ps_cost = np.abs(
        propensity_scores[:n, None] -
        propensity_scores[None, n:]
    )

    # Pairwise age difference, normalized to 10 years
    age_cases = df_cases['AGE_num'].to_numpy(dtype=float)
    age_controls = df_controls['AGE_num'].to_numpy(dtype=float)

    age_cost = np.abs(age_cases[:, None] - age_controls[None, :]) / 10.0

    # Missing age: assign a penalty
    age_cost = np.where(
        np.isnan(age_cost),
        0.5,
        age_cost
    )

    # Combined matching cost
    cost_matrix = ps_cost + 1.0 * age_cost

    # Solve the linear sum assignment problem
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Create new dataframes for matched cases and controls
    matched_cases = df_cases.iloc[row_ind]
    matched_controls = df_controls.iloc[col_ind]

    return matched_cases, matched_controls


# Load metadata
df = pd.read_excel('../data/Metadata_EMseqsamples.xlsx')

# Join cancer info
ctdna_df = pd.read_excel('../data/ctDNA pre surgery cohort 2.xlsx')
ctdna_data = ctdna_df.set_index("Patient ID").to_dict(orient="index")
locations = []
msi_statuses = []
metastasized = []
for code, subtype, other in zip(df['CODE'], df['Subtype'], df['other_info']):
    location, msi, has_metastasis = np.nan, np.nan, np.nan
    if subtype in SUBTYPE_MSI:
        location, msi = SUBTYPE_MSI[subtype]
    if code in ctdna_data:
        location = ctdna_data[code]['Location']
        msi = ctdna_data[code]['MSS/MSI status']
        has_metastasis = ('meta' in ctdna_data[code]['Reason for Surgery'].lower())
    if not isinstance(has_metastasis, bool):
        has_metastasis = False
    if other in {'pT3N0 (MSI-H)', 'pT2N0 (MSI-H)', 'pT4aN1 (MSI-H)', 'pT3N1b (MSI-H)'}:
        msi = 'MSI-High'
    locations.append(location)
    msi_statuses.append(msi)
    metastasized.append(has_metastasis)
df['Cancer-location'] = locations
df['MSS/MSI-status'] = msi_statuses
df['Metastasis'] = metastasized

# Standardize some columns
df['stage'] = df['stage'].apply(lambda s: STANDARDIZE_STAGE.get(s, s))
df['MSS/MSI-status'] = df['MSS/MSI-status'].apply(lambda s: STANDARDIZE_MSI_STATUS.get(s, s))
df['Cancer-location'] = df['Cancer-location'].apply(lambda s: STANDARDIZE_LOCATION.get(s, s))

# Format some columns numerically
df['Run_num'] = df['Run'].apply(lambda s: RUN_TO_DATETIME.get(s, np.nan))
df['sex_num'] = df['sex'].apply(lambda s: SEX_TO_FLOAT.get(s, np.nan))
df['AGE_num'] = df['AGE'].apply(lambda s: age_to_float(s))
df['stage_num'] = df['stage'].apply(lambda s: STAGE_TO_FLOAT.get(s, np.nan))

# Get colorectal cancer cases (COLCA)
mask = df['category'].apply(lambda s: s in {'cancer', 'Cancer'})
mask &= df['Cancer Type'].apply(lambda s: s in {'Colorectal cancer', 'colorectal cancer', 'Colorectal Cancer'})
mask &= df['Emseq_GCCODE'].apply(lambda s: s != 'NA')
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['stage'].apply(lambda s: s not in {'0', 0})
mask &= df['other info'].apply(lambda s: s not in {'Nivolumab - 12m', 'Nivolumab - 6w', 'Nivolumab - 2w', 'Nivolumab - 12w', 'Nivolumab - 4w', 'Nivolumab - 9m', 'Nivolumab - 6m', 'replicate', 'repl', 'Do not upload data, thienpont lab library', 'Do not upload data'})
mask &= df['library notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
df_cases = df[mask]

# Get controls
mask = df['category'].apply(lambda s: s in {'Healthy controls', 'GRP_healthy'})
mask &= df['Emseq_GCCODE'].apply(lambda s: s != 'NA')
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl', 'Do not upload data, thienpont lab library', 'Do not upload data'})
mask &= df['library notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
df_controls = df[mask]

# Match cases and controls
matched_cases, matched_controls = match(df_cases, df_controls, vars_cols=['AGE_num', 'sex_num', 'Run_num'], id_col='Emseq_GCCODE')
print(len(matched_cases), len(matched_controls))
print('Age difference', np.nanmean(matched_cases['AGE_num']) - np.nanmean(matched_controls['AGE_num']))
print('Sex difference', np.nanmean(matched_cases['sex_num']) - np.nanmean(matched_controls['sex_num']))
print('Run difference', np.nanmean(matched_cases['Run_num']) - np.nanmean(matched_controls['Run_num']))
print(set(matched_cases['other info']))
print(set(matched_controls['other info']))

# Reset indices to ensure positional correspondence
matched_cases = matched_cases.reset_index(drop=True)
matched_controls = matched_controls.reset_index(drop=True)
matched_df = pd.concat([matched_cases, matched_controls], ignore_index=True)

print(collections.Counter(matched_df['sex']))
print(collections.Counter(matched_cases['Metastasis']))
print(collections.Counter(matched_cases['Cancer-location']))
print(collections.Counter(matched_cases['MSS/MSI-status']))
exit()

# Define split based on matched control category
is_grp = (matched_controls["category"] == "GRP_healthy").to_numpy()

# Training set
train_df = pd.concat(
    [
        matched_cases.iloc[is_grp],
        matched_controls.iloc[is_grp]
    ],
    ignore_index=True
)

print(collections.Counter(train_df["category"]))
train_df.to_csv("training_set.csv", index=False)

# Test set
test_df = pd.concat(
    [
        matched_cases.iloc[~is_grp],
        matched_controls.iloc[~is_grp]
    ],
    ignore_index=True
)

print(collections.Counter(test_df["category"]))
test_df.to_csv("test_set.csv", index=False)