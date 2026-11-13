from typing import List, Tuple

import os
import numpy as np
import datetime
import collections
import pandas as pd
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedKFold
import lightgbm
from scipy.optimize import linear_sum_assignment


CONFIG_FOLDER = "../config"
os.makedirs(CONFIG_FOLDER, exist_ok=True)


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

SAMPLE_IDS_NO_PUC = {'GC140036', 'GC147297_GTTTGCTC-ACCAACAG', 'GC147297_AAGGCGTA-CCATCCGC', 'GC142065', 'GC146766_TACGGTCT-TCCATTGC', 'GC146769_GCCTTAAC-TAGGAGCT', 'GC137547', 'GC146664_GATCTTGC-CATGAGCA', 'GC146373', 'GC137540', 'GC147297_ATTACCCA-CACTGTAG', 'GC146387', 'GC146338', 'GC146769_GTCGATTG-AAGACCGT', 'GC146343', 'GC146328', 'GC142074', 'GC146365', 'GC146664_CACAGACT-CTGTACCA', 'GC146664_GTCATCGT-CCGCTTAA', 'GC146376', 'GC137541', 'GC147297_ATGGCTGT-CCTGTCAA', 'GC146384', 'GC146664_CGCAACTA-TCGGATTC', 'GC146665_CGAATTGC-GTAAGGTG', 'GC146665_CTCAGAAG-GGCAAGTT', 'GC146391', 'GC146345', 'GC146766_TTCCAGGT-CAGTGCTT', 'GC146351', 'GC146665_TTGCGAGA-TATGGCAC', 'GC146355', 'GC146665_GAACGAAG-GAATCACC', 'GC146665_CCAGTATC-AGCGAGAT', 'GC146329', 'GC145418', 'GC146347', 'GC161787_TCCTCATG-GTACACCT', 'GC146361', 'GC146664_GAACGGTT-GGAGGAAT', 'GC137544', 'GC145416', 'GC146665_GAGAGTAC-TGGAAGCA', 'GC146389', 'GC142078', 'GC146664_TGGTGAAG-CTTACAGC', 'GC146665_CTTAGGAC-CTTCACTG', 'GC146665_TCTACGCA-CAATAGCC', 'GC146333', 'GC146665_CATGAGCA-GATCTTGC', 'GC137553', 'GC146340', 'GC146665_GAATCCGT-CAACTCCA', 'GC146664_AGCTAAGC-GTCGTTAC', 'GC146664_TTACGTGC-AGACCTTG', 'GC146330', 'GC137546', 'GC146664_AGCTACCA-GTCCTAAG', 'GC146665_CAGGTTCA-ATAACGCC', 'GC161787_CCGCTTAA-GTCATCGT', 'GC161787_CTCAGAAG-GGCAAGTT', 'GC146390', 'GC140049', 'GC137543', 'GC142073', 'GC161786_TATGGCAC-TTGCGAGA', 'GC146665_TCATCTCC-CTTGCTAG', 'GC146344', 'GC146339', 'GC146364', 'GC146665_GCGTTAGA-CAGCATAC', 'GC146769_TCTAGGAG-AGTGACCT', 'GC137538', 'GC142062', 'GC146665_GTAACCGA-AAGACACC', 'GC146664_GGCAAGTT-CTCAGAAG', 'GC146381', 'GC147298_ATTACCCA-CACTGTAG', 'GC146769_AGCGAGAT-CCAGTATC', 'GC146337', 'GC146665_TACGGTCT-TCCATTGC', 'GC146356', 'GC146357', 'GC146766_TCGTCTGA-GTCCTTGA', 'GC146352', 'GC137550', 'GC146664_CTCGACTT-ATCTGACC', 'GC146769_TCCATTGC-TACGGTCT', 'GC146664_GATAGCCA-CAACACAG', 'GC145427', 'GC146769_ATAACGCC-CAGGTTCA', 'GC146353', 'GC146664_GGTGTACA-AGACGCTA', 'GC146354', 'GC145425', 'GC146386', 'GC146366', 'GC146378', 'GC137549', 'GC146766_TACTCCAG-GGTATAGG', 'GC146359', 'GC146380', 'GC146665_TAGGAGCT-GCCTTAAC', 'GC146382', 'GC137555', 'GC146362', 'GC137551', 'GC146377', 'GC137556', 'GC146766_AAGACCGT-GTCGATTG', 'GC146368', 'GC146665_GTCCTAAG-AGCTACCA', 'GC140053', 'GC146371', 'GC146664_CTCGAACA-GCAATTCC', 'GC161787_GTAAGGTG-CGAATTGC', 'GC146383', 'GC147298_AAGCGACT-ACGAATCC', 'GC146769_GAATCACC-GAACGAAG', 'GC146769_TATGGCAC-TTGCGAGA', 'GC146363', 'GC137548', 'GC146350', 'GC146385', 'GC146664_GGACATCA-TACCTGCA', 'GC146367', 'GC146769_GTCCTTGA-TCGTCTGA', 'GC146341', 'GC146665_GGAAGAGA-CGAGAGAA', 'GC146332', 'GC146664_CAATGCGA-GACGAACT', 'GC146664_CAACTCCA-GAATCCGT', 'GC146336', 'GC145423', 'GC146335', 'GC137539', 'GC146766_TAGGAGCT-GCCTTAAC', 'GC146664_GGTTGAAC-TCCACGTT', 'GC146370', 'GC146374', 'GC146665_GTCGTTAC-AGCTAAGC', 'GC140052', 'GC146375', 'GC145435', 'GC146769_GGTATAGG-TACTCCAG', 'GC146769_GTAAGGTG-CGAATTGC', 'GC146664_CAGCATAC-GCGTTAGA', 'GC142076', 'GC146665_GCAATTCC-CTCGAACA', 'GC145430', 'GC145429', 'GC146360', 'GC137557', 'GC146342', 'GC146334', 'GC146769_TGCGTAAC-AGCCTATC', 'GC146769_CGAGAGAA-GGAAGAGA', 'GC161787_CGAATTGC-GTAAGGTG', 'GC146664_CAATAGCC-TCTACGCA', 'GC146664_CGTATCTC-CAAGGTAC', 'GC147297_TCCCACGA-TTCAATAG', 'GC146369', 'GC146665_TACTCCAG-GGTATAGG', 'GC146372', 'GC142056', 'GC146769_CTTGCTAG-TCATCTCC', 'GC137552', 'GC147297_GTCTGAGT-GACTTGTG', 'GC137545', 'GC146331', 'GC146664_CCACAACA-GTACCACA', 'GC146766_CAGGTTCA-ATAACGCC', 'GC146348', 'GC146664_TGGAAGCA-GAGAGTAC', 'GC146349', 'GC146346', 'GC146665_AGACCTTG-TTACGTGC', 'GC161787_GTCCTAAG-AGCTACCA', 'GC146665_ATCTGACC-CTCGACTT', 'GC146769_CAGTGCTT-TTCCAGGT', 'GC146664_CTTCGGTT-AACCAGAG', 'GC161787_CAGCATAC-GCGTTAGA', 'GC147297_TGTTCGCC-AGGTAGGA', 'GC146664_AAGACACC-GTAACCGA', 'GC146379', 'GC146766_TTACCGAC-CGAATACG', 'GC140042', 'GC142079', 'GC147297_ACAAGGCA-TCGCGCAA', 'GC146769_CGAATACG-TTACCGAC', 'GC146358', 'GC137542', 'GC146664_CTTCACTG-CTTAGGAC', 'GC137554', 'GC146388'}

STANDARDIZE_CATEGORY = {
    "Healthy pregnant": "pregnant_control", "Cancer": "cancer"
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


def filename_from_sample_id(s: str) -> str:
    if s in {
        'GC121327', 'GC121330', 'GC118481', 'GC118484', 'GC121334',
        'GC121335', 'GC119384', 'GC114137', 'GC114131', 'GC118465', 'GC145430, GC161787_AACCAGAG-CTTCGGTT',
        'GC121344', 'GC118493', 'GC121329', 'GC121333', 'GC119367', 'GC121341', 'GC118496', 'GC114134',
        'GC119406', 'GC114127', 'GC114183', 'GC119408', 'GC145429, GC161786_TCCTCATG-GTACACCT', 'GC119380',
        'GC121346', 'GC114140', 'GC119414', 'GC119361', 'GC114133', 'GC118473', 'GC118466', 'GC118483',
        'GC119430', 'GC119394', 'GC119395', 'GC119390', 'GC119377', 'GC121339', 'GC145423, GC161786_TTGCGAGA-TATGGCAC',
        'GC119376', 'GC131894 ', 'GC121326', 'GC114135', 'GC114132', 'GC114142', 'GC118468', 'GC121348',
        'GC119405', 'GC119409', 'GC119344', 'GC119429', 'GC114139', 'GC119381', 'GC119389', 'GC114186',
        'GC118472', 'GC119369', 'GC121331', 'GC114129', 'GC121338', 'GC121342', 'GC119403', 'GC119382',
        'GC119423', 'GC119418', 'GC119420', 'GC145416, GC161786_TCATCTCC-CTTGCTAG', 'GC119410', 'GC118482',
        'GC145418, GC161786_TACTCCAG-GGTATAGG',
        'GC118486', 'GC118463', 'GC114138', 'GC121343', 'GC114185',
        'GC110688', 'GC112846', 'GC112847', 'GC112866', 'GC112867', 'GC112868', 'GC112869', 'GC115809',
        'GC110687', 'GC116718', 'GC116719', 'GC116720', 'GC116721', 'GC116722', 'GC116723', 'GC116724', 'GC116725',
        'GC116726', 'GC116727', 'GC116728', 'GC116729', 'GC116730', 'GC116731', 'GC116732', 'GC116733', 'GC116734'
    }:
        return s + ".bam"
    elif (len(s) == 9) and s.startswith("GC100"):
        return s + ".bam"
    return s + ".merged.rmdup.bam"


def reference_from_sample_id(s: str):
    if s in SAMPLE_IDS_NO_PUC:
        return 'ref1'
    else:
        return 'ref2'


def clean_dataframe(df):
    df["category"] = df["category"].apply(lambda s: STANDARDIZE_CATEGORY.get(s, s))
    df["sample_id"] = df["Emseq_GCCODE"].apply(lambda s: s.split(',')[0].strip())
    df["bam"] = df["sample_id"].apply(filename_from_sample_id)
    df["batch"] = df["Run"]
    df["reference"] = df["sample_id"].apply(reference_from_sample_id)
    df = df[["sample_id", "bam", "category", "cancer_type", "cancer_stage", "cancer_location", "batch", "reference"]]
    return df


# Load metadata
df = pd.read_excel('../data/Metadata_EMseqsamples.xlsx')

# Join cancer info
ctdna_df = pd.read_excel('../data/ctDNA pre surgery cohort 2.xlsx')
ctdna_data = ctdna_df.set_index("Patient ID").to_dict(orient="index")
locations = []
msi_statuses = []
metastasized = []
for code, subtype, other in zip(df['CODE'], df['Subtype'], df['other info']):
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
df['cancer_stage'] = df['stage'].apply(lambda s: STANDARDIZE_STAGE.get(s, s))
df['msi_status'] = df['MSS/MSI-status'].apply(lambda s: STANDARDIZE_MSI_STATUS.get(s, s))
df['cancer_location'] = df['Cancer-location'].apply(lambda s: STANDARDIZE_LOCATION.get(s, s))
df['cancer_type'] = df['Cancer Type']

# Format some columns numerically
df['Run_num'] = df['Run'].apply(lambda s: RUN_TO_DATETIME.get(s, np.nan))
df['sex_num'] = df['sex'].apply(lambda s: SEX_TO_FLOAT.get(s, np.nan))
df['AGE_num'] = df['AGE'].apply(lambda s: age_to_float(s))
df['stage_num'] = df['stage'].apply(lambda s: STAGE_TO_FLOAT.get(s, np.nan))

# Get pregnant controls
mask = df['category'].apply(lambda s: s == 'Healthy pregnant')
mask &= df['GA at sampling'].apply(lambda s: str(s) != 'nan')
mask &= df['Emseq_GCCODE'].apply(lambda s: s != 'NA')
mask &= df['Sample type'].apply(lambda s: s == 'cfDNA')
mask &= df['status'].apply(lambda s: s == 'sequenced')
mask &= df['other info'].apply(lambda s: s not in {'replicate', 'repl', 'Do not upload data, thienpont lab library', 'Do not upload data'})
mask &= df['library notes'].apply(lambda s: s not in {'Heitzer', 'NO PROCESS ARCHIVED CONSENSUS', 'NaOH, low library yield', 'BISULFITE', 'low library yield', 'Low library yield', 'NaOH, Low library yield'})
df_preg = df[mask]
df_preg_clean = clean_dataframe(df_preg)
df_preg_clean['GA'] = df[mask]['GA at sampling']
df_preg_clean.to_csv(os.path.join(CONFIG_FOLDER, "pregnant_controls.csv"), index=False)

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
train_df = clean_dataframe(train_df)
train_df.to_csv(os.path.join(CONFIG_FOLDER, "training_set.csv"), index=False)

# Test set
test_df = pd.concat(
    [
        matched_cases.iloc[~is_grp],
        matched_controls.iloc[~is_grp]
    ],
    ignore_index=True
)
print(collections.Counter(test_df["category"]))
test_df = clean_dataframe(test_df)
test_df.to_csv(os.path.join(CONFIG_FOLDER, "test_set.csv"), index=False)


df = pd.concat([df_preg_clean, train_df, test_df], axis=0, join="inner", ignore_index=True)
df.to_csv(os.path.join(CONFIG_FOLDER, "samples.csv"), index=False)