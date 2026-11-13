from typing import Tuple, List
import os
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.preprocessing import LabelEncoder

from modalibpy.deconvolution.metdecode import MetDecode, load_input_file
from modalibpy.stages import Annotation, BreastCancerStage
from modalibpy.dataset import Dataset
from modalibpy.gene_expression.data import GeneExpressionData


def load_legacy_dataset(filepaths: List[str], root: str = ".") -> Dataset:
    cohort = pd.read_excel(os.path.join(root, "paper", "metadata", "ctDNA pre surgery cohort 2.xlsx"))
    extra_age = dict(zip(cohort["Patient ID"].astype(str), cohort["Age"]))
    extra_sex = {str(code): sex[0] for code, sex in zip(cohort["Patient ID"], cohort["Gender"])}
    extra_msi = dict(zip(cohort["Patient ID"].astype(str), cohort["MSS/MSI status"]))
    extra_location = dict(zip(cohort["Patient ID"].astype(str), cohort["Location"]))
    extra_reason = dict(zip(cohort["Patient ID"].astype(str), cohort["Reason for Surgery"]))
    metadata = pd.read_excel(os.path.join(root, "paper", "metadata", "Metadata_EMseqsamples.xlsx"))

    ichor, annotations, ages, sexes, batches, notes = {}, {}, {}, {}, {}, {}
    for _, row in metadata.iterrows():
        gc_code = row["GC CODE"]
        patient_id = row["CODE"]
        patient_id = patient_id.strip() if isinstance(patient_id, str) else str(patient_id) if isinstance(patient_id, int) else ""

        if row["Cancer Type"] in {"Breast cancer", "breast cancer"}:
            subtype = row["Subtype"] if isinstance(row["Subtype"], str) else ""
            annotation = BreastCancerStage(True, row["stage"], subtype)
        elif isinstance(row["Cancer Type"], str) and row["Cancer Type"].lower() == "colorectal cancer":
            location = extra_location.get(patient_id, "")
            location = location if isinstance(location, str) else ""
            cancer_stage = row["stage"]
            reason = extra_reason.get(patient_id, "")
            if not location and isinstance(reason, str) and "liver" in reason.lower().replace(" ", ""):
                location, cancer_stage = "Liver metastasis", "IV"
            annotation = Annotation(True, cancer_stage, location, cancer_type="colorectal")
            subtype = row["Subtype"]
            if isinstance(subtype, str):
                if "MSS" in subtype:
                    annotation.msi = "MSS"
                elif "MSI" in subtype:
                    annotation.msi = "MSI-H"
            if patient_id in extra_msi:
                annotation.msi = extra_msi[patient_id]
                if not isinstance(annotation.msi, str):
                    annotation.msi = ""
        elif isinstance(row["category"], str) and row["category"] in {
            "cancer", "Cancer", "GRP_cancer", "grp_cancer"
        }:
            annotation = Annotation(True, row["stage"], cancer_type=row["Cancer Type"].lower())
        else:
            annotation = Annotation()
        annotations[gc_code] = annotation

        age = row["AGE"]
        if isinstance(age, str):
            if "-" in age:
                age = np.mean([float(value) for value in age.split("-")])
            elif age == "?":
                age = np.nan
            else:
                age = float(age)
        elif np.isnan(age):
            age = np.nan
        else:
            age = float(age)
        if np.isnan(age) and str(row["CODE"]) in extra_age:
            age = extra_age[str(row["CODE"])]
        ages[gc_code] = age

        sex = row["sex"]
        if isinstance(sex, str):
            if sex.upper() == "F":
                sex = 1
            elif sex.upper() == "M":
                sex = 0
            else:
                raise NotImplementedError(f'Unknown sex "{sex}"')
        else:
            sex = np.nan
        if np.isnan(sex) and str(row["CODE"]) in extra_sex:
            sex = extra_sex[str(row["CODE"])] == "F"
        sexes[gc_code] = sex

        batches[gc_code] = row["Run"]
        notes[gc_code] = row["notes"]
        ichor[gc_code] = row["Ichore - methbased"]

    gc_codes = [os.path.splitext(os.path.basename(filepath))[0] for filepath in filepaths]

    def values_for_sample(values: dict, gc_code: str):
        return values[gc_code] if gc_code in values else values[gc_code.split("_")[0]]

    dataset_sexes = np.asarray([values_for_sample(sexes, code) for code in gc_codes], dtype=float)
    dataset_ages = np.asarray([values_for_sample(ages, code) for code in gc_codes], dtype=float)
    dataset_ichor = np.asarray([values_for_sample(ichor, code) for code in gc_codes], dtype=float)
    dataset_annotations = np.asarray([values_for_sample(annotations, code) for code in gc_codes], dtype=object)
    dataset_batches = np.asarray([values_for_sample(batches, code) for code in gc_codes], dtype=str)
    dataset_notes = np.asarray([values_for_sample(notes, code) for code in gc_codes], dtype=object)
    groups = np.squeeze(LabelEncoder().fit_transform(dataset_batches))

    dataset = Dataset(filepaths, dataset_annotations, groups=groups)
    dataset.sexes = dataset_sexes
    dataset.ages = dataset_ages
    dataset.ichor = dataset_ichor
    dataset.notes = dataset_notes
    dataset.gc_codes = np.asarray(gc_codes, dtype=object)
    return dataset


def _load_dataset_from_metadata(metadata_df: pd.DataFrame, region_name: str, results_folder: str) -> Dataset:

    filepaths = []
    annotations = []
    batches = []
    for sample_id, cancer_type, cancer_stage in zip(metadata_df['sample_id'], metadata_df['cancer_type'], metadata_df['cancer_stage']):
        filepaths.append(os.path.join(results_folder, "features", sample_id, f"{region_name}.npz"))
        if cancer_type == 'Healthy':
            annotation = Annotation(is_cancer=False)
        else:
            annotation = Annotation(is_cancer=True, grade=cancer_stage, cancer_type=cancer_type)
        annotations.append(annotation)
    
    groups = LabelEncoder().fit_transform(metadata_df['batch'].to_numpy())

    dataset = Dataset(filepaths, annotations, groups=groups)
    #dataset.sexes = sexes
    #dataset.ages = ages
    #dataset.ichor = ichor
    #dataset.batches = batches
    #dataset.notes = notes
    dataset.gc_codes = np.asarray(list(metadata_df['sample_id'].to_numpy()), dtype=object)
    
    return dataset


def load_expression_reference(results_folder: str) -> Tuple[np.ndarray, List[str], List[str], List[str]]:
    RE = np.load(os.path.join(results_folder, 'reference_expression/mat.npz'), allow_pickle=True)
    mask = np.asarray([cell_type not in {'Unknown'} for cell_type in RE['cell_types']], dtype=bool)
    cell_types = RE['cell_types'][mask]
    compartments = RE['compartments'][mask]
    tissues = RE['tissues'][mask]
    mu = RE['mu'][mask]
    RE = GeneExpressionData(mu, RE['gene_symbols'])
    return RE, cell_types, compartments, tissues


def load_metdecode_reference(regions_folder: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dmrs = []
    with open(os.path.join(regions_folder, 'metdecode/dmrs.csv'), 'r') as f:
        for line in f:
            elements = line.rstrip().split(',')
            if len(elements) < 3:
                continue
            dmrs.append(f'chr{elements[0]}-{elements[1]}-{elements[2]}')

    M_atlas, D_atlas, cell_types, marker_names = load_input_file(os.path.join(regions_folder, 'metdecode/atlas.tsv'))
    metdecode_cell_types = np.asarray(cell_types, dtype=object)

    mask = np.asarray([1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1], dtype=bool)
    metdecode_cell_types = metdecode_cell_types[mask]
    M_atlas = M_atlas[mask]
    D_atlas = D_atlas[mask]

    marker_dict = {marker_name: i for i, marker_name in enumerate(marker_names)}
    M_atlas, D_atlas = MetDecode.add_pseudo_counts(M_atlas, D_atlas)
    R_atlas = M_atlas / D_atlas

    M_sorted, D_sorted, R_sorted = [], [], []
    for dmr in dmrs:
        k = marker_dict[dmr]
        R_sorted.append(R_atlas[:, k])
        M_sorted.append(M_atlas[:, k])
        D_sorted.append(D_atlas[:, k])
    R_atlas = np.asarray(R_sorted).T
    M_atlas = np.asarray(M_sorted).T
    D_atlas = np.asarray(D_sorted).T
    return R_atlas, M_atlas, D_atlas


def load_dataset(metadata_df: pd.DataFrame, results_folder: str, regions_folder: str) -> Dataset:

    # Load data
    dataset = _load_dataset_from_metadata(metadata_df, "50kb-noXY", results_folder)
    lib_size = np.sum(dataset.data['Coverage'], axis=1)
    dataset.data['Library size'] = lib_size

    # Methylation-based cell type deconvolution
    R_atlas, _, _ = load_metdecode_reference(regions_folder)
    dataset_meth = _load_dataset_from_metadata(metadata_df, "meth", results_folder)
    for key in ['CpG Methylation (first moment)', 'CpG Methylation (second moment)', 'CH Methylation (first moment)', 'CpG count']:
        dataset.data[key] = dataset_meth.data[key]
    dataset_deconvolution = _load_dataset_from_metadata(metadata_df, "metdecode-noXY", results_folder)
    dataset.data['CpG Methylation (deconvolution, cfDNA)'] = dataset_deconvolution.data['CpG Methylation (first moment)']
    dataset.data['CpG count (deconvolution, cfDNA)'] = dataset_deconvolution.data['CpG count']
    dataset.data['CpG Methylation (deconvolution, reference)'] = R_atlas
    dataset.data['GC-content (deconvolution, cfDNA)'] = dataset_deconvolution.data['GC-content']

    # Fragmentation-based cell type deconvolution
    RE, cell_types, compartments, tissues = load_expression_reference(results_folder)
    filepath = os.path.join(regions_folder, 'gene-bodies-no-XY.csv')
    df = pd.read_csv(filepath, sep=',', names=['chr_id', 'start', 'end', 'gene'])
    gene_names = np.asarray(list(df['gene'].to_numpy()), dtype=object)
    dataset_genes = _load_dataset_from_metadata(metadata_df, "genebodies-noXY", results_folder)
    X = dataset_genes.data['Proportion of long fragments']
    n = dataset_genes.data['Coverage']
    idx = np.argsort(np.min(n, axis=0))[-5000:]
    XE = GeneExpressionData(X[:, idx], gene_names[idx])
    XE_, RE_ = XE.match(RE)
    dataset.data['Genes'] = XE_.gene_names
    dataset.data['Fragmentation-based-gene-expression'] = XE_.X
    dataset.data['Gene-expression-reference'] = RE_.X
    dataset.data['Gene-expression-reference-cell-types'] = cell_types
    dataset.data['Gene-expression-reference-compartments'] = compartments
    dataset.data['Gene-expression-reference-tissues'] = tissues

    # Side information
    dataset.y = np.asarray([annot.is_cancer for annot in dataset.annotations], dtype=int)
    
    return dataset
