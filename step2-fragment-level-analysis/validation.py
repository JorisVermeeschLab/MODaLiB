import collections
from abc import ABCMeta, abstractmethod
from typing import Optional, List, Tuple
import tqdm

import matplotlib.pyplot as plt
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import pandas as pd
import seaborn
import tqdm
from sklearn.decomposition import KernelPCA
from sklearn.manifold import TSNE
from sklearn.metrics import matthews_corrcoef, roc_auc_score, roc_curve, confusion_matrix
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.model_selection import LeaveOneOut, PredefinedSplit
from sklearn.preprocessing import StandardScaler, LabelEncoder, KernelCenterer
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

BINARY = True

labels = ['Healthy', 'Breast', 'Cervical', 'Colorectal', 'Ovarian', 'Hematological', 'Solid (misc)', 'Skin']
N_CLASSES = len(labels)


class Dataset(metaclass=ABCMeta):

    @abstractmethod
    def to_kernel_matrix(self, kernel_idx: np.ndarray) -> np.ndarray:
        pass


class KernelDataset(Dataset):

    def __init__(self, filepaths: List[str]):
        self.filepaths: List[str] = filepaths
        self.n_kernel_types: int = 5  # TODO
        self.pca: bool = True

        self.tmp = {}

    def to_kernel_matrix(self, kernel_idx: np.ndarray) -> np.ndarray:

        if self.pca:
            X = []
            for filepath in self.filepaths:
                data = np.load(filepath, allow_pickle=True)
                Ks = data['K']
                n_bins = len(data['K']) // self.n_kernel_types
                for k in kernel_idx:
                    for j in range(n_bins):
                        key = (filepath, k, j)
                        if key not in self.tmp:
                            K = Ks[j * self.n_kernel_types + k, :, :]
                            try:
                                X.append(KernelPCA(n_components=50, kernel='precomputed').fit_transform(K))
                                X[-1] -= np.mean(X[-1], axis=0)[np.newaxis, :]
                                variance = np.std(X[-1][:, 0])
                                if variance > 0:
                                    X[-1] /= variance
                                self.tmp[key] = X[-1]

                            except ValueError:
                                pass
                        else:
                            X.append(self.tmp[key])
            X = np.concatenate(X, axis=1)
            # X = StandardScaler().fit_transform(X)
            return pairwise_kernels(X, metric='rbf')
        else:
            X = []
            for filepath in self.filepaths:
                data = np.load(filepath, allow_pickle=True)
                Ks = data['K']
                n_bins = len(data['K']) // self.n_kernel_types
                for k in kernel_idx:
                    for j in range(n_bins):
                        key = (filepath, k, j)
                        if key not in self.tmp:
                            K = Ks[j * self.n_kernel_types + k, :, :]
                            n = len(K)
                            #V = np.eye(n) - (1. / n)
                            #K = V @ K @ V
                            K = KernelCenterer().fit_transform(K)

                            # TODO: manage complex eigenvalues
                            variance = np.sum(np.abs(np.clip(np.linalg.eigvals(K), 0, None)))
                            if variance > 0:
                                K /= variance
                            X.append(K)
                            self.tmp[key] = X[-1]
                        else:
                            X.append(self.tmp[key])
            K = np.sum(X, axis=0)
            D = np.diagonal(K)[:, np.newaxis] + np.diagonal(K)[np.newaxis, :] - 2 * K
            np.fill_diagonal(D, 0)
            gamma = 1. / 2000
            return np.exp(-gamma * (D ** 2))
            return K


    def get_sample_names(self) -> np.ndarray:
        data = np.load(self.filepaths[-1], allow_pickle=True)
        return data['sample_names']


class ChromosomeBounds:

    @staticmethod
    def get_10kb():
        return np.asarray([
            0,        # chr1 (24896)
            24896,    # chr2 (24220)
            49116,    # chr3 (19830)
            68946,    # chr4 (19022)
            87968,    # chr5 (18154)
            106122,   # chr6 (17081)
            123203,   # chr7 (15935)
            139138,   # chr8 (14514)
            153652,   # chr9 (13840)
            167492,   # chr10 (13380)
            180872,   # chr11 (13509)
            194381,   # chr12 (13328)
            207709,   # chr13 (11437)
            219146,   # chr14 (10705)
            229851,   # chr15 (10200)
            240052,   # chr16 (9034)
            249085,   # chr17 (8326)
            257411,   # chr18 (8038)
            265449,   # chr19 (5862)
            271311,   # chr20 (6445)
            277756,   # chr21 (4671)
            282427,   # chr22 (5082)
            287509])  # end


# Load GC-content of reference genome
gc_content = np.load('gc-content.npy')
mappability = np.load('mappability.npy')

df = pd.read_csv('../samples-details.tsv', delimiter='\t')
categories = df['category'].to_numpy()
gc_codes = df['GC CODE'].to_numpy()
cat = {}
metadata = {}
for gc_code, category in zip(gc_codes, categories):
    if isinstance(gc_code, str) and gc_code.startswith('GC'):
        row = df[df['GC CODE'] == gc_code]
        cat[gc_code] = category
        metadata[gc_code] = row.iloc[0]

# Load feature data
data = np.load('X.npz', allow_pickle=True)
sample_names = data['sample_names']
X = data['X'].item()

"""
gc_content = np.load('mappability-old.npy')
print(gc_content.shape)
new_gc_content = []
bounds = ChromosomeBounds.get_10kb()
for i in range(len(bounds) - 1):
    start = bounds[i]
    end = bounds[i + 1]
    gc = gc_content[start:end]
    new_gc = []
    tmp = []
    for j in range(len(gc)):
        tmp.append(gc[j])
        if len(tmp) == 5:
            new_gc.append(np.mean(tmp))
            tmp = []
    if len(tmp) > 0:
        new_gc.append(np.mean(tmp))
    target_length = X[i + 1][0].shape[1]
    while len(new_gc) < target_length:
        new_gc.append(new_gc[-1])
    new_gc_content += new_gc[:target_length]
new_gc_content = np.asarray(new_gc_content)
"""

# Load kernel data
dataset2 = KernelDataset([f'../kernels/K-{k}.npz' for k in range(1, 26)])
sample_names2 = dataset2.get_sample_names()

# Make sure the samples in K and X actually match
sample_dict = {sample_name: i for i, sample_name in enumerate(sample_names2)}
sample_dict['GC1002975'] = sample_dict['GC110694']  # TODO
sample_dict['GC1002979'] = sample_dict['GC110695']  # TODO
X_new = {c: [[] for s in range(5)] for c in range(1, 23)}
for sample_name in sample_names2:
    for c in range(1, 23):
        for s in range(5):
            if sample_name in sample_dict:
                X_new[c][s].append(X[c][s][sample_dict[sample_name]])
            else:
                import random  # TODO
                X_new[c][s].append(X[c][s][random.choice(list(sample_dict.values()))])
for c in range(1, 23):
    for s in range(5):
        X_new[c][s] = np.asarray(X_new[c][s])
X = X_new

sequencing_run = []
sequencing_depth = []
age = []
sex = []
y = []
is_grp = []
cfdna_mask = []
categories = []
for sample_name in sample_names:
    row = metadata[sample_name]

    # Age
    if isinstance(row['AGE'], str):
        if '-' in row['AGE']:
            age.append(0.5 * float(row['AGE'].split('-')[0]) + 0.5 * float(row['AGE'].split('-')[1]))
        elif np.isnan(float(row['AGE'])):
            age.append(0)
        else:
            age.append(float(row['AGE']))
    elif np.isnan(float(row['AGE'])):
        age.append(0)
    else:
        age.append(float(row['AGE']))

    # Sex
    if row['sex'] == 'M':
        sex.append(0)
    elif row['sex'] in {'V', 'F'}:
        sex.append(1)
    else:
        sex.append(2)

    # Run
    if isinstance(row['Run'], float) and np.isnan(row['Run']):
        sequencing_run.append(-1)
    else:
        sequencing_run.append(row['Run'])

    # Sequencing depth
    if row['Depth required'] == '10X':
        sequencing_depth.append(0)
    elif row['Depth required'] == '30X':
        sequencing_depth.append(1)
    else:
        sequencing_depth.append(2)

    if sample_name not in cat:
        categories.append('healthy')
        is_grp.append(False)
    else:
        category = cat[sample_name].lower()
        is_grp.append('grp' in category)
        if 'healthy' in category:
            categories.append('healthy')
        else:
            categories.append('cancer')
    cfdna_mask.append(row['Sample type'] == 'cfDNA')
    if isinstance(row['Cancer Type'], float) or (row['Cancer Type'] == 'healthy'):
        y.append(0)
        assert 'healthy' in cat[sample_name].lower()
    elif 'breast' in row['Cancer Type'].lower():
        y.append(1)
    elif 'cervical' in row['Cancer Type'].lower():
        y.append(2)
    elif 'colorectal' in row['Cancer Type'].lower():
        y.append(3)
    elif 'ovarian' in row['Cancer Type'].lower():
        y.append(4)
    elif ('dlbcl' in row['Cancer Type'].lower()) or ('lymfeklier' in row['Cancer Type'].lower()) or ('haematological' in row['Cancer Type'].lower()) or ('myelodysplasitc' in row['Cancer Type'].lower()) or ('myeloma' in row['Cancer Type'].lower()) or ('lymphocytosis' in row['Cancer Type'].lower()) or ('hodgkin' in row['Cancer Type'].lower()):
        y.append(5)
    elif 'prostate' in row['Cancer Type'].lower():
        y.append(6)
    elif 'lip cancer' in row['Cancer Type'].lower():
        y.append(6)
    elif ('carcinoma' in row['Cancer Type'].lower()) or ('carcinoom' in row['Cancer Type'].lower()) or ('skin' in row['Cancer Type'].lower()) or ('melanoma' in row['Cancer Type'].lower()):
        y.append(7)
    elif 'endometrium' in row['Cancer Type'].lower():
        y.append(6)
    elif 'bladder' in row['Cancer Type'].lower():
        y.append(6)
    elif 'brain' in row['Cancer Type'].lower():
        y.append(6)
    elif 'pancreas' in row['Cancer Type'].lower():
        y.append(6)
    elif 'tongue' in row['Cancer Type'].lower():
        y.append(6)
    elif 'thyroid' in row['Cancer Type'].lower():
        y.append(6)
    elif 'kidney' in row['Cancer Type'].lower():
        y.append(6)
    else:
        raise NotImplementedError(row['Cancer Type'])

sequencing_depth = np.asarray(sequencing_depth, dtype=float)
sequencing_run = LabelEncoder().fit_transform(sequencing_run)
sex = np.asarray(sex, dtype=int)
age = np.asarray(age, dtype=float)
confounders = np.asarray([
    sequencing_depth,
    # sequencing_run,
    sex, age
]).T
categories = np.asarray(categories, dtype=object)
cfdna_mask = np.asarray(cfdna_mask, dtype=bool)
is_grp = np.asarray(is_grp, dtype=bool)
y = np.asarray(y, dtype=int)
print('Labels', collections.Counter(y))
if BINARY:
    y = (y > 0).astype(int)


def lowess_correction(X: np.ndarray, exog: np.ndarray, desc: str, frac: float) -> np.ndarray:
    exog = np.round(exog * 1000).astype(int) // 10
    X_adapted = np.copy(X)
    for j in range(X.shape[1]):
        for i in tqdm.tqdm(range(len(X)), desc=desc):
            endog = X[i, j, :]

            y_pred = lowess(endog, exog, frac=frac, return_sorted=False)
            mask = (y_pred > 0)
            X_adapted[i, j, mask] = X[i, j, mask] / y_pred[mask]
            mask = (exog == 0)
            X_adapted[i, j, mask] = X[i, j, mask]

    return X_adapted


def gc_correction(X: np.ndarray, mappability: np.ndarray) -> np.ndarray:
    return lowess_correction(X, mappability, 'GC correction', 0.3)


def evaluate(kernel_idx):

    X_ = np.concatenate([[X[chr_id][i] for i in kernel_idx] for chr_id in range(1, 23)], axis=2)  # TODO

    X_ = np.transpose(X_, (1, 0, 2))

    print(X_.shape)

    # X_ = gc_correction(X_, gc_content)
    # X_ = gc_correction(X_, mappability)

    X_ = X_.reshape(len(X_), -1)


    X_ = StandardScaler().fit_transform(X_)
    K_f = pairwise_kernels(X_, metric='rbf')


    K_agg = dataset2.to_kernel_matrix(kernel_idx)


    #alpha = 0.2
    #np.fill_diagonal(K_agg, alpha * np.diagonal(K_agg) + (1 - alpha) * np.ones(len(K_agg)))

    # Cross-validation
    y_target, y_pred1, y_pred2 = [], [], []

    cv = PredefinedSplit(sequencing_run)
    # cv = LeaveOneOut()
    #cv = KFold(n_splits=5)
    for i, (train_index, test_index) in tqdm.tqdm(enumerate(cv.split(X_))):
        X_train = X_[train_index, :]
        X_test = X_[test_index, :]
        K_train = K_agg[train_index, :][:, train_index]
        K_test = K_agg[test_index, :][:, train_index]
        Kf_train = K_f[train_index, :][:, train_index]
        y_train = y[train_index]
        Kf_test = K_f[test_index, :][:, train_index]
        y_test = y[test_index]

        model = SVC(kernel='precomputed', probability=True, C=1000, class_weight='balanced')

        model.fit(Kf_train, y_train)
        y_pred1.append(model.predict_proba(Kf_test))
        y_pred1[-1] = np.concatenate((y_pred1[-1], np.zeros((y_pred1[-1].shape[0], N_CLASSES - y_pred1[-1].shape[1]))), axis=1)

        model.fit(K_train, y_train)  # TODO
        y_pred2.append(model.predict_proba(K_test))
        y_pred2[-1] = np.concatenate((y_pred2[-1], np.zeros((y_pred2[-1].shape[0], N_CLASSES - y_pred2[-1].shape[1]))), axis=1)
        y_target.append(y_test)
    y_target = np.squeeze(np.concatenate(y_target, axis=0))
    y_pred1 = np.squeeze(np.concatenate(y_pred1, axis=0))
    y_pred2 = np.squeeze(np.concatenate(y_pred2, axis=0))
    print(f'MCC 1: {matthews_corrcoef(y_target, np.argmax(y_pred1, axis=1))}')
    print(f'MCC 2: {matthews_corrcoef(y_target, np.argmax(y_pred2, axis=1))}')
    if BINARY:
        auc = roc_auc_score(y_target, y_pred1[:, 1])
        roc = roc_curve(y_target, y_pred1[:, 1])
        roc2 = roc_curve(y_target, y_pred2[:, 1])
        print(f'AUC 1: {auc}')
        print(f'AUC 2: {roc_auc_score(y_target, y_pred2[:, 1])}')
    else:
        auc = roc_auc_score(y_target, y_pred1, multi_class='ovr')
        roc = roc_curve(y_target > 0, 1 - y_pred1[:, 0])
        roc2 = roc_curve(y_target > 0, 1 - y_pred2[:, 0])
        print(f'AUC 1: {auc}')
        print(f'AUC 2: {roc_auc_score(y_target, y_pred2, multi_class="ovr")}')

    if not BINARY:
        cm = confusion_matrix(y_target, np.argmax(y_pred2, axis=1))
        plt.figure(figsize=(6, 6))
        ax = plt.subplot()
        seaborn.heatmap(cm, annot=True, fmt='g', ax=ax, cbar=False)
        ax.set_xlabel('Predictions')
        ax.set_ylabel('Annotations')
        ax.xaxis.set_ticklabels(labels, rotation=90)
        ax.yaxis.set_ticklabels(labels, rotation=0)
        plt.tight_layout()
        plt.savefig('cm.png', dpi=200)
        plt.clf()
        plt.cla()
        plt.close()

    return auc, roc, roc2


settings = [
    ([0], 'Nucleosomes', 'purple', '-'),
    #([1, 2, 3, 4], None, 'purple', '--'),
    ([1], 'Fragment lengths', 'royalblue', '-'),
    #([0, 2, 3, 4], None, 'royalblue', '--'),
    ([2], 'Methylation', 'gold', '-'),
    #([0, 1, 3, 4], None, 'gold', '--'),
    ([3], 'Read counts', 'forestgreen', '-'),
    #([0, 1, 2, 4], None, 'forestgreen', '--'),
    ([4], 'Mutations', 'crimson', '-'),
    #([0, 1, 2, 3], None, 'crimson', '--'),
    ([0, 1, 2, 3, 4], 'All', 'black', '-'),
]


fig, ax = plt.subplots(1, 2, figsize=(16, 8))
for kernel_idx, label, color, linestyle in settings:
    auc, roc, roc2 = evaluate(kernel_idx)
    if linestyle == '-':
        ax[0].plot(roc[0], roc[1], alpha=0.6, label=label, color=color, linestyle=linestyle)
        ax[1].plot(roc2[0], roc2[1], alpha=0.6, label=label, color=color, linestyle=linestyle)
# plt.tight_layout()
plt.legend()
ax[0].set_xlabel('False positive rate')
ax[0].set_ylabel('True positive rate')
ax[1].set_xlabel('False positive rate')
ax[1].set_ylabel('True positive rate')
ax[0].set_title('RBF kernel')
ax[1].set_title('Custom fragment-level kernel')
ax[0].grid(alpha=0.4, linestyle='--', color='grey', linewidth=0.5)
ax[1].grid(alpha=0.4, linestyle='--', color='grey', linewidth=0.5)
ax[0].spines[['right', 'top']].set_visible(False)
ax[1].spines[['right', 'top']].set_visible(False)
plt.savefig('roc.png', dpi=200)
plt.show()


kernel_idx = [0, 1, 2, 3, 4]  # TODO
X_ = np.concatenate([[X[chr_id][i] for i in kernel_idx] for chr_id in range(1, 23)], axis=2)
X_ = np.transpose(X_, (1, 0, 2))
"""
print(X_.shape)
p_values = np.empty((4, X_.shape[2]))
for k in range(4):
    for j in tqdm.tqdm(range(X_.shape[2])):
        mask = (y > 0)
        p_values[k, j] = ks_2samp(X_[mask, k, j], X_[~mask, k, j]).pvalue
np.save('pvalues-10x.npy', p_values)
"""
X_ = X_.reshape(len(X_), -1)
X_ = StandardScaler().fit_transform(X_)
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
coords = KernelPCA(kernel='rbf').fit_transform(X_)
#coords = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=2).fit_transform(X)
mask = (categories == 'healthy')
colors = ['royalblue', 'orange', 'crimson', 'violet', 'pink', 'purple', 'green', 'brown']
labels = ['Healthy', 'Breast', 'Cervical', 'Colorectal', 'Ovarian', 'Hematological', 'Solid (misc)', 'Skin']
for i in range(len(labels)):
    mask = (y == i)
    artist = ax.scatter(coords[mask, 0], coords[mask, 1], color=colors[i], picker=5, alpha=0.4, label=labels[i])
    # artist.obj = sample_names[i]
ax.grid(linestyle='--', color='grey', alpha=0.4, linewidth=0.5)
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('First kernel PC')
ax.set_ylabel('Second kernel PC')
ax.set_title('RBF kernel')
ax.legend(prop={'size': 6})
plt.tight_layout()
plt.savefig('kpca-0.1x.png', dpi=200)



# model = AffinityPropagation()
"""
model = KMeans(n_clusters=2)
y_pred = model.fit_predict(X)
print(f'Adjusted Rand score 1: {adjusted_rand_score(y, y_pred)}')
print(f'Davies-Bouldin score 1: {davies_bouldin_score(X, y)}')
print(f'Calinski-Harabasz score 1: {calinski_harabasz_score(X, y)}')
y_pred = model.fit_predict(KernelPCA(kernel='precomputed').fit_transform(K_agg))
print(f'Adjusted Rand score 2: {adjusted_rand_score(y, y_pred)}')
print(f'Davies-Bouldin score 2: {davies_bouldin_score(KernelPCA(kernel="precomputed").fit_transform(K_agg), y)}')
print(f'Calinski-Harabasz score 2: {calinski_harabasz_score(KernelPCA(kernel="precomputed").fit_transform(K_agg), y)}')
"""

def on_pick(event):
    gc_code = event.artist.obj
    row = metadata[gc_code]
    print(row['GC CODE'], row['category'], row['Cancer Type'], row['Subtype'], row['stage'], row['AGE'], row['sex'], row['Sample type'])

fig, ax = plt.subplots(1, 2, figsize=(16, 8))
# coords = KernelPCA(kernel='rbf').fit_transform(X_)
coords = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=2).fit_transform(X_)
mask = (categories == 'healthy')
for i in range(len(mask)):
    color = 'royalblue' if mask[i] else 'orange'
    label = 'Healthy' if mask[i] else 'Cancer'
    artist = ax[0].scatter(coords[i, 0], coords[i, 1], color=color, picker=5, alpha=0.4, label=label)
    artist.obj = sample_names[i]
ax[0].grid(linestyle='--', color='grey', alpha=0.4)
ax[0].set_xlabel('First kernel PC')
ax[0].set_ylabel('Second kernel PC')
ax[0].set_title('RBF kernel')
ax[0].spines[['right', 'top']].set_visible(False)
# ax[0].legend()
w = np.array([24896, 24220, 19830, 19022, 18154, 17081, 15935, 14514, 13840,
              13380, 13509, 13328, 11437, 10705, 10201, 9033, 8326, 8038,
              5862, 6445, 4671, 5082])
w = w / np.sum(w)
K_agg = np.sum(w[:, np.newaxis, np.newaxis, np.newaxis] * K, axis=0)
K_agg = K_agg[kernel_idx, ...]
K_agg = np.mean(K_agg, axis=0)

print(np.mean(np.isnan(K_agg)))
plt.imshow(K_agg)
plt.show()



alpha = 1
#np.fill_diagonal(K_agg, (1 - alpha) * np.diagonal(K_agg) + alpha)
#coords = KernelPCA(kernel='precomputed').fit_transform(K_agg)
D = np.sqrt(np.diagonal(K_agg)[np.newaxis, :] + np.diagonal(K_agg)[:, np.newaxis] - 2 * K_agg)
coords = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=2, metric='precomputed').fit_transform(D)
mask = (y == 0)

in_legend = set()
for i in range(len(mask)):
    color = 'royalblue' if mask[i] else 'orange'
    label = 'Healthy' if mask[i] else 'Cancer'
    kwargs = dict(color=color, picker=5, alpha=0.4)
    if color not in in_legend:
        in_legend.add(color)
        kwargs['label'] = label
    artist = ax[1].scatter(coords[i, 0], coords[i, 1], **kwargs)
    artist.obj = sample_names[i]
ax[1].grid(linestyle='--', color='grey', alpha=0.4)
ax[1].set_xlabel('First kernel PC')
ax[1].set_ylabel('Second kernel PC')
ax[1].set_title('Custom fragment-level kernel')
ax[1].spines[['right', 'top']].set_visible(False)
ax[1].legend()
# fig.canvas.callbacks.connect('pick_event', on_pick)
plt.savefig('kpca.png', dpi=200)
plt.show()
