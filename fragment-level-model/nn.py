from typing import Optional, List, Tuple

import os
import argparse
import sys
import tqdm
import random
import time
import json

import numpy as np
import scipy.special
import torch
import captum.attr

from modalib.core import Sample, ReferenceNucMap, chr_id_to_chr_name


ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, '..', 'step1-bam-processing', 'out_')
OUT_DIR = os.path.join(ROOT, '..', 'outputs')
REF_GENOME_PATH = os.path.join(ROOT, '..', 'references', 'hg38.fa')
nuc_filepath = os.path.join(ROOT, '..', 'nucleosome-peaks', 'GSE71378_CA01.peaks.hg38.txt')

TF_DIR = os.path.join(ROOT, 'tfs')
os.makedirs(TF_DIR, exist_ok=True)

GC_BREAST_CANCERS = ['GC1002311', 'GC116964', 'GC116969', 'GC116967', 'GC116965', 'GC116976', 'GC116968', 'GC116960', 'GC116971', 'GC116966', 'GC116973', 'GC116970', 'GC110693', 'GC116972', 'GC116975', 'GC116974', 'GC116963', 'GC116962', 'GC110680', 'GC110696', 'GC110681']
GC_FEMALE_CONTROLS = ['GC119424', 'GC119419', 'GC119414', 'GC119413', 'GC119353', 'GC119412', 'GC119418', 'GC119422', 'GC119428', 'GC119426', 'GC119416', 'GC119355', 'GC119356', 'GC119423', 'GC119430', 'GC119420', 'GC119354', 'GC114128', 'GC119427', 'GC114129', 'GC119429', 'GC118494', 'GC116991', 'GC119351', 'GC118492', 'GC119388', 'GC114130', 'GC114133', 'GC118498', 'GC118464', 'GC119357', 'GC118485', 'GC118467', 'GC119343', 'GC118497', 'GC114132', 'GC118468', 'GC119344', 'GC114127', 'GC118466', 'GC119340', 'GC119359', 'GC119342', 'GC114184', 'GC118463', 'GC119373', 'GC118469', 'GC119406', 'GC119379', 'GC119392', 'GC119358', 'GC119368', 'GC119402', 'GC119403', 'GC119408', 'GC119369', 'GC119376', 'GC119386', 'GC119374', 'GC119395', 'GC119400', 'GC114131', 'GC119411', 'GC119391', 'GC119401', 'GC119370', 'GC119399', 'GC119394', 'GC119387', 'GC119380', 'GC119404', 'GC119378', 'GC119384', 'GC119405', 'GC118479', 'GC119377', 'GC118488', 'GC118475', 'GC119381', 'GC118478', 'GC118481', 'GC118483', 'GC114183', 'GC114115', 'GC114134', 'GC118477']
GC_BREAST_DATASET = GC_FEMALE_CONTROLS + GC_BREAST_CANCERS

GC_COLORECTAL_DATASET_TRAIN = ['GC121326', 'GC121327', 'GC121328', 'GC121329', 'GC121330', 'GC137538', 'GC137539', 'GC137540', 'GC137541', 'GC137542', 'GC137543', 'GC137544', 'GC137545', 'GC137546', 'GC137547', 'GC137548', 'GC137549', 'GC137550', 'GC137551', 'GC137552', 'GC137553', 'GC137554', 'GC137555', 'GC137556', 'GC146328', 'GC146329', 'GC146330', 'GC146331', 'GC146332', 'GC146333', 'GC146334', 'GC146335', 'GC146336', 'GC146337', 'GC146338', 'GC146339', 'GC146340', 'GC146341', 'GC146342', 'GC146343', 'GC146344', 'GC146345', 'GC146346', 'GC146347', 'GC146348', 'GC146349', 'GC146350', 'GC146351', 'GC146352', 'GC146353', 'GC146354', 'GC146355', 'GC146356', 'GC146357', 'GC146358', 'GC146359', 'GC146360', 'GC146361', 'GC146362', 'GC146363', 'GC146364', 'GC146365', 'GC146366', 'GC146367', 'GC146368', 'GC146369', 'GC146370', 'GC146371', 'GC146372', 'GC146373', 'GC146374', 'GC146375', 'GC146376', 'GC146377', 'GC146378', 'GC146379', 'GC146380', 'GC146383', 'GC146384', 'GC146385', 'GC146386', 'GC146387', 'GC146388', 'GC146389', 'GC146390', 'GC146391', 'GC146664_AAGACACC-GTAACCGA', 'GC146664_AGCTAAGC-GTCGTTAC', 'GC146664_CAACTCCA-GAATCCGT', 'GC146664_CAATAGCC-TCTACGCA', 'GC146664_CAATGCGA-GACGAACT', 'GC146664_CACAGACT-CTGTACCA', 'GC146664_CAGCATAC-GCGTTAGA', 'GC146664_CCACAACA-GTACCACA', 'GC146664_CGCAACTA-TCGGATTC', 'GC146664_CGTATCTC-CAAGGTAC', 'GC146664_CTCGAACA-GCAATTCC', 'GC146664_CTCGACTT-ATCTGACC', 'GC146664_CTTCACTG-CTTAGGAC', 'GC146664_CTTCGGTT-AACCAGAG', 'GC146664_GATAGCCA-CAACACAG', 'GC146664_GATCTTGC-CATGAGCA', 'GC146664_GGACATCA-TACCTGCA', 'GC146664_GGCAAGTT-CTCAGAAG', 'GC146664_GGTGTACA-AGACGCTA', 'GC146664_GGTTGAAC-TCCACGTT', 'GC146664_GTCATCGT-CCGCTTAA', 'GC146664_TGGAAGCA-GAGAGTAC', 'GC146664_TTACGTGC-AGACCTTG', 'GC146665_ATCTGACC-CTCGACTT', 'GC146665_CAGGTTCA-ATAACGCC', 'GC146665_CATGAGCA-GATCTTGC', 'GC146665_CGAATTGC-GTAAGGTG', 'GC146665_CTCAGAAG-GGCAAGTT', 'GC146665_CTTAGGAC-CTTCACTG', 'GC146665_GAACGAAG-GAATCACC', 'GC146665_GAATCCGT-CAACTCCA', 'GC146665_GCAATTCC-CTCGAACA', 'GC146665_GCGTTAGA-CAGCATAC', 'GC146665_GGAAGAGA-CGAGAGAA', 'GC146665_GTCCTAAG-AGCTACCA', 'GC146665_TACTCCAG-GGTATAGG', 'GC146665_TTGCGAGA-TATGGCAC', 'GC146766_AAGACCGT-GTCGATTG', 'GC146766_CAGGTTCA-ATAACGCC', 'GC146766_TACGGTCT-TCCATTGC', 'GC146766_TACTCCAG-GGTATAGG', 'GC146766_TAGGAGCT-GCCTTAAC', 'GC146766_TCGTCTGA-GTCCTTGA', 'GC146766_TTACCGAC-CGAATACG', 'GC146766_TTCCAGGT-CAGTGCTT', 'GC146769_CGAGAGAA-GGAAGAGA', 'GC131878', 'GC131880', 'GC131883', 'GC131887', 'GC131888', 'GC131890', 'GC131891', 'GC131894', 'GC131901', 'GC131902', 'GC131910', 'GC131915', 'GC131916', 'GC142073', 'GC142074', 'GC145416', 'GC145418', 'GC145419', 'GC145422', 'GC145423', 'GC145424', 'GC145425', 'GC145426', 'GC145427', 'GC145432', 'GC145433', 'GC145434', 'GC145435', 'GC145436', 'GC145437', 'GC145439', 'GC146664_GAACGGTT-GGAGGAAT', 'GC146664_GGAATGTC-AGGAACAC', 'GC146664_TGGTGAAG-CTTACAGC', 'GC146665_AGACCTTG-TTACGTGC', 'GC146665_CCAGTATC-AGCGAGAT', 'GC146665_GAGAGTAC-TGGAAGCA', 'GC146665_GTAACCGA-AAGACACC', 'GC146665_GTCGTTAC-AGCTAAGC', 'GC146665_TCTACGCA-CAATAGCC', 'GC146769_AGCGAGAT-CCAGTATC', 'GC146769_ATAACGCC-CAGGTTCA', 'GC146769_CAGTGCTT-TTCCAGGT', 'GC146769_CGAATACG-TTACCGAC', 'GC146769_CTTGCTAG-TCATCTCC', 'GC146769_GAATCACC-GAACGAAG', 'GC146769_GCCTTAAC-TAGGAGCT', 'GC146769_GGTATAGG-TACTCCAG', 'GC146769_GTAAGGTG-CGAATTGC', 'GC146769_GTCCTTGA-TCGTCTGA', 'GC146769_GTCGATTG-AAGACCGT', 'GC146769_TATGGCAC-TTGCGAGA', 'GC146769_TCCATTGC-TACGGTCT', 'GC146769_TCTAGGAG-AGTGACCT', 'GC146769_TGCGTAAC-AGCCTATC', 'GC147297_AAGCGACT-ACGAATCC', 'GC147297_AAGGCGTA-CCATCCGC', 'GC147297_ACAAGGCA-TCGCGCAA', 'GC147297_ACCGGAGT-AAGGAAGG', 'GC147297_ACGAATCC-ATTACCCA', 'GC147297_ATGGCTGT-CCTGTCAA', 'GC147297_CACTGTAG-AAGCGACT', 'GC147297_CCATCCGC-AAGGCGTA', 'GC147297_CTTGACGA-GCACACAA', 'GC147297_GACTTGTG-GTGCACGA', 'GC147297_GTCTGAGT-GACTTGTG', 'GC147297_GTGCACGA-TGATAGGC', 'GC147297_GTTTGCTC-ACCAACAG', 'GC147297_TCCCACGA-TTCAATAG', 'GC147297_TCGCGCAA-ACAAGGCA', 'GC147297_TGATAGGC-GTCTGAGT', 'GC147297_TGTTCGCC-AGGTAGGA', 'GC147297_TTCAATAG-TCCCACGA', 'GC147298_AAGCGACT-ACGAATCC', 'GC147298_ACGAATCC-ATTACCCA', 'GC147298_GTCTGAGT-GACTTGTG']
GC_COLORECTAL_DATASET_TEST = ['GC137557', 'GC140036', 'GC140037', 'GC140038', 'GC140039', 'GC140040', 'GC140041', 'GC140042', 'GC140043', 'GC140044', 'GC140045', 'GC140046', 'GC140047', 'GC140048', 'GC140049', 'GC140050', 'GC140051', 'GC140053', 'GC140054', 'GC140055', 'GC140056', 'GC140057', 'GC140058', 'GC140059', 'GC142059', 'GC142063', 'GC142067', 'GC142071', 'GC142075', 'GC142079', 'GC146381', 'GC146382', 'GC146664_AGCTACCA-GTCCTAAG', 'GC146665_TAGGAGCT-GCCTTAAC', 'GC131877', 'GC131879', 'GC131881', 'GC131882', 'GC131884', 'GC131885', 'GC131886', 'GC131889', 'GC131892', 'GC131893', 'GC131895', 'GC131896', 'GC131897', 'GC131898', 'GC131899', 'GC131900', 'GC131903', 'GC131904', 'GC131905', 'GC131906', 'GC131907', 'GC131908', 'GC131909', 'GC131911', 'GC131912', 'GC131913', 'GC131914', 'GC131917', 'GC131918', 'GC140052', 'GC142056', 'GC142062', 'GC142065', 'GC142076', 'GC142078', 'GC146665_TACGGTCT-TCCATTGC', 'GC147297_ACCAACAG-GTTTGCTC', 'GC147297_AGGTAGGA-TGTTCGCC', 'GC147297_ATTACCCA-CACTGTAG', 'GC147297_CCTGTCAA-ATGGCTGT', 'GC147297_GCACACAA-CTTGACGA', 'GC147298_ATTACCCA-CACTGTAG', 'GC147298_CACTGTAG-AAGCGACT', 'GC147298_GACTTGTG-GTGCACGA', 'GC147298_GTGCACGA-TGATAGGC', 'GC147298_TGATAGGC-GTCTGAGT']
GC_COLORECTAL_DATASET = GC_COLORECTAL_DATASET_TRAIN + GC_COLORECTAL_DATASET_TEST

random.shuffle(GC_BREAST_DATASET)
random.shuffle(GC_COLORECTAL_DATASET)


def end_motif_ohe(motif: str) -> List[int]:
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    features = []
    for c in motif:
        ohe = [0] * 4
        ohe[d[c]] = 1
        features = features + ohe
    return features


def dinucleotide_frequencies(sequence: str) -> List[float]:
    d = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7, 'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
    freqs = np.zeros(16)
    for i in range(len(sequence) - 1):
        di = sequence[i:i+2]
        if di in d:
            freqs[d[di]] += 1
    if np.sum(freqs) > 0:
        freqs /= np.sum(freqs)
    return list(freqs)


def gc_content(sequence: str) -> float:
    n, m = 0., 0.
    for c in sequence:
        if c in 'ACGT':
            n += 1.
        if c in 'GC':
            m += 1.
    if n > 0:
        return m / n
    else:
        return 0.41


ref = None
nucmap = None
def frags_to_tensors(fragments_info: List) -> Tuple[torch.Tensor, torch.Tensor]:
    global ref, nucmap
    if ref is None:
        import pysam
        ref = pysam.FastaFile(REF_GENOME_PATH)
    if nucmap is None:
        nucmap = ReferenceNucMap(nuc_filepath)
    n = len(fragments_info)
    X, w = [], []
    for i in range(len(fragments_info)):
        frag = fragments_info[i]

        chr_id = frag['start_chr_id']
        chr_name = chr_id_to_chr_name(chr_id)
        if chr_name in {'chr0', 'EBV'}:
            chr_name = 'chr1'
        sequence = ref.fetch(chr_name, frag['start_pos'] - 1, frag['end_pos']).upper()

        score_start = nucmap.score(frag['start_chr_id'], frag['start_pos'])
        score_end = nucmap.score(frag['start_chr_id'], frag['end_pos'])
        if (score_start < 0) and (score_end < 0):
            nuc_score = 0
        elif (score_start < 0) and (score_end >= 0):
            nuc_score = score_end
        elif (score_start >= 0) and (score_end < 0):
            nuc_score = score_start
        else:
            nuc_score = 0.5 * (score_start + score_end)

        n_cpg = frag['n_methylated_cpg'] + frag['n_unmethylated_cpg']
        if n_cpg == 0:
            meth_cpg = 0.78
        else:
            meth_cpg = frag['n_methylated_cpg'] / n_cpg

        n_cpg_r1 = frag['n_methylated_cpg_r1'] + frag['n_unmethylated_cpg_r1']
        n_cpg_r2 = frag['n_methylated_cpg_r2'] + frag['n_unmethylated_cpg_r2']
        meth_r1 = 0.78 if (n_cpg_r1 == 0) else (frag['n_methylated_cpg_r1'] / n_cpg_r1)
        meth_r2 = 0.78 if (n_cpg_r2 == 0) else (frag['n_methylated_cpg_r2'] / n_cpg_r2)
        if (n_cpg_r1 == 0) or (n_cpg_r2 == 0):
            jag_index = 0
        else:
            jag_index = (meth_r1 - meth_r2) / np.clip(meth_r1, 0.01, None)

        eps = 1e-5
        features = [
            nuc_score,
            frag['size'] / 100.,
            frag['insert_size'] / 100.,
            (frag['size'] - frag['insert_size']) / 100.,
            meth_cpg,
            np.log1p(n_cpg),
            frag['n_methylated_ch'] / (frag['n_methylated_ch'] + frag['n_unmethylated_ch'] + eps),
            np.log1p(frag['n_methylated_ch'] + frag['n_unmethylated_ch']),
            meth_r1,
            meth_r2,
            np.log1p(n_cpg_r1),
            np.log1p(n_cpg_r2),
            jag_index,
            frag['gc_count'] / (frag['s_count'] + eps),
            frag['strand'] - 1,
            gc_content(sequence),
        ]
        features = features + end_motif_ohe(frag['end_motif_5p'])
        features = features + end_motif_ohe(frag['end_motif_3p'])
        chr_ohe = [0] * 26
        chr_id = (frag['start_chr_id'] - 1 + 26) % 26
        chr_ohe[chr_id] = 1
        features = features + chr_ohe
        #features = features + dinucleotide_frequencies(sequence)

        X.append(features)
        w.append(frag['weight'])
    X = np.asarray(X)
    w = np.asarray(w)
    w *= (len(w) / np.mean(w))
    return torch.FloatTensor(X), torch.FloatTensor(w)


def get_batch(sample, batch_size: int) -> Tuple[torch.Tensor, torch.Tensor]:
    n = len(sample)
    fragments_info = [sample[random.randint(0, n - 1)] for _ in range(batch_size)]
    return frags_to_tensors(fragments_info)


def init_weights(m):
    if isinstance(m, torch.nn.Linear):
        if m.weight.size()[1] == 1:
            torch.nn.init.xavier_uniform_(m.weight)
        else:
            torch.nn.init.kaiming_uniform_(m.weight, nonlinearity='leaky_relu')
        m.bias.data.fill_(1e-3)


class Model(torch.nn.Module):

    def __init__(self, input_size, latent_size=512, dropout=0, n_hidden_layers=3, layer_norm=0, activation='prelu'):
        torch.nn.Module.__init__(self)
        n_out = 1
        layers = []
        layers.append(torch.nn.LayerNorm(input_size))
        #layers.append(Scaler(input_size))
        inplace = False
        for k in range(n_hidden_layers):

            if dropout > 0:
                layers.append(torch.nn.Dropout(p=dropout, inplace=inplace))

            if k == 0:
                layers.append(torch.nn.Linear(input_size, latent_size))
            else:
                layers.append(torch.nn.Linear(latent_size, latent_size))

            if layer_norm:
                layers.append(torch.nn.LayerNorm(latent_size))

            if activation == 'relu':
                layers.append(torch.nn.ReLU(inplace=inplace))
            elif activation == 'leakyrelu':
                layers.append(torch.nn.LeakyReLU(0.2, inplace=inplace))
            elif activation == 'prelu':
                layers.append(torch.nn.PReLU(latent_size))
            elif activation == 'tanh':
                layers.append(torch.nn.Tanh())
            elif activation == 'sigmoid':
                layers.append(torch.nn.Sigmoid())
            elif activation == 'mish':
                layers.append(torch.nn.Mish(inplace=inplace))
            else:
                layers.append(torch.nn.SELU(inplace=inplace))

        layers.append(torch.nn.Linear(latent_size, n_out))

        self.layers = torch.nn.Sequential(*layers)
        self.apply(init_weights)
            
    def forward(self, x):
        return self.layers(x)


def load_model(filepath: str) -> Model:
    n_features = 74
    model = Model(n_features)
    if os.path.exists(filepath):
        model.load_state_dict(torch.load(filepath))
    return model


def train(model_filepath: str, pos_filepath: str, neg_filepath: str, timeout: float = 3600 * 12 * 3) -> None:

    sample_pos = Sample(pos_filepath)
    sample_neg = Sample(neg_filepath)

    model = load_model(model_filepath)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.5 * 1e-4)
    criterion = torch.nn.BCEWithLogitsLoss(reduction='none')

    model.train()
    batch_size = 256

    t0 = time.time()
    while time.time() - t0 < timeout:

        X_pos, w_pos = get_batch(sample_pos, batch_size)
        X_neg, w_neg = get_batch(sample_neg, batch_size)
        X = torch.cat((X_pos, X_neg), dim=0)
        w = 0.5 * torch.cat((w_pos, w_neg), dim=0)
        y = torch.FloatTensor(np.concatenate((np.ones(len(X_pos)), np.zeros(len(X_neg)))))

        optimizer.zero_grad()

        y_hat = model(X).resize(len(X))
        loss = torch.mean(w * criterion(y_hat, y))
        loss.backward()

        accuracy = np.mean(y.cpu().data.numpy() == (y_hat > 0).cpu().data.numpy())
        print(accuracy, loss.item(), flush=True)

        optimizer.step()

        torch.save(model.state_dict(), model_filepath)


def estimate_tf(model: Model, filepath: str, timeout: float = 1000.0) -> List[float]:
    model.eval()

    sample = Sample(filepath)
    n = len(sample)

    scores = np.zeros(101, dtype=float)

    t0 = time.time()
    fragments_info = []
    while True:
        if time.time() - t0 >= timeout:
            break
        i = random.randint(0, n - 1)
        frag = sample[i]
        fragments_info.append(frag)
        if len(fragments_info) == 512:
            X, w = frags_to_tensors(fragments_info)
            w = w.cpu().data.numpy()
            fragments_info = []
            y_hat = model(X).resize(len(X)).cpu().data.numpy()
            y_hat = scipy.special.expit(y_hat)
            idx = np.round(100 * y_hat).astype(int)
            scores[idx] += w
    return [float(x) for x in scores]


def write_fragment(f, frag) -> None:
    byteorder = 'little'
    f.write(frag['start_chr_id'].to_bytes(1, byteorder))
    f.write(frag['start_pos'].to_bytes(4, byteorder))
    f.write(frag['insert_size'].to_bytes(2, byteorder))
    f.write(frag['size'].to_bytes(2, byteorder))
    f.write(frag['n_methylated_cpg_r1'].to_bytes(1, byteorder))
    f.write(frag['n_unmethylated_cpg_r1'].to_bytes(1, byteorder))
    f.write(frag['n_methylated_cpg_r2'].to_bytes(1, byteorder))
    f.write(frag['n_unmethylated_cpg_r2'].to_bytes(1, byteorder))
    f.write(frag['n_methylated_cpg'].to_bytes(1, byteorder))
    f.write(frag['n_unmethylated_cpg'].to_bytes(1, byteorder))
    f.write(frag['n_methylated_chg'].to_bytes(1, byteorder))
    f.write(frag['n_unmethylated_chg'].to_bytes(1, byteorder))
    f.write(frag['n_methylated_chh'].to_bytes(1, byteorder))
    f.write(frag['n_unmethylated_chh'].to_bytes(1, byteorder))
    f.write(frag['gc_count'].to_bytes(2, byteorder))
    f.write(frag['s_count'].to_bytes(2, byteorder))
    f.write(frag['end_motif_5p_int'].to_bytes(1, byteorder))
    f.write(frag['end_motif_3p_int'].to_bytes(1, byteorder))
    f.write(frag['mapq'].to_bytes(1, byteorder))
    f.write(frag['strand'].to_bytes(1, byteorder))


def interpret_sample(model: Model, filepath: str, timeout: float = 1000.0) -> List[int]:
    model.eval()

    sample = Sample(filepath)
    n = len(sample)

    scores = np.zeros(74)

    t0 = time.time()
    fragments_info = []
    while True:
        if time.time() - t0 >= timeout:
            break
        i = random.randint(0, n - 1)
        frag = sample[i]
        fragments_info.append(frag)
        if len(fragments_info) == 512:
            X, w = frags_to_tensors(fragments_info)
            w = w.cpu().data.numpy()
            X.requires_grad_()
            ig = captum.attr.NoiseTunnel(captum.attr.Saliency(model))
            attr = ig.attribute(X, target=0, nt_samples=50, stdevs=0.05)
            scores += np.mean(w[:, np.newaxis] * np.abs(attr.cpu().data.numpy()), axis=0)
            print(list(scores))
    return scores


def train_brca() -> None:
    train(
        os.path.join(ROOT, 'model-brca.pth'),
        os.path.join(DATA_DIR, f'GC110681.modalib'),
        os.path.join(DATA_DIR, f'GC110680.modalib')
    )


def predict_brca() -> None:
    os.makedirs(os.path.join(TF_DIR, 'brca'), exist_ok=True)
    #model = load_model(os.path.join(ROOT, 'model-colca.pth'))
    model = load_model(os.path.join(ROOT, 'model-brca.pth'))
    for gc_code in GC_BREAST_DATASET:
        filepath = os.path.join(DATA_DIR, f'{gc_code}.modalib')
        if os.path.exists(os.path.join(TF_DIR, 'brca', f'{gc_code}.json')):
            continue
        scores = estimate_tf(model, filepath)
        print(gc_code, np.mean(scores))
        with open(os.path.join(TF_DIR, 'brca', f'{gc_code}.json'), 'w') as f:
            data = {'tf': scores}
            json.dump(data, f)


def train_colca() -> None:
    train(
        os.path.join(ROOT, 'model-colca.pth'),
        os.path.join(DATA_DIR, f'GC146664_CAACTCCA-GAATCCGT.modalib'),
        os.path.join(DATA_DIR, f'GC146664_CAATAGCC-TCTACGCA.modalib')
    )


def predict_colca() -> None:
    model = load_model(os.path.join(ROOT, 'model-colca.pth'))
    for gc_code in GC_COLORECTAL_DATASET:
        filepath = os.path.join(DATA_DIR, f'{gc_code}.modalib')
        if os.path.exists(os.path.join(TF_DIR, 'colca', f'{gc_code}.json')):
            continue
        scores = estimate_tf(model, filepath)
        print(gc_code, np.mean(scores))
        with open(os.path.join(TF_DIR, 'colca', f'{gc_code}.json'), 'w') as f:
            data = {'tf': scores}
            json.dump(data, f)


def interpret_brca() -> None:
    model = load_model(os.path.join(ROOT, 'model-brca.pth'))
    weights = model.layers[0].weight.cpu().data.numpy()
    interpret_sample(model, os.path.join(DATA_DIR, f'GC110677.modalib'))


def interpret_colca() -> None:
    os.makedirs(os.path.join(TF_DIR, 'colca'), exist_ok=True)
    model = load_model(os.path.join(ROOT, 'model-colca.pth'))
    weights = model.layers[0].weight.cpu().data.numpy()
    interpret_sample(model, os.path.join(DATA_DIR, f'GC146664_GATAGCCA-CAACACAG.modalib'))


def train_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'train',
        type=str,
        help='train'
    )
    parser.add_argument(
        'pos_path',
        type=str,
        help='Path to modalib file of high-TF cancer cfDNA sample'
    )
    parser.add_argument(
        'neg_path',
        type=str,
        help='Path to modalib file of healthy cfDNA sample'
    )
    parser.add_argument(
        '--o',
        type=str,
        help='Where to store the model',
        default='model.pth'
    )
    parser.add_argument(
        '--timeout',
        type=float,
        help='Training time (the longer, the better)',
        default=3600 * 12 * 3
    )
    args = parser.parse_args()
    train(
        args.o,
        args.pos_path,
        args.neg_path,
        timeout=timeout
    )


def tf_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'tf',
        type=str,
        help='tf'
    )
    parser.add_argument(
        'model_path',
        type=str,
        help='Path to model'
    )
    parser.add_argument(
        'path',
        type=str,
        help='Path to cfDNA sample'
    )
    parser.add_argument(
        '--o',
        type=str,
        help='Path to the output file',
        default='results.json'
    )
    parser.add_argument(
        '--timeout',
        type=float,
        help='TF inference time (the longer, the better)',
        default=1000.0
    )
    parser.add_argument(
        '--cutoff',
        type=float,
        help='Probability cutoff to classify a cfDNA fragment as tumor-derived',
        default=0.95
    )
    args = parser.parse_args()

    model = load_model(args.model_path)
    model.eval()
    scores = estimate_tf(
        model,
        args.path,
        timeout=args.timeout
    )
    print("Estimated TF: ", np.mean(np.asarray(scores) >= args.cutoff))
    with open(args.o, 'w') as f:
        data = {'tf': scores}
        json.dump(data, f)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide a command such as 'train' or 'tf'")
    if sys.argv[1] == 'train':
        train_main()
    elif sys.argv[1] == 'tf':
        tf_main()
    elif sys.argv[1] == 'predict-brca':
        predict_brca()
    elif sys.argv[1] == 'train-brca':
        train_brca()
    elif sys.argv[1] == 'predict-colca':
        predict_colca()
    elif sys.argv[1] == 'train-colca':
        train_colca()
    elif sys.argv[1] == 'interpret-colca':
        interpret_colca()
    elif sys.argv[1] == 'interpret-brca':
        interpret_brca()
    else:
        raise NotImplementedError(f'Unknown command "{sys.argv[1]}"')
