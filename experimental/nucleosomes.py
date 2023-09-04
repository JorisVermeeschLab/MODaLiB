from typing import List, Tuple

import numpy as np
import scipy.special
import matplotlib.pyplot as plt


def create_ref_nucleosome_map(n_nucleosomes: int) -> np.ndarray:
	nuc_distances = np.random.randint(120, 200, size=n_nucleosomes)
	nuc_distances[0] = 0
	nuc_distances[5] = 1000
	nuc_positions = np.cumsum(nuc_distances)
	return nuc_positions


def create_pdf_from_positions(nuc_positions: np.ndarray) -> np.ndarray:
	pdf = np.zeros(nuc_positions[-1])
	for i in range(1, len(nuc_positions)):
		for t in range(nuc_positions[i - 1], nuc_positions[i]):
			pdf[t] = min(abs(nuc_positions[i - 1] - t), abs(nuc_positions[i] - t)) / (0.5 * (nuc_positions[i] - nuc_positions[i - 1]))
	pdf = np.maximum(pdf, 0.1)
	pdf /= np.sum(pdf)
	return pdf


def sample_fragments(pdf: np.ndarray, n_fragments: int) -> List[Tuple[int, int]]:
	idx = np.arange(len(pdf))
	starts = np.random.choice(idx, p=pdf, size=n_fragments)
	ends = np.empty_like(starts)
	for i in range(n_fragments):
		try:
			ends[i] = starts[i] + 80 + np.random.choice(
				np.arange(120),
				p=pdf[starts[i]+80:starts[i] + 200] / np.sum(pdf[starts[i]+80:starts[i] + 200])
			)
		except ValueError:
			ends[i] = starts[i]
			starts[i] = ends[i] - 80 - np.random.choice(
				np.arange(120),
				p=pdf[ends[i]-200:ends[i]-80] / np.sum(pdf[ends[i]-200:ends[i]-80])
			)
	return [(start, end) for start, end in zip(starts, ends)]


def wps(genome_size: int, fragments: List[Tuple[int, int]]) -> np.ndarray:
	window_size = 60
	hw = window_size // 2
	scores = np.zeros(genome_size)
	spanning = np.zeros(genome_size)
	ending = np.zeros(genome_size)
	for start, end in fragments:
		spanning[start+hw:end-hw+1] += 1
		ending[start-hw:start+hw] += 1
		ending[end-hw+1:end+hw+1] += 1
	scores = (spanning - ending) / (spanning + ending)
	return scores


nuc_positions = create_ref_nucleosome_map(20)
pdf = create_pdf_from_positions(nuc_positions)

fragments = sample_fragments(pdf, 100000)

scores = wps(len(pdf), fragments)

for t in nuc_positions:
	plt.axvline(x=t)
plt.plot(scores)
plt.show()
