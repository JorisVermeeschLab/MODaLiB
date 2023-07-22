import random

import numpy as np
import matplotlib.pyplot as plt
import seaborn

LOG_SCALE = False

X_0_1x = np.load('pvalues-0.1x.npy')
X_1x = np.load('pvalues-1x.npy')
X_10x = np.load('pvalues-10x.npy')

fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
colors = ['#264653', '#2a9d8f', '#f4a261', '#e76f51']

r = ax[0].violinplot(list(X_10x), vert=False, showmeans=True, showextrema=True)
color = 'black'
r['cbars'].set_color(color)
r['cbars'].set_linewidth(0.2)
r['cmins'].set_color(color)
r['cmaxes'].set_color(color)
r['cmeans'].set_color(color)
for body, color in zip(r['bodies'], colors):
    body.set_color(color)
# ax[0, 3].set_xscale('log')
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[0].set_xlabel('Two-sample KS p-value (10x)')
ax[0].set_yticks(range(1, 5), ['Nucleosomes', 'Lengths', 'Methylation', 'Mutations'])
if LOG_SCALE:
    ax[0].set_xscale('log')
ax[0].axvline(x=0.05, color='black', linestyle='--', linewidth=0.5)

r = ax[1].violinplot(list(X_1x), vert=False, showmeans=True, showextrema=True)
color = 'black'
r['cbars'].set_color(color)
r['cbars'].set_linewidth(0.2)
r['cmins'].set_color(color)
r['cmaxes'].set_color(color)
r['cmeans'].set_color(color)
for body, color in zip(r['bodies'], colors):
    body.set_color(color)
# ax[0, 3].set_xscale('log')
ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].set_xlabel('Two-sample KS p-value (1x)')
ax[1].set_yticks(range(1, 5), ['Nucleosomes', 'Lengths', 'Methylation', 'Mutations'])
if LOG_SCALE:
    ax[1].set_xscale('log')
ax[1].axvline(x=0.05, color='black', linestyle='--', linewidth=0.5)

r = ax[2].violinplot(list(X_0_1x), vert=False, showmeans=True, showextrema=True)
color = 'black'
r['cbars'].set_color(color)
r['cbars'].set_linewidth(0.2)
r['cmins'].set_color(color)
r['cmaxes'].set_color(color)
r['cmeans'].set_color(color)
for body, color in zip(r['bodies'], colors):
    body.set_color(color)
# ax[0, 3].set_xscale('log')
ax[2].spines['right'].set_visible(False)
ax[2].spines['top'].set_visible(False)
ax[2].set_xlabel('Two-sample KS p-value (0.1x)')
ax[2].set_yticks(range(1, 5), ['Nucleosomes', 'Lengths', 'Methylation', 'Mutations'])
if LOG_SCALE:
    ax[1].set_xscale('log')
ax[2].axvline(x=0.05, color='black', linestyle='--', linewidth=0.5)

plt.tight_layout()
if LOG_SCALE:
    plt.savefig('ks-log.png', dpi=200)
else:
    plt.savefig('ks.png', dpi=200)
