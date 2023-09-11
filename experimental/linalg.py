import numpy as np
from sklearn.decomposition import KernelPCA


n_features = 30
X = np.random.rand(20, n_features)
X *= 7 * (np.random.rand(1, n_features) - 0.5)
K = X @ X.T

pcs = KernelPCA(kernel='precomputed').fit_transform(K)
variances = np.var(pcs, axis=0)
print(f'Total variance: {np.sum(variances)}')
eigenvalues = np.linalg.eigvals(K)
print(f'Total variance (eig): {np.sum(eigenvalues)}')

n = len(K)
V = np.eye(n) - (1. / n)
K = V @ K @ V
eigenvalues = np.linalg.eigvals(K)
K /= np.sum(eigenvalues)

pcs = KernelPCA(kernel='precomputed').fit_transform(K)
print(pcs.shape, K.shape)
variances = np.var(pcs, axis=0)

print(f'Total variance: {np.sum(variances)}')
eigenvalues = np.linalg.eigvals(K)
print(f'Total variance (eig): {np.sum(eigenvalues)}')
