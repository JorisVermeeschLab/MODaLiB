from typing import Tuple

import numpy as np
import tqdm
import numpy as np
import scipy.optimize
import scipy.special
from scipy.optimize import nnls
from scipy.special import logit
import torch
import torch.nn.functional
from sklearn.decomposition import NMF

from modalibpy.deconvolution.optimizer import Optimizer


def unsupervised_multinomial_deconvolution(X: np.ndarray, m: int, max_n_iter: int = 2000) -> Tuple[np.ndarray, np.ndarray]:

    lr = 1e-4

    X = torch.FloatTensor(X / np.sum(X, axis=1)[:, np.newaxis])

    # Initialization
    model = NMF(n_components=m)
    W = model.fit_transform(X)
    H = model.components_
    W = W + 1e-15
    H = H + 1e-15
    W_logit = torch.nn.Parameter(torch.FloatTensor(np.log(W / np.sum(W, axis=1)[:, np.newaxis])))
    H_logit = torch.nn.Parameter(torch.FloatTensor(np.log(H / np.sum(H, axis=1)[:, np.newaxis])))

    optimizer = Optimizer(verbose=False)
    optimizer.add([W_logit], lr)
    optimizer.add([H_logit], lr)

    patience = 50
    n_steps_without_improvement = 0
    best_loss = np.inf
    pbar = tqdm.tqdm(range(max_n_iter))
    for iteration in pbar:

        optimizer.zero_grad()

        W = torch.softmax(W_logit, dim=1)
        H = torch.softmax(H_logit, dim=1)

        loss = torch.mean(torch.square(X - torch.mm(W, H)))
        loss.backward()
        pbar.set_description(str(loss))

        # Update parameters
        optimizer.step(loss.item())

        # scheduler.step()
        if loss.item() >= best_loss:
            n_steps_without_improvement += 1
            if n_steps_without_improvement >= patience:
                break
        else:
            best_loss = loss.item()
            n_steps_without_improvement = 0

    W = torch.softmax(W_logit, dim=1).cpu().data.numpy()
    H = torch.softmax(H_logit, dim=1).cpu().data.numpy()

    return W, H
