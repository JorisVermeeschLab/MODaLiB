from typing import Tuple, List, Optional

import numpy as np
import scipy.optimize
import torch
import torch.nn.functional
import tqdm
from sklearn.preprocessing import RobustScaler

from .optimizer import Optimizer


def pearson_corrcoef(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:

    #from fast_soft_sort.pytorch_ops import soft_rank, soft_sort
    #x = soft_rank(x.unsqueeze(0)).squeeze()
    #y = soft_rank(y.unsqueeze(0)).squeeze()

    vx = x - torch.mean(x)
    vy = y - torch.mean(y)
    return torch.mean(vx * vy) / (torch.sqrt(torch.mean(vx ** 2)) * torch.sqrt(torch.mean(vy ** 2)))


class MLP(torch.nn.Module):

    def __init__(self, shape: List[int]):
        torch.nn.Module.__init__(self)
        self.shape: List[int] = shape
        assert len(self.shape) > 1

        # Define layers
        layers = torch.nn.Sequential()
        layers.append(torch.nn.LayerNorm(self.shape[0]))
        for i in range(len(self.shape) - 1):
            layers.append(torch.nn.Linear(self.shape[i], self.shape[i + 1]))
            if i < len(self.shape) - 2:
                layers.append(torch.nn.LeakyReLU())
                #layers.append(torch.nn.BatchNorm1d(self.shape[i + 1]))
        #layers.append(torch.nn.ReLU())
        self.layers = torch.nn.Sequential(*layers)

        # Xavier initialization
        def init_weights(m):
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    m.bias.data.fill_(0.001)
        self.layers.apply(init_weights)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        return self.layers.forward(X)


def interdomain_deconvolution(
        X: np.ndarray,
        R: np.ndarray,
        model: Optional[torch.nn.Module] = None
) -> Tuple[np.ndarray, torch.nn.Module]:

    #R = RobustScaler().fit_transform(R)

    n_in_channels = X.shape[1]
    n_hidden_channels = 64
    n_out_channels = 1

    X = torch.FloatTensor(X)
    R = torch.FloatTensor(R)

    m = len(R)
    p = len(X)

    # Initialization
    alpha = np.random.rand(m) + 1e-7
    alpha /= np.sum(alpha)
    alpha_logit = torch.nn.Parameter(torch.FloatTensor(np.log(alpha)))

    optimizer = Optimizer(verbose=False)
    optimizer.add([alpha_logit], 1e-2)
    train_nn = (model is None)
    if train_nn:
        model = MLP([n_in_channels, n_hidden_channels, n_hidden_channels, n_out_channels])
        optimizer.add(model.parameters(), 1e-4)
        model.train()
    else:
        model.eval()

    with torch.no_grad():
        x = model(X).squeeze()

    if not train_nn:
        x = x.cpu().data.numpy()
        R = R.cpu().data.numpy()
        alpha, _ = scipy.optimize.nnls(R.T, x)
        alpha /= np.sum(alpha)
        return alpha, model

    patience = 30
    n_steps_without_improvement = 0
    best_loss = np.inf
    max_n_iter = 5000
    pbar = tqdm.tqdm(range(max_n_iter))
    for iteration in pbar:

        optimizer.zero_grad()

        if train_nn:
            x = model(X).squeeze()

        alpha = torch.softmax(alpha_logit, dim=0)

        x_hat = torch.mv(R.t(), alpha)

        residuals = torch.square(x - x_hat)

        loss = torch.mean(residuals)
        #loss = -torch.abs(pearson_corrcoef(x, x_hat))
        rho = torch.abs(pearson_corrcoef(x, x_hat)).item()
        loss.backward()
        pbar.set_description(str(rho))

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

    alpha = torch.softmax(alpha_logit, dim=0).cpu().data.numpy()

    model.eval()

    return alpha, model
