from typing import List, Optional, Callable

import numpy as np
import torch
import pyro


def svi_optimize(
        args: List[np.ndarray],
        model: Callable,
        guide: Optional[Callable] = None,
        lr: float = 0.005,
        n_steps: int = 201,
        verbose: bool = True
) -> None:

    if guide is None:
        def guide_mle(*args):
            pass
        guide = guide_mle

    args = [torch.FloatTensor(arg) for arg in args]

    pyro.clear_param_store()
    adam = pyro.optim.SGD({'lr': lr})
    svi = pyro.infer.SVI(model, guide, adam, loss=pyro.infer.Trace_ELBO())

    for step in range(n_steps):
        loss = svi.step(*args)
        if verbose and (step % 1 == 0):
            print('[iter {}]  loss: {:.4f}'.format(step, loss))
