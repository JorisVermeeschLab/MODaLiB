from typing import Tuple
import os
import re
import tqdm
import xml.etree.ElementTree as ET

import numpy as np


def find_most_reliable_base(r1: np.ndarray, r2: np.ndarray, robust: bool = False) -> int:
    if robust:
        window_size = 9
        diff = np.abs(r1 - r2)
        diff = np.convolve(diff, np.ones(window_size), 'valid') / window_size  # Moving average
        return np.argmin(diff) + (window_size // 2)
    else:
        return len(r1) // 2


def compute_bounds_(curve: np.ndarray, i: int, alpha: float=0.5, beta: float=2.0, gamma: float=0.07, clip: int=0) -> Tuple[float, float]:

    n = len(curve)
    lb = ub = i
    beta = curve[i]
    auc = np.mean(curve[lb:ub+1])

    C = np.zeros((n, n))
    C2 = np.zeros((n, n))
    C3 = np.zeros((n, n))
    for ub in range(n):
        for lb in range(ub + 1):
            cost = abs(np.mean(curve[lb:ub+1]) - beta)
            C[lb, ub] = cost
            C[ub, lb] = cost
            C2[lb, ub] = abs(ub - lb + 1)
            C2[ub, lb] = abs(ub - lb + 1)
            if lb == ub:
                cost = 0
            else:
                cost = np.max(np.abs(curve[lb:ub] - curve[lb+1:ub+1]))
            C3[lb, ub] = cost
            C3[ub, lb] = cost

            if (lb < clip) or (ub >= n - clip) or (np.max(np.abs(curve[lb:ub+1] - beta)) > gamma):
                C[lb, ub] = np.inf
                C[ub, lb] = np.inf

    C = C + alpha * (1 - C2 / n) + beta * C3

    lb, ub = np.unravel_index(np.argmin(C), C.shape)
    if ub < lb:
        lb, ub = ub, lb

    return lb, ub


def parse_path(s: str) -> np.ndarray:
    path = []
    elements = s.split()
    i = 0
    x, y = 0, 0
    while i < len(elements):
        cmd = elements[i]
        if cmd in {'M', 'L'}:
            x = float(elements[i + 1])
            y = float(elements[i + 2])
            i += 3
        elif cmd in {'m', 'l'}:
            x += float(elements[i + 1])
            y += float(elements[i + 2])
            i += 3
        elif cmd == 'Z':
            x, y = path[0]
            i += 1
        elif cmd == 'H':
            x = float(elements[i + 1])
            i += 2
        elif cmd == 'h':
            x += float(elements[i + 1])
            i += 2
        elif cmd == 'V':
            y = float(elements[i + 1])
            i += 2
        elif cmd == 'v':
            y += float(elements[i + 1])
            i += 2
        else:
            raise NotImplementedError(f'Unknown svg path command "{cmd}"')
        path.append([x, y])
    return np.asarray(path)


def is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def compute_bounds(filepath: str, **kwargs) -> str:
    tree = ET.parse(filepath)
    root = tree.getroot()

    lb_coord, ub_coord = None, None
    lb, ub = 1.0, 0.0

    paths = []
    for element in root:
        if 'path' in element.tag:
            coords = parse_path(element.attrib['d'])
            paths.append(coords)
        elif 'text' in element.tag:
            if is_float(element.text) and ('.' in element.text):
                if float(element.text) < lb:
                    lb_coord = float(element.attrib['y'])
                    lb = float(element.text)
                elif float(element.text) > ub:
                    ub_coord = float(element.attrib['y'])
                    ub = float(element.text)
    assert lb < ub

    n_bases = np.min([len(path) for path in paths])
    paths = [path for path in paths if (len(path) == n_bases)]
    assert len(paths) == 2

    for i in range(len(paths)):
        idx = np.argsort(paths[i][:, 0])
        paths[i] = paths[i][idx, 1]
        paths[i] = ((paths[i] - lb_coord) / (ub_coord - lb_coord)) * (ub - lb) + lb

    i = find_most_reliable_base(paths[0], paths[1])
    lb1, ub1 = compute_bounds_(paths[0], i, **kwargs)
    lb2, ub2 = compute_bounds_(paths[1], i, **kwargs)

    lb1 += 1
    lb2 += 1
    ub1 += 1
    ub2 += 1

    if lb1 == 1:
        lb1 = 0
    if lb2 == 1:
        lb2 = 0
    if ub1 == len(paths[0]):
        ub1 = 0
    if ub2 == len(paths[1]):
        ub2 = 0

    bounds = [int(lb1), int(ub1), int(lb2), int(ub2)]

    paths = [
        [float(x) for x in paths[0]],
        [float(x) for x in paths[1]]
    ]

    return paths, bounds
