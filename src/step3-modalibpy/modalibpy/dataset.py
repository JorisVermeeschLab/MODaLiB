from typing import List, Optional

import numpy as np

from .io import load_profiles
from .stages import Annotation
from .modalities.base import MultimodalData


class Dataset(object):

    def __init__(
            self,
            filepaths: List[str],
            annotations: List[Annotation],
            groups: np.ndarray,
            data: Optional[MultimodalData] = None
    ) -> None:
        self.filepaths: List[str] = list(filepaths)
        self.annotations: List[Annotation] = list(annotations)
        self.groups: np.ndarray = groups

        # TODO: don't load if data is too large
        self.data: MultimodalData = load_profiles(filepaths) if (data is None) else data

    @staticmethod
    def subset_data(data: MultimodalData, key: np.ndarray) -> MultimodalData:
        blacklist = {'CpG Methylation (deconvolution, reference)', 'Gene-expression-reference'}
        new_data = {}
        for mod in data.keys():
            if (mod in blacklist) or (not isinstance(data[mod], np.ndarray) or (len(data[mod].shape) != 2)):
                new_data[mod] = data[mod]
            else:
                new_data[mod] = data[mod][key]
        return new_data

    def __getitem__(self, key: np.ndarray) ->'Dataset':
        return Dataset(
            list(np.asarray(self.filepaths, dtype=object)[key]),
            list(np.asarray(self.annotations, dtype=object)[key]),
            self.groups[key],
            data=Dataset.subset_data(self.data, key)
        )

    def __len__(self) -> int:
        return len(self.filepaths)

    def __add__(self, other: 'Dataset') -> 'Dataset':

        keys = set(self.data.keys()).intersection(set(other.data.keys()))
        if len(keys) < max(len(set(self.data.keys())), len(set(other.data.keys()))):
            print('Warning: Some modalities are incompatible between datasets')
        data = {}
        for key in keys:
            data[key] = np.concatenate((self.data[key], other.data[key]), axis=0)

        return Dataset(
            self.filepaths + other.filepaths,
            self.annotations + other.annotations,
            np.concatenate((self.groups, other.groups), axis=0),
            data=data
        )
