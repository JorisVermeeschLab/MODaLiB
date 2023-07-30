Compute features and kernel matrices
====================================

## Installation

First install Python dependencies:
```bash
cd ..
pip3 install -r requirements.txt
```

Compile the C++ and Cython code:
```bash
cd step2-fragment-level-analysis
python3 setup.py install --user
```

---

## Usage

### Create a reference nucleosome map



### Extract features

```python
from modalib.core import Sample, extract_features

samples = [
    Sample('path/to/sample/1'),
    Sample('path/to/sample/2'),
    Sample('path/to/sample/3')
]
nuc_folder = '/lustre1/project/stg_00019/research/Antoine/MODaLiB/step2-fragment-level-analysis/nucleosomes'
matrices = extract_features(samples, nuc_folder)
```

`nuc_folder` is a folder 

### Compute fragment-level kernel matrices

```python
from modalib.core import Sample, pairwise_kernels
```
