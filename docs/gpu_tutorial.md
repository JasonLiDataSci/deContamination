
# <span style="color: black;">**GPU Implementation of DeLeakage**</span>

GPU‑accelerated MCMC sampling library for gene expression data, wrapped in Python.

## Overview

**DeLeakage** provides a high‑performance CUDA backend for graph‑based MCMC sampling on gene expression data.  
A simple Python API lets you run large‑scale MCMC jobs on your GPU without worrying about CUDA details.

---

## Download & Install

You don’t need PyPI—just grab the entire repository and install locally.

### 1. Clone or Download

```bash
# Option A: clone with Git
git clone https://github.com/JasonLiDataSci/cuda_mcmc.git
cd cuda_mcmc

# Option B: download ZIP from GitHub
1. Visit https://github.com/JasonLiDataSci/cuda_mcmc
2. Click “Code” → “Download ZIP”
3. Unzip and cd into the folder:
unzip cuda_mcmc-main.zip
cd cuda_mcmc-main
````

### 2. Build the CUDA Library

Make sure you have CUDA Toolkit (≥11.0) installed and in your `PATH`.

```bash
# Compile the shared library
make build
```

> This runs `nvcc -O2 -std=c++14 -fPIC -shared mcmc.cu -o cuda_mcmc/libmcmc.so`.

If you prefer without Makefile:

```bash
nvcc -O2 -std=c++14 --compiler-options '-fPIC' -shared mcmc.cu -o cuda_mcmc/libmcmc.so
```

### 3. Install the Python Package

It’s best to use a virtual environment:

```bash
# Create & activate venv
python3 -m venv venv
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install locally
pip install .
```

> This will install the `cuda_mcmc` package and include the compiled `libmcmc.so`.

---

## Quick Start

1. Prepare your data files in the project root:

   * `Y_obs.txt`       (G×N expression matrix)
   * `nei_list.txt`    (neighbor list)
   * `cell_size.txt`   (cell size vector)
   * `Y_label.txt`     (cell labels)

2. Edit or use the provided `example_run.py`:

   ```python
   # example_run.py
   import os
   import cuda_mcmc

   def main():
       outdir = "./output/"
       os.makedirs(outdir, exist_ok=True)

       ret = cuda_mcmc.run_mcmc(
           g=300,      # feature dimension G
           k=3,        # number of cell types K
           b=0,        # number of neighbors B
           n=3000,     # total iterations N
           a=20,       # minibatch size
           c=78,       # tail sample size
           t=3000,     # record frequency / burn‑in
           s=666,      # random seed
           data_file="./data/Y_obs.txt",
           nei_file="./data/nei_list.txt",
           label_file="./data/Y_label.txt",
           cell_size="./data/cell_size.txt",
           mb_dir="./data/MB/",
           output_dir=outdir
       )

       if ret == 0:
           print("MCMC completed successfully. Results in", outdir)
       else:
           print("MCMC failed with return code", ret)

   if __name__ == "__main__":
       main()
   ```

3. Run it:

   ```bash
   python example_run.py
   ```

---

## API Reference

```python
cuda_mcmc.run_mcmc(
    g: int,          # feature count
    k: int,          # cell type count
    b: int,          # neighbor count
    n: int,          # total iterations
    a: int,          # minibatch size
    c: int,          # tail sample size
    t: int,          # record frequency / burn‑in
    s: int,          # random seed
    data_file: str,  # path to expression matrix
    nei_file: str,   # path to neighbor list
    label_file: str, # path to cell labels
    cell_size: str,  # path to cell size file
    mb_dir: str,     # intermediate files dir
    output_dir: str  # output directory
) → int
```

* Returns `0` on success, non‑zero on failure.

---



