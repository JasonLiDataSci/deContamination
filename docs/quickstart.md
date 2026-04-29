# <span style="color: black;">**Getting Started with DeLeakage**</span>

This tutorial walks you through loading simulation data, running DeLeakage, and visualizing results.


## 1. Installation
First, ensure the `deContamination` module is installed. If using a wheel:
```bash
pip install deContamination-*.whl
```


## 2. Load Simulation Data
To inspect the generated data (optional):
```python
import numpy as np

# Load observed data
Y_obs = np.loadtxt("./data/Y_obs.txt", dtype=int)
print("Y_obs shape:", Y_obs.shape)  # Should be (300, 10000)

# Load neighbor list
nei_list = np.loadtxt("./data/nei_list.txt", dtype=int)
print("nei_list shape:", nei_list.shape)  # Should be (10000, 49)
```


## 3. Run DeLeakage
Use the script `test_decontamination.py` (adapted from your example):
```python
from pathlib import Path
import traceback
import deContamination as dc


def main():
    out_dir = Path("./output/")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loaded module from:", dc.__file__)
    print("Available attributes:", [x for x in dir(dc) if not x.startswith("_")])

    # Parameters (match simulation data)
    params = {
        "G": 300,
        "K": 3,
        "N_nei": 49,
        "N": 10000,
        "N_MB": 50,
        "N_tail": 454,
        "n_record": 3000,
        "seed": 123,
        "output_list": "./output/",
        "data_name": "./data/Y_obs.txt",
        "nei_name": "./data/nei_list.txt",
        "dist_name": "./data/nei_dist.txt",
        "label_name": "./data/Y_label.txt",
        "cell_size_name": "./data/cell_size.txt",
        "MB_dir": "./data/MB/",
    }

    print("Start running...")
    try:
        ret = dc.run(
            params["G"],
            params["K"],
            params["N_nei"],
            params["N"],
            params["N_MB"],
            params["N_tail"],
            params["n_record"],
            params["seed"],
            params["output_list"],
            params["data_name"],
            params["nei_name"],
            params["dist_name"],
            params["label_name"],
            params["cell_size_name"],
            params["MB_dir"],
        )
        print("Finished. Return value:", ret)
    except Exception:
        print("Run failed with exception:")
        traceback.print_exc()


if __name__ == "__main__":
    main()
```

Run it:
```bash
python3 test_decontamination.py
```


## 4. Visualize Results
Using `plot_results.py` to compare true, raw vs. decontaminated data:

Run it:
```bash
python .\scripts\plot_results.py --gex .\True\Y_true.txt --gene .\dataDir\Gene.csv --cell .\dataDir\Cell_meta.csv --marker .\dataDir\marker.json --output .\eval_True

python .\scripts\plot_results.py --gex .\dataDir\Count.mtx --gene .\dataDir\Gene.csv --cell .\dataDir\Cell_meta.csv --marker .\dataDir\marker.json --output .\eval_Raw

python .\scripts\plot_results.py --gex .\our\decon_fact_vol.txt --gene .\dataDir\Gene.csv --cell .\dataDir\Cell_meta.csv --marker .\dataDir\marker.json --output .\eval_our

```


## 5. Input/Output Summary
### Inputs (./data/)
- `Y_obs.txt`: Observed expression matrix (genes x cells)
- `nei_list.txt`: Neighbor indices for each cell
- `nei_dist.txt`: Neighbor distances
- `Y_label.txt`: Cell type labels
- `cell_size.txt`: Cell size estimates
- `MB/`: Microball data directory

### Outputs (./output/)
- Decontaminated expression matrix (e.g., `Y_pred.txt`)
- Log files/other intermediate results


## 7. Notes
- Module name: `deContamination`
- Argument order must match the C++ binding exactly
- If rebuilding the wheel for a new Python/OS version, reinstall the updated wheel

