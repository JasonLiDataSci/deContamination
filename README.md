# deContamination

Python bindings for the `Contamination` C++ implementation, built with `pybind11`, `CMake`, and OpenMP.

## Overview

This package exposes the C++ implementation as a Python module named `deContamination`.
It is intended for use in scripted workflows and batch experiments without manual compilation.

## Installation

Install the wheel file with `pip`:

```bash
python3 -m pip install decontamination-0.1.0-cp310-cp310-linux_x86_64.whl
```

After installation, verify the module:

```bash
python3 -c "import deContamination; print(deContamination.__file__)"
```

## Usage

A minimal test script is shown below.

```python
# test_decontamination.py
from pathlib import Path
import traceback
import deContamination


def main():
    out_dir = Path("./output/")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loaded module from:", deContamination.__file__)
    print("Available attributes:")
    print([x for x in dir(deContamination) if not x.startswith("_")])

    if hasattr(deContamination, "run"):
        run = deContamination.run
    elif hasattr(deContamination, "run_contamination"):
        run = deContamination.run_contamination
    else:
        raise AttributeError("No callable function found in deContamination.")

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
        "data_name": "../../dcuda/data/Y_obs.txt",
        "nei_name": "../../dcuda/data/nei_list.txt",
        "dist_name": "../../dcuda/data/nei_dist.txt",
        "label_name": "../../dcuda/data/Y_label.txt",
        "cell_size_name": "../../dcuda/data/cell_size.txt",
        "MB_dir": "../../dcuda/data/MB/",
    }

    print("Start running...")
    try:
        ret = run(
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

Run the script with:

```bash
python3 test_decontamination.py
```

## Input files

The package expects the same input files used by the original C++ workflow.
Update the paths in the Python script to match your local directory layout.

Typical inputs include:

* observed data file
* neighbor list file
* neighbor distance file
* label file
* cell size file
* MB directory
* true Z file, if enabled in your binding

## Output

Results are written to the output directory passed to the Python function.
In the example above, all outputs are written under `./output/`.

## Notes

* The module name is `deContamination`.
* The callable function name depends on the pybind11 binding.
* The argument order must match the C++ binding exactly.
* If the wheel is rebuilt for another Python version or operating system, a new wheel is required.

## Troubleshooting

### `ModuleNotFoundError: No module named 'deContamination'`

Check that the wheel was installed in the same Python environment used to run the script.

### `No callable function found`

Check the function name exported in the `PYBIND11_MODULE(...)` binding.

### Output files are missing

Confirm that the output directory exists and that the input paths are correct.

## License

Add your license information here.
