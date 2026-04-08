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