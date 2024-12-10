from operator import index
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import os
import statistics
import sys


def get_size(file_path, unit="bytes", r=2):
    file_size = os.path.getsize(file_path)
    exponents_map = {"bytes": 0, "kb": 1, "mb": 2, "gb": 3}
    if unit not in exponents_map:
        raise ValueError(
            "Must select from \
        ['bytes', 'kb', 'mb', 'gb']"
        )
    else:
        size = file_size / 1024 ** exponents_map[unit]
        return round(size, r)


def get_cn(s):
    t = s.split("_")
    if t[1] == "0":
        return "gindex without cache"
    else:
        return f"gindex with cache of length {t[1]}"


def main(argv):
    n_q = argv[3]
    p_i = argv[4]
    kb_to_gb = 9.5367431640625e-7
    if os.path.exists("~/.matplotlib/stylelib/nord-light.mplstyle"):
        plt.style.use("~/.matplotlib/stylelib/nord-light.mplstyle")
    index_df = pd.read_csv(argv[0])
    index_df["max_mem"] = index_df["max_mem"] * kb_to_gb
    query_df = pd.read_csv(argv[1])
    query_df["max_mem"] = query_df["max_mem"] * kb_to_gb
    out_pref = argv[2]
    if out_pref[-1] == "/":
        out_pref = out_pref[0:-1]
    ## print(query_df)
    ## print(index_df)
    df = query_df
    plt.figure(figsize=(10, 6))
    for tool in df["tool"].unique():
        subset = df[df["tool"] == tool]
        plt.plot(subset["l_query"], subset["wall_clock"], marker="o", label=tool)

    plt.xlabel("Query length")
    plt.ylabel("Time (seconds)")
    plt.title(f"Querying time ({n_q} queries) varying cache length")
    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels = [get_cn(label) for label in labels]

    # plt.legend(title="Tool")
    plt.tight_layout()
    plt.legend(
        handles, new_labels, loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=4
    )

    plt.savefig(f"{out_pref}/query_time.pdf", dpi=500, bbox_inches="tight")
    df = query_df
    plt.figure(figsize=(10, 6))
    for tool in df["tool"].unique():
        subset = df[df["tool"] == tool]
        plt.plot(subset["l_query"], subset["max_mem"], marker="o", label=tool)

    plt.xlabel("Query length")
    plt.ylabel("Memory (gigabytes)")
    plt.title(f"Querying memory usage ({n_q} queries) varying cache length")
    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels = [get_cn(label) for label in labels]
    plt.tight_layout()
    plt.legend(
        handles, new_labels, loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=4
    )
    # plt.legend(title="Tool")
    plt.savefig(f"{out_pref}/query_mem.pdf", dpi=500, bbox_inches="tight")

    df = index_df
    plt.figure(figsize=(10, 6))
    for tool in df["tool"].unique():
        subset = df[df["tool"] == tool]
        plt.scatter(subset["wall_clock"], subset["max_mem"], label=tool, s=100)

    plt.xlabel("Time (seconds)")
    plt.ylabel("Memory (gigabytes)")
    plt.title("Indexing performances varying cache length")
    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels = [get_cn(label) for label in labels]

    plt.tight_layout()
    plt.legend(
        handles, new_labels, loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=4
    )
    # plt.legend(handles, new_labels, title="Tool")
    # plt.legend(title="Tool")
    plt.savefig(f"{out_pref}/index.pdf", dpi=500, bbox_inches="tight")

    print(p_i)
    name_files = []
    size_files = []
    for root, dirs, files in os.walk(p_i):
        for file in files:
            if file.endswith(".cache"):
                file_path = os.path.join(root, file)
                ca = file.split(".")[0]
                name_files.append(f"{ca}")
                size_files.append(get_size(file_path, "gb", 3))
    print(name_files)
    print(size_files)
    with open(f"{out_pref}/size.tex", "w") as f:
        f.write(f"Cache & Size\\\\\n\\hline\n")
        for i in range(len(size_files)):
            f.write(f"{name_files[i]} & {size_files[i]}\\\\\n")


if __name__ == "__main__":
    main(sys.argv[1:])
