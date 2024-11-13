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


def main(argv):
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

    # Etichette e titolo
    plt.xlabel("l_query")
    plt.ylabel("wall_clock")
    plt.title("Query wall_clock")
    plt.legend(title="Tool")
    plt.savefig(f"{out_pref}/query_time.pdf", dpi=500)
    df = query_df
    plt.figure(figsize=(10, 6))
    for tool in df["tool"].unique():
        subset = df[df["tool"] == tool]
        plt.plot(subset["l_query"], subset["max_mem"], marker="o", label=tool)

    # Etichette e titolo
    plt.xlabel("l_query")
    plt.ylabel("max_mem")
    plt.title("Query max_mem")
    plt.legend(title="Tool")
    plt.savefig(f"{out_pref}/query_mem.pdf", dpi=500)

    df = index_df
    plt.figure(figsize=(10, 6))
    for tool in df["tool"].unique():
        subset = df[df["tool"] == tool]
        plt.scatter(
            subset["wall_clock"], subset["max_mem"], label=tool, s=100
        )  # S = dimensione punti

    # Etichette e titolo
    plt.xlabel("wall_clock")
    plt.ylabel("max_mem")
    plt.title("index")
    plt.legend(title="Tool")
    plt.savefig(f"{out_pref}/index.pdf", dpi=500)


if __name__ == "__main__":
    main(sys.argv[1:])
