import sys
import csv


def main(argv):
    fieldnames = [
        "tool",
        "task",
        "n_query",
        "l_query",
        "user_time",
        "sys_time",
        "max_mem",
        "wall_clock",
    ]
    res = {
        "tool": argv[0],
        "task": argv[1],
        "n_query": int(argv[2]),
        "l_query": int(argv[3]),
    }
    # print(sys.stdin)
    for line in sys.stdin:
        line = line[1:-1]
        tokens = line.split(sep=":")
        if tokens[0] == "User time (seconds)":
            res["user_time"] = float(tokens[1].lstrip())
        if tokens[0] == "System time (seconds)":
            res["sys_time"] = float(tokens[1].lstrip())
        if tokens[0] == "Maximum resident set size (kbytes)":
            res["max_mem"] = int(tokens[1].lstrip())
        if tokens[0] == "Elapsed (wall clock) time (h":
            tot = 0.0
            for x in tokens[4:]:
                tot = tot * 60 + float(x.lstrip())
                res["wall_clock"] = tot

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(res)


if __name__ == "__main__":
    main(sys.argv[1:])
