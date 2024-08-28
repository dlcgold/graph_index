#!/usr/bin/env python3
import sys


def main(argv):
    match_file = argv[0]

    my_m = {}
    m = {}
    with open(match_file, "r") as f:
        c = 0
        for line in f:
            line = line[1:]
            tokens = line.split()
            matches = tokens[0].split(">")
            for i in range(len(matches)):
                t = matches[i].split("-")
                matches[i] = t[0]
            m[tokens[0]] = matches
            my_matches = tokens[1].split(">")
            if tokens[0] in my_m.keys():
                my_m[tokens[0]].append(my_matches)
            else:
                my_m[tokens[0]] = my_matches
            # if matches != my_matches:
            #     print(f"EM at line {c} could mismatch")
            #     print(matches)
            #     print(my_matches)
            # c = c + 1

    for read, mm in my_m.items():
        found = False
        for mmm in mm:
            if mmm in m[read]:
                found = True
        if not found:
            print(f"read {read} mismatch")


if __name__ == "__main__":
    main(sys.argv[1:])
