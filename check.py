#!/usr/bin/env python3
import sys


def main(argv):
    gfa_file = argv[0]
    match_file = argv[1]
    gaf_file = argv[2]

    node_map = {}

    my_matches = {}
    ga_matches = {}
    reads = set()
    with open(gfa_file, "r") as f:
        i = 0
        for line in f:
            if line.startswith("S"):
                tokens = line.split()
                node_map[str(i)] = tokens[1]
                i += 1

    with open(match_file, "r") as f:
        curr_token = "token"
        for line in f:
            if line.startswith("read"):
                tokens = line.split()
                curr_token = tokens[1]
                my_matches[curr_token] = []
                reads.add(curr_token)
            else:
                tokens = line.split()
                my_matches[curr_token].append(node_map[tokens[len(tokens) - 1]])

    # print(my_matches)

    with open(gaf_file, "r") as f:
        for i, line in enumerate(f):
            tokens = line.split()
            curr_token = tokens[0]
            if tokens[5][0] == "<":
                continue
            if "cg:Z:151=" not in line:
                continue
            node = tokens[5].split(">")[1]
            reads.add(curr_token)
            if curr_token not in ga_matches.keys():
                ga_matches[curr_token] = [node]
            else:
                ga_matches[curr_token].append(node)

    # print(ga_matches)
    for read in reads:
        print(read)
        if read not in my_matches.keys():
            print(read, "not in my")
            my_m = []
        else:
            my_m = my_matches[read]
            my_m.sort()
        if read not in ga_matches.keys():
            print(read, "not in ga")
            ga_m = []
        else:
            ga_m = ga_matches[read]
            ga_m.sort()
        print(my_m, ga_m)
        # my_m = my_matches[read].sort()
        # ga_m = ga_matches[read].sort()
        # if my_m != ga_m:
        #     print(my_m, ga_m)
        #     continue
        if my_m == ga_m:
            print("good")
        else:
            print("bad")


if __name__ == "__main__":
    main(sys.argv[1:])
