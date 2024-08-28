#!/usr/bin/env python3
import sys

##import gfapy
import random
import networkx as nx


def main(argv):
    gfa_file = argv[0]
    n_reads = int(argv[1])
    node_map = {}
    # gfa = gfapy.Gfa.from_file(gfa_file)
    G = nx.DiGraph()
    nn = []
    with open(gfa_file, "r") as f:
        inn = 0
        for line in f:
            if line.startswith("S"):
                tokens = line.split()
                G.add_node(
                    int(tokens[1]), ids=int(tokens[1]), seq=tokens[2], l=len(tokens[2])
                )
                nn.append(int(tokens[1]))
                node_map[int(tokens[1])] = inn
                inn = inn + 1
            elif line.startswith("L"):
                tokens = line.split()
                G.add_edge(int(tokens[1]), int(tokens[3]))
    # for s in gfa.segments:
    #     tokens = str(s).split()
    #     G.add_node(int(tokens[1]), ids=int(tokens[1]), seq=tokens[2], l=len(tokens[2]))

    # for e in gfa.edges:
    #     tokens = str(e).split()
    #     G.add_edge(int(tokens[1]), int(tokens[3]))
    # print(node_map)
    for i in range(n_reads):
        id_node = random.choice(nn)
        start_node = G.nodes()[id_node]
        start_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
        n_adds_nodes = random.choice(range(500))
        if n_adds_nodes < 50:
            i = i - 1
            continue
        if n_adds_nodes == 0:
            end_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
            s = start_node["seq"][start_offset : start_node["l"] - end_offset]
            n = start_node["ids"]
            m = node_map[n]
            if len(s) >= 50:
                if len(s) >= 10000:
                    print(f">{n}-{m}\n{s[:10000]}")
                else:
                    print(f">{n}-{m}\n{s}")
        else:
            list_nodes = [start_node["ids"]]
            end_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
            s = start_node["seq"][start_offset:]
            n = start_node["ids"]
            if len(s) >= 10000:
                s = s[len(s) - 1000 :]
            j = 1
            end = False
            while j < n_adds_nodes and not end:
                nb = list(G.neighbors(start_node["ids"]))
                if len(nb) == 0:
                    break
                id_node = random.choice(range(len(nb)))
                start_node = G.nodes()[nb[id_node]]
                list_nodes.append(start_node["ids"])
                if j != n_adds_nodes - 1:
                    if (len(s) + len(start_node["seq"])) <= 10000:
                        s = s + start_node["seq"]
                    else:
                        off = 500 - len(s)
                        s = s + start_node["seq"][:off]
                        end = True
                else:
                    end_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
                    if (
                        len(s) + len(start_node["seq"][: start_node["l"] - end_offset])
                    ) <= 10000:
                        s = s + start_node["seq"][: start_node["l"] - end_offset]
                    else:
                        off = 500 - len(s)
                        s = s + start_node["seq"][:off]
                        end = True
                j = j + 1
            ns = ">".join([f"{x}-{node_map[x]}" for x in list_nodes])
            if len(s) >= 256 and len(s) <= 10000:
                print(f">{ns}\n{s}")


if __name__ == "__main__":
    main(sys.argv[1:])
