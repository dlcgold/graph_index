#!/usr/bin/env python3
import sys

##import gfapy
import random
import networkx as nx


def main(argv):
    gfa_file = argv[0]
    n_reads = int(argv[1])
    m_reads = int(argv[2])
    l_reads = int(argv[3])
    if m_reads > l_reads:
        m_reads = int(l_reads / 10)
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
    i = 0
    while i < n_reads:
        id_node = random.choice(nn)
        start_node = G.nodes()[id_node]
        start_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
        n_adds_nodes = random.choice(range(3 * int(l_reads / 32)))
        # if n_adds_nodes < 3:
        #     continue
        if n_adds_nodes == 0:
            end_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
            s = start_node["seq"][start_offset : start_node["l"] - end_offset]
            n = start_node["ids"]
            m = node_map[n]
            i = i + 1
            if len(s) >= l_reads:
                print(f">{n}-{m}\n{s[:l_reads]}")
            elif len(s) == l_reads:
                print(f">{n}-{m}\n{s}")
            else:
                i = i - 1
                # continue
        else:
            list_nodes = [start_node["ids"]]
            end_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
            s = start_node["seq"][start_offset:]
            n = start_node["ids"]
            if len(s) >= l_reads:
                s = s[len(s) - int(0.9 * len(s)) :]
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
                    if (len(s) + len(start_node["seq"])) <= l_reads:
                        s = s + start_node["seq"]
                    else:
                        off = l_reads - len(s) + 2
                        s = s + start_node["seq"][:off]
                        end = True
                else:
                    end_offset = random.choice(range(int((start_node["l"] + 1) / 2)))
                    if (
                        len(s) + len(start_node["seq"][: start_node["l"] - end_offset])
                    ) <= l_reads:
                        s = s + start_node["seq"][: start_node["l"] - end_offset]
                    else:
                        off = l_reads - len(s) + 2
                        s = s + start_node["seq"][:off]
                        end = True
                j = j + 1
            # ns = ">".join([f"{x}-{node_map[x]}" for x in list_nodes])
            ns = f"{list_nodes[0]}-{node_map[list_nodes[0]]}"
            ##print(len(s))
            i = i + 1
            if len(s) == l_reads:
                print(f">{ns}\n{s}")
            elif len(s) >= l_reads:
                print(f">{ns}\n{s[:l_reads]}")
            else:
                i = i - 1
                # continue


if __name__ == "__main__":
    main(sys.argv[1:])
