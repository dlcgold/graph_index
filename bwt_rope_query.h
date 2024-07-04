//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_H
#define GRAPHINDEX_BWT_ROPE_H

#include <algorithm>
#include <cstdio>
#include <future>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/types.h>
#include <thread>
#include <utility>
#include <vector>
#include <zlib.h>

#include "kseq.h"
#include "rld0.h"
#include "utils.h"

#include <cstdint>
KSEQ_INIT(gzFile, gzread)

bool verbose = false;
unsigned long long int interval_size = 0;
unsigned long long int interval_size_avg = 1;

struct node_sai {
  rldintv_t sai;
  unsigned int curr_node;
  std::vector<unsigned int> path;
};

uint8_t get_bwt_symb(const rld_t *index, unsigned int pos) {
  rldintv_t sai;
  uint8_t symb = 5;
  sai.x[0] = pos;
  sai.x[1] = pos;
  sai.x[2] = 1;
  rldintv_t osai[6];
  rld_extend(index, &sai, osai, 1);
  for (uint8_t s = 0; s < 6; s++) {
    if (osai[s].x[2] == 1) {
      symb = s;
      break;
    }
  }
  return symb;
}
unsigned int findnode(const rld_t *index, unsigned int pos, unsigned end,
                      std::vector<std::vector<unsigned int>> &tags) {
  // fprintf(stderr, "pos %d", pos);
  uint8_t symb = get_bwt_symb(index, pos);
  rldintv_t sai;
  sai.x[0] = pos;
  sai.x[1] = pos;
  sai.x[2] = 1;
  rldintv_t osai[6];
  while (symb != 0) {
    rld_extend(index, &sai, osai, 1);
    sai = osai[symb];
    symb = get_bwt_symb(index, sai.x[0]);
  }
  // std::cout << "dollar: " << rld_rank11(index, sai.x[0], 0) << " at "
  //           << sai.x[0] << "\n";
  auto end_node = tags[0][rld_rank11(index, sai.x[0], 0)];
  if (end_node == end - 1) {
    return 0;
  } else {
    return end_node + 1;
  }
}

std::vector<unsigned int>
get_end_nodes(const rld_t *index, unsigned int b_i, unsigned int l_i,
              std::vector<std::vector<unsigned int>> &tags,
              std::vector<std::vector<unsigned int>> &adj) {
  std::vector<unsigned int> nodes;
  auto d_b = rld_rank11(index, b_i, 0);
  auto d_e = rld_rank11(index, b_i + l_i, 0);
  for (unsigned int j = d_b; j < d_e; j++) {
    auto end_node = 0;
    if (tags[0][j] != tags[0].size() - 1)
      end_node = tags[0][j] + 1;

    nodes.push_back(end_node);
  }
  return nodes;
}

std::vector<std::pair<unsigned int, unsigned int>>
get_adj_nodes_f(const rld_t *index, unsigned int end_node,
                std::vector<std::vector<unsigned int>> &tags,
                std::vector<std::vector<unsigned int>> &adj) {
  std::vector<std::pair<unsigned int, unsigned int>> nodes;
  if (verbose) {
    fprintf(stderr, "adding siblings of %d: ", end_node);
  }
  for (unsigned int k = 0; k < adj[end_node].size(); k++) {
    if (verbose) {
      fprintf(stderr, " sibling of %d: %d\n", end_node, adj[end_node][k]);
    }
    nodes.push_back(
        std::make_pair(tags[1][adj[end_node][k]], adj[end_node][k]));
  }
  if (verbose) {
    fprintf(stderr, "\n");
  }
  return nodes;
}

std::vector<node_sai>
get_intervals(const rld_t *index, unsigned int b_i, unsigned int l_i,
              std::vector<std::vector<unsigned int>> &tags,
              std::vector<std::vector<unsigned int>> &adj, uint8_t symb,
              std::vector<unsigned int> p, unsigned int e_n = INT_MAX) {
  std::vector<node_sai> intervals;
  if (e_n == INT_MAX) {
    auto end_nodes = get_end_nodes(index, b_i, l_i, tags, adj);
    if (verbose && end_nodes.size() >= 1) {
      std::cerr << "analyzing end nodes for interval [" << b_i << ", "
                << b_i + l_i << "]: ";
      auto d_b = rld_rank11(index, b_i, 0);
      auto d_e = rld_rank11(index, b_i + l_i, 0);
      std::cerr << "  (" << d_b << ", " << d_e << "): ";
      for (auto n : end_nodes) {
        std::cerr << n << " ";
      }
      std::cerr << "\n";
    }
    for (auto n : end_nodes) {
      // fprintf(stderr, "considering node %ld\n", n);
      auto adj_nodes = get_adj_nodes_f(index, n, tags, adj);
      for (auto an : adj_nodes) {
        rldintv_t sai_t;
        sai_t.x[0] = an.first;
        sai_t.x[1] = an.first;
        sai_t.x[2] = 1;
        if (get_bwt_symb(index, an.first) == symb) {
          std::vector<unsigned int> nn = {n, an.second};
          intervals.push_back({sai_t, an.second, nn});
          // fprintf(stderr, "size %ld \n",
          //         intervals[intervals.size() - 1].path[0]);
        }
      }
    }
  } else {
    auto adj_nodes = get_adj_nodes_f(index, e_n, tags, adj);
    for (auto an : adj_nodes) {
      rldintv_t sai_t;
      sai_t.x[0] = an.first;
      sai_t.x[1] = an.first;
      sai_t.x[2] = 1;
      if (get_bwt_symb(index, an.first) == symb) {
        auto n_p = p;
        n_p.push_back(an.second);
        intervals.push_back({sai_t, an.second, n_p});
      }
    }
  }
  return intervals;
}

void ext(const rld_t *index, const uint8_t *s, int i, unsigned int l,
         rldintv_t &sai, std::vector<std::vector<unsigned int>> &tags,
         std::vector<std::vector<unsigned int>> &adj,
         std::vector<unsigned int> labels_map, unsigned int curr_node,
         char *read_name) {

  unsigned int tollerance = l / 5;
  tollerance = 0;

  uint8_t symb = i > 0 ? s[i - 1] : 5;
  // std::vector<unsigned int> match_nodes;
  std::vector<node_sai> int_curr =
      get_intervals(index, sai.x[0], sai.x[2], tags, adj, symb, {});

  int_curr.push_back({sai, curr_node});
  if (verbose) {
    fprintf(stderr, "at %d for char %d:\n", i, s[i]);
    for (auto ic : int_curr) {
      fprintf(stderr, "[%ld, %ld, %ld]\n", ic.sai.x[0], ic.sai.x[1],
              ic.sai.x[2]);
    }
  }
  std::vector<node_sai> int_next;

  for (; i >= 0; --i) {
    if (verbose) {
      fprintf(stderr, "\n");
    }
    for (auto ic : int_curr) {
      if (verbose) {
        fprintf(stderr, "analyzing [%ld, %ld, %ld] with char %d at pos %d\n",
                ic.sai.x[0], ic.sai.x[1], ic.sai.x[2], s[i], i);
        fprintf(stderr, "in %d we have %d\n", ic.sai.x[0],
                get_bwt_symb(index, ic.sai.x[0]));
      }
      rldintv_t osai[6];
      rld_extend(index, &ic.sai, osai, 1);
      // for (int c = 0; c < 6; ++c) {
      //   osai[c].x[1] = ic.x[1];
      // }
      auto sai = osai[s[i]];
      if (verbose) {
        fprintf(stderr, "extended into [%ld, %ld, %ld]\n", sai.x[0], sai.x[1],
                sai.x[2]);
        for (int c = 0; c < 6; ++c) {
          fprintf(stderr, "other with %d: [%ld, %ld, %ld]\n", c, osai[c].x[0],
                  osai[c].x[1], osai[c].x[2]);
        }
      }
      if (sai.x[2] > 0) {
        if (i != 0 && i <= l - tollerance) {
          symb = i > 0 ? s[i - 1] : 5;
          if (sai.x[2] == 1 && ic.curr_node != tags[0].size() &&
              get_bwt_symb(index, sai.x[0]) == 0) {
            auto tmp_int = get_intervals(index, sai.x[0], sai.x[2], tags, adj,
                                         symb, ic.path, ic.curr_node);
            if (!tmp_int.empty()) {
              if (verbose) {
                for (auto ic2 : tmp_int) {
                  fprintf(stderr, "adding [%ld, %ld, %ld]\n", ic2.sai.x[0],
                          ic2.sai.x[1], ic2.sai.x[2]);
                }
              }
              int_next.insert(int_next.end(), tmp_int.begin(), tmp_int.end());
            }
          } else {
            auto tmp_int =
                get_intervals(index, sai.x[0], sai.x[2], tags, adj, symb, {});
            if (!tmp_int.empty()) {
              // for (auto nnn : tmp_int[0].path) {
              //   fprintf(stderr, "path with %ld \n", nnn);
              // }
              if (verbose) {
                for (auto ic2 : tmp_int) {
                  fprintf(stderr, "adding [%ld, %ld, %ld]\n", ic2.sai.x[0],
                          ic2.sai.x[1], ic2.sai.x[2]);
                }
              }
              int_next.insert(int_next.end(), tmp_int.begin(), tmp_int.end());
            }
          }
        }
        int_next.push_back({sai, ic.curr_node, ic.path});

        interval_size += int_next.size();
        // fprintf(stderr, "INTERVALS: %ld %ld\n", interval_size,
        // int_next.size());

        if (verbose) {
          for (auto ic : int_next) {
            fprintf(stderr, "tmp next [%ld, %ld, %ld]\n", ic.sai.x[0],
                    ic.sai.x[1], ic.sai.x[2]);
          }
        }
      }
    }
    int_curr = int_next;
    if (int_curr.size() == 0) {
      return;
    }
    int_next.clear();
    if (verbose) {
      fprintf(stderr, "at %d for char %d finals intervals are:\n", i, s[i]);
      for (auto ic2 : int_curr) {
        fprintf(stderr, "[%ld, %ld, %ld]\n", ic2.sai.x[0], ic2.sai.x[1],
                ic2.sai.x[2]);
      }
    }
  }
  for (auto ic : int_curr) {
    if (ic.sai.x[2] >= 1) {
      // fprintf(stderr, "Final match [%ld, %ld, %ld]\n", ic.sai.x[0],
      // ic.sai.x[1],
      //         ic.sai.x[2]);

      if (ic.curr_node != tags[0].size()) {
        if (ic.path.size() > 1) {
          // fprintf(stderr, "Match at nodes: ");
          std::string matches = "";
          for (int i = ic.path.size() - 1; i >= 0; i--) {
            // fprintf(stderr, "%ld ", ic.path[i]);
            if (i != 0) {
              std::ostringstream s;
              s << labels_map[ic.path[i]] << ">";
              matches = matches + s.str();
            } else {
              std::ostringstream s;
              s << labels_map[ic.path[i]];
              matches = matches + s.str();
            }
          }
          std::cerr << "@" << read_name << "\t" << matches << "\n";
          // fprintf(stderr, "%s\t%s\n", read_name, matches);
          //  fprintf(stderr, "\n");
        } else {
          // fprintf(stderr, "Match at node: %ld\n", ic.curr_node);
          fprintf(stderr, "@%s\t%d\n", read_name, labels_map[ic.curr_node]);
        }
      } else {
        for (unsigned int k = ic.sai.x[0]; k < ic.sai.x[0] + ic.sai.x[2]; k++) {
          auto node_f = findnode(index, k, tags[0].size(), tags);
          // fprintf(stderr, "Match at node: %ld\n", node_f);
          fprintf(stderr, "@%s\t%d\n", read_name, labels_map[node_f]);
        }
      }
    } else {
      fprintf(stderr, "@%s\tNO MATCH\n", read_name);
    }
  }
}

void query_bwt_rope(std::string index_pre, const char *query_file,
                    int threads) {
  auto i_file = index_pre + ".bwt";
  auto t_file = index_pre + ".dollars";
  auto g_file = index_pre + ".graph";
  auto l_file = index_pre + ".labels";

  rld_t *index = rld_restore(i_file.c_str());

  std::vector<std::vector<unsigned int>> tags = uintmat_load(t_file.c_str());

  // for (auto n : tags[0]) {
  //   std::cerr << n << " ";
  // }
  // std::cout << "\n-----------\n";
  // for (auto n : tags[1]) {
  //   std::cerr << n << " ";
  // }
  // std::cerr << "\n";

  std::vector<std::vector<unsigned int>> adj = uintmat_load(g_file.c_str());
  std::vector<unsigned int> labels_map = uintvec_load(l_file.c_str());
  int cc = 0;
  // for (auto v : adj) {
  //   std::cout << cc << ": ";
  //   for (auto n : v) {
  //     std::cout << n << " ";
  //   }
  //   std::cout << "\n";
  //   cc++;
  // }
  uint8_t *s;
  // q interval
  rldintv_t sai;
  for (int c = 0; c < 6; ++c) {
    fm6_set_intv(index, c, sai);
    // printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  }

  int errors = 0;
  gzFile fp = gzopen(query_file, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  int r_c = 0;
  while ((l = kseq_read(ks)) >= 0) {
    // fprintf(stderr, "read %s\n", ks->name.s);

    r_c++;
    s = (uint8_t *)ks->seq.s;
    // change encoding
    for (i = 0; i < l; ++i) {
      s[i] = fm6_i(s[i]);
    }
    // for (i = 0; i < l; ++i) {
    //   printf("%d ", s[i]);
    // }
    // std::cout << "\n";
    rldintv_t osai[6];
    i = l - 1;
    fm6_set_intv(index, s[i], sai);
    if (i == 0) {

      for (unsigned int k = sai.x[0]; k < sai.x[0] + sai.x[2]; k++) {
        auto node_f = findnode(index, k, tags[0].size(), tags);
        // fprintf(stderr, "Match at node %d\n", node_f);
        fprintf(stderr, "%s\t%d\n", ks->name.s, labels_map[node_f]);
      }

    } else {
      --i;
      auto s_b = sai.x[1];
      sai.x[1] = tags[0].size();
      // auto nodes = ext(index, s, i, l, sai, tags, adj, tags[0].size());
      ext(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(),
          ks->name.s);
      // fprintf(stderr, "finished");
      sai.x[1] = s_b;
    }
    r_c++;
    // fprintf(stderr, "\n-----------------\n");
  }
  fprintf(stderr, "~Avg total intervals considered for %lld reads: %lld\n", r_c,
          (unsigned long long int)interval_size / r_c);
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);
}

#endif // GRAPHINDEX_BWT_ROPE_H
