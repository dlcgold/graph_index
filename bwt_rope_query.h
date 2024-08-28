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

struct node_match {
  unsigned int node = 0;
  unsigned int l_off = 0;
  unsigned int r_off = 0;

  bool operator<(const node_match &other) const {
    if (node == other.node) {
      if (l_off == other.l_off) {
        return r_off < other.r_off;
      } else {
        return l_off < other.l_off;
      }
    } else {
      return node < other.node;
    }
  }
};

std::ostream &operator<<(std::ostream &os, node_match const &n) {
  os << "[" << n.node << ", " << n.l_off << ", " << n.r_off << "]";
  return os;
}

struct node_sai {
  rldintv_t sai;
  unsigned int curr_node;
  std::vector<unsigned int> path;

  int end() const { return sai.x[0] + sai.x[2]; }
  bool operator<(const node_sai &other) const {
    return sai.x[0] < other.sai.x[0];
  }
  bool operator==(const node_sai &other) const {
    return sai.x[0] == other.sai.x[0] && sai.x[2] == other.sai.x[2];
  }
};

std::vector<node_sai> merge(std::vector<node_sai> intervals) {

  std::sort(intervals.begin(), intervals.end());

  intervals.erase(std::unique(intervals.begin(), intervals.end()),
                  intervals.end());
  std::vector<node_sai> merged;
  for (const auto &interval : intervals) {
    if (merged.empty() || merged.back().end() < interval.sai.x[0]) {
      merged.push_back(interval);
    } else {
      merged.back().sai.x[2] =
          max(merged.back().end(), interval.end()) - merged.back().sai.x[0];
    }
  }
  return merged;
}
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

node_match findnode_off(const rld_t *index, unsigned int pos, unsigned end,
                        std::vector<std::vector<unsigned int>> &tags,
                        unsigned int l) {
  // std::fprintf(stderr, "pos %d", pos);
  node_match r;
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
    r.l_off++;
  }
  // std::cout << "dollar: " << rld_rank11(index, sai.x[0], 0) << " at "
  //           << sai.x[0] << "\n";
  auto end_node = tags[0][rld_rank11(index, sai.x[0], 0)];
  if (end_node != end - 1) {
    r.node = end_node + 1;
  }
  r.r_off = tags[2][r.node] - r.l_off - l;
  // sai.x[0] = pos;
  // sai.x[1] = pos;
  // sai.x[2] = 1;
  // rldintv_t osair[6];
  // while (symb != 0) {
  //   rld_extend(index, &sai, osair, 0);
  //   sai = osair[fm6_comp(symb)];
  //   symb = get_bwt_symb(index, sai.x[0]);
  //   r.r_off++;
  // }
  return r;
}

unsigned int findnode(const rld_t *index, unsigned int pos, unsigned end,
                      std::vector<std::vector<unsigned int>> &tags) {
  // std::fprintf(stderr, "pos %d", pos);
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
    std::fprintf(stderr, "adding siblings of %d: ", end_node);
  }
  for (unsigned int k = 0; k < adj[end_node].size(); k++) {
    if (verbose) {
      std::fprintf(stderr, " sibling of %d: %d\n", end_node, adj[end_node][k]);
    }
    nodes.push_back(
        std::make_pair(tags[1][adj[end_node][k]], adj[end_node][k]));
  }
  if (verbose) {
    std::fprintf(stderr, "\n");
  }
  return nodes;
}

std::vector<node_sai>
get_intervals(const rld_t *index, unsigned int b_i, unsigned int l_i,
              std::vector<std::vector<unsigned int>> &tags,
              std::vector<std::vector<unsigned int>> &adj, uint8_t symb,
              std::vector<unsigned int> p, unsigned int e_n = INT_MAX) {
  std::vector<node_sai> intervals;
  // if (e_n == INT_MAX) {
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
    auto adj_nodes = get_adj_nodes_f(index, n, tags, adj);
    for (auto an : adj_nodes) {
      rldintv_t sai_t;
      sai_t.x[0] = an.first;
      sai_t.x[1] = an.first;
      sai_t.x[2] = 1;
      if (verbose) {
        fprintf(stderr, "comparing %d vs %d\n", get_bwt_symb(index, an.first),
                symb);
      }
      if (get_bwt_symb(index, an.first) == symb) {
        //  std::vector<unsigned int> nn = {n, an.second};
        std::vector<unsigned int> nn = {n};
        intervals.push_back({sai_t, an.second, nn});
        // std::fprintf(stderr, "size %ld \n",
        //         intervals[intervals.size() - 1].path[0]);
      }
    }
  }
  // } else {
  //   auto adj_nodes = get_adj_nodes_f(index, e_n, tags, adj);
  //   for (auto an : adj_nodes) {
  //     rldintv_t sai_t;
  //     sai_t.x[0] = an.first;
  //     sai_t.x[1] = an.first;
  //     sai_t.x[2] = 1;
  //     if (get_bwt_symb(index, an.first) == symb) {
  //       auto n_p = p;
  //       // n_p.push_back(an.second);
  //       intervals.push_back({sai_t, an.second, nn});
  //     }
  //   }
  // }

  // std::sort(intervals.begin(), intervals.end(),
  //           [](const node_sai &a, const node_sai &b) {
  //             return a.sai.x[0] < b.sai.x[0];
  //           });
  // std::vector<node_sai> merged;
  // for (const auto &interval : intervals) {
  //   if (merged.empty() || merged.back().end() < interval.sai.x[0]) {
  //     merged.push_back(interval);
  //   } else {
  //     merged.back().sai.x[2] =
  //         max(merged.back().end(), interval.end()) - merged.back().sai.x[0];
  //   }
  // }
  return intervals;
}

void ext(const rld_t *index, const uint8_t *s, int i, unsigned int l,
         rldintv_t &sai, std::vector<std::vector<unsigned int>> &tags,
         std::vector<std::vector<unsigned int>> &adj,
         std::vector<unsigned int> labels_map, unsigned int curr_node,
         char *read_name) {

  unsigned int tollerance = l / 5;
  tollerance = 0;

  uint8_t symb = i > 0 ? s[i] : 5;
  // std::vector<unsigned int> match_nodes;
  std::vector<node_sai> int_curr =
      get_intervals(index, sai.x[0], sai.x[2], tags, adj, symb, {});

  int_curr.push_back({sai, curr_node});
  if (verbose) {
    std::fprintf(stderr, "at %d for char %d:\n", i, s[i]);
    for (auto ic : int_curr) {
      std::fprintf(stderr, "[%ld, %ld, %ld]\n", ic.sai.x[0], ic.sai.x[1],
                   ic.sai.x[2]);
    }
  }
  std::vector<node_sai> int_next;
  bool start = true;
  for (; i >= 0; --i) {
    if (verbose) {
      std::fprintf(stderr, "\n");
    }
    for (auto ic : int_curr) {
      if (verbose) {
        std::fprintf(stderr,
                     "analyzing [%ld, %ld, %ld] with char %d at pos %d\n",
                     ic.sai.x[0], ic.sai.x[1], ic.sai.x[2], s[i], i);
        std::fprintf(stderr, "in %ld we have %d\n", ic.sai.x[0],
                     get_bwt_symb(index, ic.sai.x[0]));
      }
      rldintv_t osai[6];
      rld_extend(index, &ic.sai, osai, 1);
      // for (int c = 0; c < 6; ++c) {
      //   osai[c].x[1] = ic.x[1];
      // }
      auto sai = osai[s[i]];
      if (verbose) {
        std::fprintf(stderr, "extended into [%ld, %ld, %ld]\n", sai.x[0],
                     sai.x[1], sai.x[2]);
        for (int c = 0; c < 6; ++c) {
          std::fprintf(stderr, "other with %d: [%ld, %ld, %ld]\n", c,
                       osai[c].x[0], osai[c].x[1], osai[c].x[2]);
        }
      }
      if (sai.x[2] > 0) {
        symb = i > 0 ? s[i - 1] : 5;

        auto tmp_int =
            get_intervals(index, sai.x[0], sai.x[2], tags, adj, symb, {});
        if (!tmp_int.empty()) {
          // for (auto nnn : tmp_int[0].path) {
          //   std::fprintf(stderr, "path with %ld \n", nnn);
          // }
          if (verbose) {

            for (auto ic2 : tmp_int) {
              std::fprintf(stderr, "adding [%ld, %ld, %ld]\n", ic2.sai.x[0],
                           ic2.sai.x[1], ic2.sai.x[2]);
            }
          }
          int_next.insert(int_next.end(), tmp_int.begin(), tmp_int.end());
        }

        int_next.push_back({sai, ic.curr_node, ic.path});

        // std::fprintf(stderr, "INTERVALS: %ld %ld\n", interval_size,
        // int_next.size());

        if (verbose) {
          for (auto ic : int_next) {
            std::fprintf(stderr, "tmp next [%ld, %ld, %ld]\n", ic.sai.x[0],
                         ic.sai.x[1], ic.sai.x[2]);
          }
          std::cout << "---------------\n";
        }
      }
    }
    // int_curr = int_next;
    int_curr = merge(int_next);
    interval_size += int_curr.size();
    // if (i == l - 2) {
    //   for (auto ic2 : int_curr) {
    //     std::fprintf(stderr, "adding [%ld, %ld, %ld]\n", ic2.sai.x[0],
    //                  ic2.sai.x[1], ic2.sai.x[2]);
    //   }
    // }
    if (int_curr.size() == 0) {
      return;
    }
    int_next.clear();
    if (verbose) {
      std::fprintf(stderr, "at %d for char %d finals intervals are:\n", i,
                   s[i]);
      for (auto ic2 : int_curr) {
        std::fprintf(stderr, "[%ld, %ld, %ld]\n", ic2.sai.x[0], ic2.sai.x[1],
                     ic2.sai.x[2]);
      }
    }
  }
  for (auto ic : int_curr) {
    if (ic.sai.x[2] >= 1) {
      // std::fprintf(stderr, "Final match [%ld, %ld, %ld]\n", ic.sai.x[0],
      // ic.sai.x[1],
      //         ic.sai.x[2]);

      if (ic.curr_node != tags[0].size()) {
        if (ic.path.size() > 1) {
          // std::fprintf(stderr, "Match at nodes: ");
          std::string matches = "";
          for (int i = ic.path.size() - 1; i >= 0; i--) {
            // std::fprintf(stderr, "%ld ", ic.path[i]);
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
          std::cout << "@" << read_name << "\t" << matches << "\n";
          // std::fprintf(stderr, "%s\t%s\n", read_name, matches);
          //  std::fprintf(stderr, "\n");
        } else {
          // std::fprintf(stderr, "Match at node: %ld\n", ic.curr_node);
          std::fprintf(stdout, "@%s\t%d\n", read_name,
                       labels_map[ic.curr_node]);
        }
      } else {
        for (unsigned int k = ic.sai.x[0]; k < ic.sai.x[0] + ic.sai.x[2]; k++) {
          auto node_f = findnode(index, k, tags[0].size(), tags);
          // std::fprintf(stderr, "Match at node: %ld\n", node_f);
          std::fprintf(stdout, "@%s\t%d\n", read_name, labels_map[node_f]);
        }
      }
    } else {
      std::fprintf(stdout, "@%s\tNO MATCH\n", read_name);
    }
  }
}

std::pair<std::vector<uint8_t *>, unsigned int>
split_s(const uint8_t *s, unsigned int l, unsigned int s_p) {
  std::vector<uint8_t *> r;

  size_t start = l;
  size_t f_l = 0;

  while (start > 0) {
    size_t end = start;
    start = (start > s_p) ? start - s_p : 0;

    uint8_t *segment = new uint8_t[end - start];
    std::copy(s + start, s + end, segment);

    if (start == 0) {
      f_l = end - start;
    }
    r.insert(r.begin(), segment);
  }
  return std::make_pair(r, f_l);
}

std::vector<node_match>
locate_nodes(const rld_t *index, const uint8_t *s, unsigned int l,
             std::vector<std::vector<unsigned int>> &tags) {
  std::vector<node_match> nodes;
  // for (int i = 0; i < l; i++) {
  //   std::cout << "$ACTG"[s[i]];
  // }
  // std::cout << "\n";
  rldintv_t sai;

  int i = l - 1;
  rldintv_t osai[6];
  fm6_set_intv(index, s[i], sai);
  --i;

  for (; i >= 0; --i) {
    rld_extend(index, &sai, osai, 1);
    sai = osai[s[i]];
    if (sai.x[2] <= 0) {
      break;
    }
  }

  for (unsigned int k = sai.x[0]; k < sai.x[0] + sai.x[2]; k++) {
    auto node_f = findnode_off(index, k, tags[0].size(), tags, l);
    // std::fprintf(stderr, "Match at node: %ld\n", node_f);
    nodes.push_back(node_f);
  }
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

void match(const rld_t *index, const uint8_t *s, unsigned int l,
           unsigned int s_p, std::vector<std::vector<unsigned int>> &tags) {
  auto sub_s = split_s(s, l, s_p);
  for (int i = 0; i < l; i++) {
    std::cout << "$ACGT"[s[i]];
  }
  std::cout << "\n";

  for (int i = 0; i < sub_s.first.size(); i++) {
    auto ll = (i == 0) ? sub_s.second : s_p;
    for (int j = 0; j < ll; j++) {
      std::cout << "$ACGT"[sub_s.first[i][j]];
    }
    std::cout << "\n";
  }

  // std::cout << "\n";
  auto f_n = std::vector<std::vector<node_match>>(sub_s.first.size());
  for (int i = 0; i < sub_s.first.size(); i++) {

    auto ll = (i == 0) ? sub_s.second : s_p;
    // std::cout << ll << "\n";
    // for (int j = 0; j < ll; j++) {
    //   std::cout << "$ACTG"[sub_s.first[i][j]];
    // }
    // std::cout << "\n";
    // fflush(stdout);
    f_n[i] = locate_nodes(index, sub_s.first[i], ll, tags);
    // for (auto n : f_n[i]) {
    //   std::cout << n;
    // }
    // std::cout << "\n";
  }
  int i = 0;
  std::cout << f_n.size() << "\n";
  for (auto nn : f_n) {
    std::cout << i << ": ";
    for (auto n : nn) {
      std::cout << n << " ";
    }
    std::cout << "\n";
    i++;
  }
}

void query_bwt_rope(std::string index_pre, const char *query_file,
                    int threads) {
  auto i_file = index_pre + ".bwt";
  auto t_file = index_pre + ".dollars";
  auto g_file = index_pre + ".graphi";
  auto g_file2 = index_pre + ".grapho";
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
  std::vector<std::vector<unsigned int>> adj2 = uintmat_load(g_file2.c_str());
  unsigned int max_l = adj[adj.size() - 1][0];
  unsigned int avg_l = adj[adj.size() - 1][1];
  // unsigned int s_p = 0;
  // if (max_l >= 32) {
  //   s_p = 32;
  // } else {
  //   s_p = max_l;
  // }
  // s_p = avg_l;

  auto s_p = threads;
  if (s_p > avg_l) {
    s_p = avg_l;
  }
  std::cerr << "max_l = " << max_l << " avg_l = " << avg_l << " s_p = " << s_p
            << std::endl;
  std::vector<unsigned int> labels_map = uintvec_load(l_file.c_str());
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

  gzFile fp = gzopen(query_file, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  int r_c = 0;
  while ((l = kseq_read(ks)) >= 0) {
    std::fprintf(stderr, "%d\r", r_c);
    s = (uint8_t *)ks->seq.s;
    for (i = 0; i < l; ++i) {
      printf("%c", ks->seq.s[i]);
    }
    printf("\n");

    // change encoding
    for (i = 0; i < l; ++i) {
      s[i] = fm6_i(s[i]);
    }
    for (i = 0; i < l; ++i) {
      printf("%d", s[i]);
    }

    printf("\n");
    // i = l - 1;
    // fm6_set_intv(index, s[i], sai);
    // if (i == 0) {

    //   for (unsigned int k = sai.x[0]; k < sai.x[0] + sai.x[2]; k++) {
    //     auto node_f = findnode(index, k, tags[0].size(), tags);
    //     // std::fprintf(stderr, "Match at node %d\n", node_f);
    //     std::fprintf(stderr, "%s\t%d\n", ks->name.s, labels_map[node_f]);
    //   }

    // } else {
    //   --i;
    //   auto s_b = sai.x[1];
    //   sai.x[1] = tags[0].size();
    //   ext(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(),
    //       ks->name.s);
    //   sai.x[1] = s_b;
    // }
    // r_c++;

    match(index, s, l, s_p, tags);
  }
  // std::fprintf(stderr, "~Avg total intervals considered for %d reads:
  // %lld\n",
  //              r_c, (unsigned long long int)interval_size / r_c);
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);
}

#endif // GRAPHINDEX_BWT_ROPE_H
