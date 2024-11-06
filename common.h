#ifndef COMMON_H_
#define COMMON_H_
#include "gfa.h"
#include "kseq.h"
#include "mfmi/rlcsa.hpp"
#include "rlcsa.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "utils.h"

bool verbose = false;

int log4(int x) { return log(x) / log(4); }

struct node_sai {
  rldintv_t sai;
  uint64_t curr_node;
  std::vector<uint64_t> path;
  uint64_t end() const { return sai.x[0] + sai.x[2]; }
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
      merged.back().sai.x[2] = std::max(merged.back().end(), interval.end()) -
                               merged.back().sai.x[0];
    }
  }
  return merged;
}
uint8_t get_bwt_symb(const rld_t *index, uint64_t pos) {
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

uint64_t findnode(const rld_t *index, uint64_t pos, uint64_t end,
                  std::vector<std::vector<uint64_t>> &tags) {
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
  uint64_t end_node = tags[0][rld_rank11(index, sai.x[0], 0)];
  if (end_node == end - 1) {
    return 0;
  } else {
    return end_node + 1;
  }
}

std::vector<uint64_t> get_end_nodes(const rld_t *index, uint64_t b_i,
                                    uint64_t l_i,
                                    std::vector<std::vector<uint64_t>> &tags,
                                    std::vector<std::vector<uint64_t>> &adj) {
  std::vector<uint64_t> nodes;
  uint64_t d_b = rld_rank11(index, b_i, 0);
  uint64_t d_e = rld_rank11(index, b_i + l_i, 0);
  for (uint64_t j = d_b; j < d_e; j++) {
    uint64_t end_node = 0;
    if (tags[0][j] != tags[0].size() - 1)
      end_node = tags[0][j] + 1;
    nodes.push_back(end_node);
  }
  return nodes;
}

std::vector<std::pair<uint64_t, uint64_t>>
get_adj_nodes_f(const rld_t *index, uint64_t end_node,
                std::vector<std::vector<uint64_t>> &tags,
                std::vector<std::vector<uint64_t>> &adj) {
  std::vector<std::pair<uint64_t, uint64_t>> nodes;
  if (verbose) {
    std::fprintf(stderr, "adding siblings of %ld: ", end_node);
  }
  for (uint64_t k = 0; k < adj[end_node].size(); k++) {
    if (verbose) {
      std::fprintf(stderr, " sibling of %ld: %ld\n", end_node, adj[end_node][k]);
    }
    nodes.push_back(
        std::make_pair(tags[1][adj[end_node][k]], adj[end_node][k]));
  }
  if (verbose) {
    std::fprintf(stderr, "\n");
  }
  return nodes;
}

std::vector<node_sai> get_intervals(const rld_t *index, uint64_t b_i,
                                    uint64_t l_i,
                                    std::vector<std::vector<uint64_t>> &tags,
                                    std::vector<std::vector<uint64_t>> &adj,
                                    uint8_t symb, std::vector<uint64_t> p,
                                    uint64_t e_n = UINT64_MAX) {
  std::vector<node_sai> intervals;
  // if (e_n == INT_MAX) {
  std::vector<uint64_t> end_nodes = get_end_nodes(index, b_i, l_i, tags, adj);
  //if (verbose && end_nodes.size() >= 1) {
  //if(end_nodes.size()>=1){  
  // // std::cerr << "analyzing end nodes for interval [" << b_i << ", "
  // //           << b_i + l_i << "]: ";
  //  uint64_t d_b = rld_rank11(index, b_i, 0);
  //  uint64_t d_e = rld_rank11(index, b_i + l_i, 0);
  //  //std::cerr << "  (" << d_b << ", " << d_e << "): ";
  //  //for (auto n : end_nodes) {
  //  //  std::cerr << n << " ";
  //  //}
  //  //std::cerr << "\n";
  //}
  for (uint64_t n : end_nodes) {
    std::vector<std::pair<uint64_t, uint64_t>> adj_nodes = get_adj_nodes_f(index, n, tags, adj);
    //std::cerr << "sibling: ";
    for (auto an : adj_nodes) {
     if (verbose){
      std::cerr << "(" << an.first <<  ", " << an.second << ") ";
     }
    //std::cerr << "\n";
      rldintv_t sai_t;
      sai_t.x[0] = an.first;
      sai_t.x[1] = an.first;
      sai_t.x[2] = 1;
      if (verbose) {
        fprintf(stderr, "comparing %d vs %d\n", get_bwt_symb(index, an.first),
                symb);
      }
      /* std::vector<uint64_t> nn = {n}; */
      /* intervals.push_back({sai_t, an.second, nn}); */
      if (symb != 6) {
       //fprintf(stderr, "%d vs %d \n", get_bwt_symb(index, an.first), symb);
        if (get_bwt_symb(index, an.first) == symb) {
          //  std::vector<uint64_t> nn = {n, an.second};
          std::vector<uint64_t> nn = {n};
          intervals.push_back({sai_t, an.second, nn});
          // std::fprintf(stderr, "size %ld \n",
          //         intervals[intervals.size() - 1].path[0]);
        }
      } else {
        std::vector<uint64_t> nn = {n};
        intervals.push_back({sai_t, an.second, nn});
      }
    }
  }
  //fprintf(stderr, "with %d ", symb);
  //std::cerr << " interval size: " << intervals.size() << "\n";
  return intervals;
}

std::vector<node_sai> ext_int(const rld_t *index, const uint8_t *s, int i,
                              uint64_t l, rldintv_t &sai,
                              std::vector<std::vector<uint64_t>> &tags,
                              std::vector<std::vector<uint64_t>> &adj,
                              std::vector<uint64_t> labels_map,
                              uint64_t curr_node) {

  // uint64_t tollerance = 0;

  uint8_t symb = i > 0 ? s[i] : 5;
  // std::vector<uint64_t> match_nodes;
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
  // bool start = true;
  int db = 0;
  for (; i >= 0; --i) {
    db++;
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

    // if (i == l - 2) {
    //   for (auto ic2 : int_curr) {
    //     std::fprintf(stderr, "adding [%ld, %ld, %ld]\n", ic2.sai.x[0],
    //                  ic2.sai.x[1], ic2.sai.x[2]);
    //   }
    // }
    if (int_curr.size() == 0) {
      return {};
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
  return int_curr;
}

#endif // COMMON_H_
