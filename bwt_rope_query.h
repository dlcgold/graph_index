//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_QUERY_H
#define GRAPHINDEX_BWT_ROPE_QUERY_H
#include "utils.h"

#include "kseq.h"
#include "rld0.h"
#include "utils.h"
#include <algorithm>
#include <cstdio>
#include <future>
#include <iostream>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/types.h>
#include <thread>
#include <utility>
#include <vector>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

bool verbose = false;
unsigned long long int interval_size = 0;
unsigned long long int interval_size_avg = 1;

std::vector<std::vector<unsigned int>> interval_spectrum;

int log4(int x) { return log(x) / log(4); }

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
std::vector<std::vector<int>> debug;

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

std::vector<node_sai> ext_int(const rld_t *index, const uint8_t *s, int i,
                              unsigned int l, rldintv_t &sai,
                              std::vector<std::vector<unsigned int>> &tags,
                              std::vector<std::vector<unsigned int>> &adj,
                              std::vector<unsigned int> labels_map,
                              unsigned int curr_node) {

  interval_spectrum.push_back({});
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

void ext(const rld_t *index, const uint8_t *s, int i, unsigned int l,
         rldintv_t &sai, std::vector<std::vector<unsigned int>> &tags,
         std::vector<std::vector<unsigned int>> &adj,
         std::vector<unsigned int> labels_map, unsigned int curr_node,
         char *read_name, unsigned int r_c, std::vector<node_sai> int_s) {

  interval_spectrum.push_back({});
  unsigned int tollerance = l / 5;
  tollerance = 0;

  uint8_t symb = i > 0 ? s[i] : 5;
  i--;
  // std::vector<unsigned int> match_nodes;
  // std::vector<node_sai> int_curr =
  //     get_intervals(index, sai.x[0], sai.x[2], tags, adj, symb, {});

  // int_curr.push_back({sai, curr_node});
  //
  std::vector<node_sai> int_curr = int_s;

  for (auto t : int_s) {
    auto tmp_int =
        get_intervals(index, t.sai.x[0], t.sai.x[2], tags, adj, symb, {});
    int_curr.insert(int_curr.end(), tmp_int.begin(), tmp_int.end());
  }
  int_curr = merge(int_curr);
  if (verbose) {
    std::fprintf(stderr, "at %d for char %d:\n", i, s[i]);
    for (auto ic : int_curr) {
      std::fprintf(stderr, "[%ld, %ld, %ld]\n", ic.sai.x[0], ic.sai.x[1],
                   ic.sai.x[2]);
    }
  }
  std::vector<node_sai> int_next;
  bool start = true;
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
    // interval_size += int_curr.size();
    interval_spectrum[r_c].push_back(int_curr.size());
    // if (db == 5) {
    //   for (auto in : int_curr) {
    //     debug.push_back(
    //         {in.sai.x[0], in.sai.x[1], in.sai.x[2], in.sai.info,
    //         in.curr_node});
    //   }
    //   int n[5] = {1, 1, 1, 1, 1};
    //   std::cout << debug.size() << " " << debug.size() * sizeof(n) <<
    //   std::endl; std::cout << debug.size() << " " << debug.size() *
    //   sizeof(debug.at(0))
    //             << std::endl;
    //   for (auto n : debug.at(0)) {
    //     std::cout << n << std::endl;
    //   }
    // }
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
          uint8_t symb = get_bwt_symb(index, k);
          if (symb == 0) {
            auto end_node = tags[0][rld_rank11(index, sai.x[0], 0)];
            if (end_node == tags[0].size() - 1) {
              std::fprintf(stdout, "@%s\t%d\n", read_name, labels_map[0]);
            } else {
              std::fprintf(stdout, "@%s\t%d\n", read_name,
                           labels_map[end_node + 1]);
            }
          } else {
            auto node_f = findnode(index, k, tags[0].size(), tags);
            // std::fprintf(stderr, "Match at node: %ld\n", node_f);
            std::fprintf(stdout, "@%s\t%d\n", read_name, labels_map[node_f]);
          }
        }
      }
    } else {
      std::fprintf(stdout, "@%s\tNO MATCH\n", read_name);
    }
  }
}

void query_bwt_rope(std::string index_pre, const char *query_file,
                    int threads) {
  auto i_file = index_pre + ".bwt";
  auto t_file = index_pre + ".dollars";
  auto g_file = index_pre + ".graph";
  auto l_file = index_pre + ".labels";
  auto c_file = index_pre + ".cache";

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

  std::vector<std::vector<std::vector<unsigned int>>> c_int =
      uintmat3_load(c_file.c_str());

  std::vector<std::vector<node_sai>> cache(c_int.size());

  std::vector<std::string> suff(cache.size());
  std::map<std::string, unsigned int> suff_map;

  unsigned int ls = log4(cache.size());
  std::cerr << "seeds of size " << cache.size() << " " << ls << "\n";
  for (int i = 0; i < cache.size(); ++i) {
    std::string suff_t(ls, ' ');
    int cur = i;

    for (int j = 0; j < ls; j++) {
      suff_t[j] = "ACGT"[cur % 4];
      cur /= 4;
    }
    // std::cout << suff_t << "\n";
    suff_map[suff_t] = i;
    // std::cout << suff[i] << "\n";
  }
  //
  int k = 0;
  for (auto s : c_int) {
    std::vector<node_sai> tv;
    rldintv_t sai;
    unsigned int curr_node;
    for (auto v : s) {
      sai.info = v[0];
      sai.x[0] = v[1];
      sai.x[1] = v[2];
      sai.x[2] = v[3];
      tv.push_back({sai, v[4], {v[4]}});
      // std::cout << sai.x[0] << " " << sai.x[1] << " " << sai.x[2] << " " <<
      // v[4]
      //           << "\n";
    }
    cache[k] = tv;
    k++;
  }

  // for (int i = 0; i < 1024; ++i) {
  //   std::cout << i << " " << c_int[i].size() << "\n";
  // }

  //  for (auto v : adj) {
  //    std::cout << cc << ": ";
  //    for (auto n : v) {
  //      std::cout << n << " ";
  //    }
  //    std::cout << "\n";
  //    cc++;
  //  }
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

  std::fflush(stdout);
  while ((l = kseq_read(ks)) >= 0) {
    std::fprintf(stderr, "%d\r", r_c);

    s = (uint8_t *)ks->seq.s;

    std::string suf(ls, ' ');
    int sufc = 0;
    // change encoding
    for (i = 0; i < l; ++i) {
      s[i] = fm6_i(s[i]);
      if (i >= l - ls) {
        suf[sufc] = map[s[i]];
        sufc++;
      }
    }

    // std::cout << suf << "\n";
    auto int_v = cache[suff_map[suf]];
    if (int_v.empty())
      continue;
    // for (auto i : int_v) {
    //   std::cout << i.sai.x[0] << " " << i.sai.x[1] << " " << i.sai.x[2] << "
    //   "
    //             << i.curr_node << "\n";
    // }
    // for (i = 0; i < l; ++i) {
    //   printf("%d ", s[i]);
    // }
    // std::cout << "\n";
    // rldintv_t osai[6];
    // l -= (ls-1);

    l -= (ls - 1);
    i = l - 1;
    // i = l;
    // fm6_set_intv(index, s[i], sai);
    // --i;
    // auto s_b = sai.x[1];
    // sai.x[1] = tags[0].size();

    // auto nodes = ext(index, s, i, l, sai, tags, adj, tags[0].size());
    ext(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(), ks->name.s,
        r_c, int_v);

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
    //   // auto nodes = ext(index, s, i, l, sai, tags, adj, tags[0].size());
    //   ext(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(),
    //       ks->name.s, r_c);
    //   // std::fprintf(stderr, "finished");
    //   sai.x[1] = s_b;
    // }
    r_c++;
    // std::fprintf(stderr, "\n-----------------\n");
  }
  // std::fprintf(stderr, "~Avg total intervals considered for %d reads:
  // %lld\n",
  //              r_c, (unsigned long long int)interval_size / r_c);
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);
  // for (auto r : interval_spectrum) {
  //   for (auto rr : r) {
  //     std::cerr << rr << " ";
  //   }
  //   std::cerr << std::endl;
  // }
}

#endif // GRAPHINDEX_BWT_ROPE_QUERY_H
