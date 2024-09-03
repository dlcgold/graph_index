//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_H
#define GRAPHINDEX_BWT_ROPE_H

// #include "bwt_rope_query.h"

#include "gfa.h"
#include "kseq.h"
#include "mfmi/rlcsa.hpp"
#include "rlcsa.hpp"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "utils.h"

// #include "common.h"
#include <cstdint>

bool verbose = false;

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
      std::vector<unsigned int> nn = {n};
      intervals.push_back({sai_t, an.second, nn});
      // if (get_bwt_symb(index, an.first) == symb) {
      //   //  std::vector<unsigned int> nn = {n, an.second};
      //   std::vector<unsigned int> nn = {n};
      //   intervals.push_back({sai_t, an.second, nn});
      //   // std::fprintf(stderr, "size %ld \n",
      //   //         intervals[intervals.size() - 1].path[0]);
      // }
    }
  }

  return intervals;
}

std::vector<node_sai> ext_int(const rld_t *index, const uint8_t *s, int i,
                              unsigned int l, rldintv_t &sai,
                              std::vector<std::vector<unsigned int>> &tags,
                              std::vector<std::vector<unsigned int>> &adj,
                              std::vector<unsigned int> labels_map,
                              unsigned int curr_node) {

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

std::vector<node_sai> ext_alph(const rld_t *index, const uint8_t symb,
                               std::vector<std::vector<unsigned int>> &tags,
                               std::vector<std::vector<unsigned int>> &adj,
                               std::vector<unsigned int> labels_map,
                               std::vector<node_sai> int_s) {

  std::vector<node_sai> int_curr = int_s;

  std::vector<node_sai> int_next;

  for (auto ic : int_curr) {
    rldintv_t osai[6];
    rld_extend(index, &ic.sai, osai, 1);
    auto sai = osai[symb];
    // printf("int detail with %ld: ", symb);
    // printf("[%ld, %ld, %ld] to [%ld, %ld, %ld]\n", ic.sai.x[0], ic.sai.x[1],
    //        ic.sai.x[2], sai.x[0], sai.x[1], sai.x[2]);
    if (sai.x[2] > 0) {
      auto tmp_int =
          get_intervals(index, sai.x[0], sai.x[2], tags, adj, symb, {});
      int_next.insert(int_next.end(), tmp_int.begin(), tmp_int.end());
      int_next.push_back({sai, ic.curr_node, ic.path});
    }
  }
  int_curr = merge(int_next);
  if (int_curr.size() == 0) {
    return {};
  }
  return int_curr;
}

std::vector<std::vector<node_sai>>
ext_by_alph(const rld_t *index, std::vector<std::vector<unsigned int>> &tags,
            std::vector<std::vector<unsigned int>> &adj,
            std::vector<unsigned int> labels_map, std::vector<node_sai> int_s) {

  std::vector<std::vector<node_sai>> res(4);

  for (uint8_t i = 0; i < 4; i++) {
    uint8_t s = i + 1;
    // printf("extending with %ld: ", s);
    // for (auto eee : int_s) {
    //   printf("[%ld, %ld, %ld] ", eee.sai.x[0], eee.sai.x[1], eee.sai.x[2]);
    // }
    // std::cout << "\n";
    if (int_s.empty()) {
      res[i] = {};
      continue;
    }
    res[i] = ext_alph(index, s, tags, adj, labels_map, int_s);
    // printf("otaining ", s);
    // for (auto eee : res[i]) {
    //   printf("[%ld, %ld, %ld] ", eee.sai.x[0], eee.sai.x[1], eee.sai.x[2]);
    // }
    // std::cout << "\n";
  }
  return res;
}

auto findDuplicates(std::vector<std::string> &v) {
  std::sort(v.begin(), v.end());
  std::vector<std::string> res;
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i] == v[i - 1] and (res.empty() or res.back() != v[i]))
      res.push_back(v[i]);
  return res;
}

void build_bwt_rope(const char *gfa_file, std::string out_prefix, int threads) {
  std::string rope_out = out_prefix + ".bwt";
  std::string tag_out = out_prefix + ".dollars";
  std::string graph_out = out_prefix + ".graph";
  std::string labels_out = out_prefix + ".labels";
  std::string int_out = out_prefix + ".cache";
  gfa_t *ingfa;
  ingfa = gfa_read(gfa_file);
  std::vector<std::pair<std::string, unsigned int>> labels;

  double t_start = realtime();
  rlcsa_t *rlc = rlc_init();

  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  fflush(stderr);
  std::vector<std::vector<unsigned int>> adj(ingfa->n_seg);
  std::vector<unsigned int> labels_map(ingfa->n_seg);
  int insn = 0;
  // std::cout << "nodes: " << ingfa->n_seg << "\n";
  for (unsigned int i = 0; i < ingfa->n_seg; ++i) {
    auto l = ingfa->seg[i].len;
    std::string seq_t = ingfa->seg[i].seq;
    labels.push_back(std::make_pair(seq_t, i));

    // auto seq = seq_t.c_str();
    // s = (uint8_t *)seq;

    // //  change encoding
    // for (int j = 0; j < l; ++j)
    //   s[j] = fm6_i(s[j]);

    // kputsn((char *)s, l + 1, &buf);
    // ++insn;
    // if (insn == INT_MAX) { // FIXME: hardcoded
    //   // NOTE: insert works only if string is "long enough" (at least ~13)
    //   rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
    //   buf.l = 0;
    //   insn = 0;
    // }

    char *segname = ingfa->seg[i].name;
    int32_t segid = i;
    uint32_t vid = segid << 1;
    ++vid;
    gfa_arc_t *inca_vid = gfa_arc_a(ingfa, vid);
    uint32_t n_inca_vid = gfa_arc_n(ingfa, vid);
    for (unsigned int j = 0; j < n_inca_vid; ++j) {
      adj[i].push_back(gfa_name2id(ingfa, ingfa->seg[inca_vid[j].w >> 1].name));
    }
    labels_map[i] = std::stoi(segname);
  }

  std::vector<std::vector<unsigned int>> adj_f(ingfa->n_seg);
  std::vector<unsigned int> labels_map_f(ingfa->n_seg);

  std::stable_sort(labels.begin(), labels.end(),
                   [](const std::pair<std::string, unsigned int> &a,
                      const std::pair<std::string, unsigned int> &b) {
                     auto aa(a.first);
                     for (auto &x : aa)
                       x = (x != 'N' ? x : 'Z');
                     auto bb(b.first);
                     for (auto &x : bb)
                       x = (x != 'N' ? x : 'Z');
                     return aa < bb;
                   });
  // std::stable_sort(labels.begin(), labels.end(),
  //                  [](const std::pair<std::string, unsigned int> &a,
  //                     const std::pair<std::string, unsigned int> &b) {
  //                    // auto aa(a.first);
  //                    auto aa = std::string(a.first.rbegin(), a.first.rend());
  //                    for (auto &x : aa)
  //                      x = (x != 'N' ? x : 'Z');
  //                    // auto bb(b.first);
  //                    auto bb = std::string(b.first.rbegin(), b.first.rend());
  //                    for (auto &x : bb)
  //                      x = (x != 'N' ? x : 'Z');
  //                    return aa < bb;
  //                  });

  auto map_l = std::vector<unsigned int>(labels.size());
  auto j = 0;
  for (auto l : labels) {
    map_l[l.second] = j;
    j++;
  }
  for (unsigned int i = 0; i < adj.size(); i++) {
    for (unsigned int j = 0; j < adj[i].size(); j++) {
      adj_f[map_l[i]].push_back(map_l[adj[i][j]]);
    }
  }

  for (unsigned int i = 0; i < labels_map.size(); i++) {
    labels_map_f[map_l[i]] = labels_map[i];
  }

  for (auto ll : labels) {
    auto seq = ll.first.c_str();
    auto l = ll.first.size();
    s = (uint8_t *)seq;

    //  change encoding
    for (int j = 0; j < l; ++j)
      s[j] = fm6_i(s[j]);

    kputsn((char *)s, l + 1, &buf);
    ++insn;
    if (insn == INT_MAX) { // FIXME: hardcoded
      // NOTE: insert works only if string is "long enough" (at least ~13)
      rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
      buf.l = 0;
      insn = 0;
    }
  }
  // std::cout << "nodes: " << labels.size() << "\n";
  // int sl = 0;
  // for (auto l : labels) {
  //   sl += l.first.size();
  // }
  // std::cout << "nodes size labels " << sl << "\n";
  // std::cout << "nodes size buffer " << buf.l << "\n";
  if (insn) {
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
    // buf.l = 0;
  }
  int dol = 0;

  // std::fprintf(stderr, "stringbuffer: ");
  // for (unsigned int i = 0; i < buf.l; i++) {
  //   std::fprintf(stderr, "%c", "$ACGTN"[buf.s[i]]);
  // }

  // std::cout << "nodes size buffer dol " << dol << "\n";
  uintmat_dump(adj_f, graph_out.c_str());
  uintvec_dump(labels_map_f, labels_out.c_str());
  rlc_dump(rlc, rope_out.c_str());

  std::vector<std::vector<unsigned int>> tags(2);
  tags[0].resize(ingfa->n_seg);
  tags[1].resize(ingfa->n_seg);

  // rlc_print_bwt(rlc);

  std::fprintf(stderr,
               "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
               __func__, realtime() - t_start, cputime());
  rlc_destroy(rlc);

  // int c = 0;
  // for (auto &l : labels) {
  //   l.second = c;
  //   // std::cout << c << " " << l.first << " " << l.second << "\n";
  //   c++;
  // }

  // std::stable_sort(labels.begin(), labels.end(),
  //                  [](const std::pair<std::string, unsigned int> &a,
  //                     const std::pair<std::string, unsigned int> &b) {
  //                    auto aa(a.first);
  //                    for (auto &x : aa)
  //                      x = (x != 'N' ? x : 'Z');
  //                    auto bb(b.first);
  //                    for (auto &x : bb)
  //                      x = (x != 'N' ? x : 'Z');
  //                    return aa < bb;
  //                  });

  // c = 0;
  // for (auto l : labels) {
  //   std::cout << c << " " << l.first << " " << l.second << "\n";
  //   c++;
  // }
  //
  // int off = 0;
  // for (auto l : labels) {
  //   if (l.second == 0) {
  //     tags[0][off] = ingfa->n_seg - 1;
  //     tags[1][off] = off;
  //   } else {
  //     tags[0][off] = l.second - 1;
  //     tags[1][off] = off;
  //   }
  //   off++;
  // }
  for (int off = 0; off < labels.size(); off++) {
    if (off == 0) {
      tags[0][off] = ingfa->n_seg - 1;
      tags[1][off] = off;
    } else {
      tags[0][off] = off - 1;
      tags[1][off] = off;
    }
  }
  // for (auto t : tags[0]) {
  //   std::cout << t << "\n";
  // }
  uintmat_dump(tags, tag_out.c_str());
  gfa_destroy(ingfa);
  free(buf.s);

  // std::vector<std::string> st = {};
  // for (auto l : labels) {
  //   st.push_back(l.first);
  // }

  // auto d = findDuplicates(st);
  // for (auto dd : d) {
  //   std::cout << dd << "\n";
  // }
  //
  std::vector<std::vector<node_sai>> intervals(1024);
  std::vector<std::string> suff(1024);
  unsigned int ls = 5;

  for (int i = 0; i < 1024; ++i) {
    std::string suff_t(5, ' ');
    int cur = i;

    for (int j = 0; j < 5; j++) {
      suff_t[j] = "ACGT"[cur % 4];
      cur /= 4;
    }

    suff[i] = suff_t;
    // std::cout << suff[i] << "\n";
  }

  rld_t *index = rld_restore(rope_out.c_str());
  rldintv_t sai;
  for (int c = 0; c < 6; ++c) {
    fm6_set_intv(index, c, sai);
  }
  auto suff_b = suff;
  // #pragma omp parallel for
  //   for (int j = 0; j < suff.size(); j++) {
  //     s = (uint8_t *)suff[j].c_str();
  //     // change encoding
  //     for (int i = 0; i < 5; ++i) {
  //       s[i] = fm6_i(s[i]);
  //     }
  //     // for (int i = 0; i < 5; ++i) {
  //     //   printf("%d ", s[i]);
  //     // }
  //     std::cerr << j << "\r";
  //     // rldintv_t osai[6];
  //     int i = 4;
  //     fm6_set_intv(index, s[i], sai);
  //     --i;
  //     auto s_b = sai.x[1];
  //     sai.x[1] = tags[0].size();
  //     // auto nodes = ext(index, s, i, l, sai, tags, adj, tags[0].size());
  //     intervals[j] =
  //         ext_int(index, s, i, 5, sai, tags, adj, labels_map,
  //         tags[0].size());
  //     // std::fprintf(stderr, "finished");
  //   }
  //   std::vector<std::vector<std::vector<unsigned int>>> intervals_f(1024);
  //   int k = 0;
  //   for (auto s : intervals) {
  //     std::cout << k << " " << suff_b[k] << ": ";
  //     std::vector<std::vector<unsigned int>> tmp_v;
  //     for (auto i : s) {
  //       std::vector<unsigned int> t = {0, i.sai.x[0], i.sai.x[1], i.sai.x[2],
  //                                      i.curr_node};
  //       tmp_v.push_back(t);
  //       std::cout << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << " "
  //                 << t[4] << "\n";
  //     }
  //     // if (!tmp_v.empty())
  //     std::cout << std::endl;
  //     intervals_f[k] = tmp_v;
  //     k++;
  //   }

  //   for (auto s : intervals_f) {
  //     // std::cout << s.size() << "\n";
  //     for (auto t : s) {
  //       std::cout << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << " "
  //                 << t[4] << "\n";
  //     }
  //     if (!s.empty())
  //       std::cout << std::endl;
  //   }

  //   uintmat3_dump(intervals_f, int_out.c_str());

  std::vector<std::vector<std::vector<node_sai>>> intervals_fast(5);

  for (int c = 1; c < 5; ++c) {
    fm6_set_intv(index, c, sai);
    // printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);

    std::vector<node_sai> e = {{sai, tags[0].size()}};
    intervals_fast[4].push_back(e);
    auto tmp_int = get_intervals(index, sai.x[0], sai.x[2], tags, adj, c, {});
    intervals_fast[4][c - 1].insert(intervals_fast[4][c - 1].end(),
                                    tmp_int.begin(), tmp_int.end());
    // printf("%d: [%ld, %ld, %ld]\n", c,
    // intervals_fast[4].back()[0].sai.x[0],
    //        intervals_fast[4].back()[0].sai.x[1],
    //        intervals_fast[4].back()[0].sai.x[2]);
  }
  // std::cerr << "Preparing cache, this will take time\n";
  for (int i = 3; i >= 0; i--) {

    // for (int j = intervals_fast[i + 1].size() - 1; j >= 0; j--) {
    for (int j = 0; j < intervals_fast[i + 1].size(); j++) {
      // if (intervals_fast[i + 1][j].empty()) {
      //   intervals_fast[i].push_back({});
      //   continue;
      // }
      auto e =
          ext_by_alph(index, tags, adj, labels_map_f, intervals_fast[i + 1][j]);
      // std::cerr << e.size() << "\n";
      //  intervals_fast[i].insert(intervals_fast[i].end, e.begin(),e.end());
      for (auto ee : e) {
        intervals_fast[i].push_back(ee);
        // for (auto eee : ee)
        //   printf("[%ld, %ld, %ld] ", eee.sai.x[0], eee.sai.x[1],
        //   eee.sai.x[2]);
        // std::cout << "\n";
      }
    }
    std::cerr << intervals_fast[i].size() << "\n";
  }
  std::vector<std::vector<std::vector<unsigned int>>> intervals_f(1024);
  int k = 0;
  for (auto s : intervals_fast[0]) {
    // std::cout << k << " " << suff[k] << ": ";
    std::vector<std::vector<unsigned int>> tmp_v;
    for (auto i : s) {
      std::vector<unsigned int> t = {0, i.sai.x[0], i.sai.x[1], i.sai.x[2],
                                     i.curr_node};
      tmp_v.push_back(t);
      // std::cout << t[0] << " " << t[1] << " " << t[2] << " " << t[3]
      //           << " "
      //           //
      //           << t[4] << "\n";
    }
    // if (!tmp_v.empty())
    // std::cout << std::endl;
    intervals_f[k] = tmp_v;
    k++;
  }
  uintmat3_dump(intervals_f, int_out.c_str());

  // std::cout << "\n---------------------\n";
  // std::vector<std::vector<std::vector<unsigned int>>> c_int =
  //     uintmat3_load(int_out.c_str());

  // std::vector<std::vector<node_sai>> cache;

  // for (auto s : c_int) {

  //   for (auto t : s) {
  //     std::cout << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << " "
  //               << t[4] << "\n";
  //   }
  //   if (!s.empty())
  //     std::cout << std::endl;
  // }
  // for (auto s : c_int) {
  //   std::vector<node_sai> tv;
  //   rldintv_t sai;
  //   unsigned int curr_node;
  //   for (auto v : s) {
  //     sai.info = v[0];
  //     sai.x[0] = v[1];
  //     sai.x[1] = v[2];
  //     sai.x[2] = v[3];
  //     // tv.push_back({sai, v[4], {v[4]}});
  //     std::cout << sai.x[0] << " " << sai.x[1] << " " << sai.x[2] << " " <<
  //     v[4]
  //               << "\n";
  //   }
  //   if (!s.empty())
  //     std::cout << std::endl;
  //   cache.push_back(tv);
  // }
}

#endif // GRAPHINDEX_BWT_ROPE_H
