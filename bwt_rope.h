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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "utils.h"

#include "common.h"
#include <cstdint>

std::vector<node_sai> ext_alph(const rld_t *index, const uint8_t symb,
                               std::vector<std::vector<unsigned int>> &tags,
                               std::vector<std::vector<unsigned int>> &adj,
                               std::vector<unsigned int> labels_map,
                               std::vector<node_sai> int_s) {

  std::vector<node_sai> int_next;
  for (auto ic : int_s) {
    rldintv_t osai[6];
    rld_extend(index, &ic.sai, osai, 1);
    auto sai = osai[symb];

    if (sai.x[2] > 0) {
      auto tmp_int = get_intervals(index, sai.x[0], sai.x[2], tags, adj, 6, {});
      if (!tmp_int.empty()) {
        int_next.insert(int_next.end(), tmp_int.begin(), tmp_int.end());
      }
      int_next.push_back({sai, ic.curr_node, ic.path});
    }
  }
  auto int_curr = merge(int_next);
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
  // #pragma omp parallel
  // #pragma omp for
#pragma omp parallel for num_threads(4)
  for (uint8_t i = 0; i < 4; i++) {
    uint8_t s = i + 1;

    if (int_s.empty()) {
      res[i] = {};
      continue;
    }
    res[i] = ext_alph(index, s, tags, adj, labels_map, int_s);
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

void build_bwt_rope(const char *gfa_file, std::string out_prefix, int threads,
                    int cache) {
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

  if (insn) {
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
  }
  int dol = 0;

  uintmat_dump(adj_f, graph_out.c_str());
  uintvec_dump(labels_map_f, labels_out.c_str());
  rlc_dump(rlc, rope_out.c_str());

  std::vector<std::vector<unsigned int>> tags(2);
  tags[0].resize(ingfa->n_seg);
  tags[1].resize(ingfa->n_seg);

  std::fprintf(stderr,
               "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
               __func__, realtime() - t_start, cputime());
  rlc_destroy(rlc);

  for (int off = 0; off < labels.size(); off++) {
    if (off == 0) {
      tags[0][off] = ingfa->n_seg - 1;
      tags[1][off] = off;
    } else {
      tags[0][off] = off - 1;
      tags[1][off] = off;
    }
  }

  uintmat_dump(tags, tag_out.c_str());
  gfa_destroy(ingfa);
  free(buf.s);

  rld_t *index = rld_restore(rope_out.c_str());
  rldintv_t sai;
  for (int c = 0; c < 6; ++c) {
    fm6_set_intv(index, c, sai);
  }

  std::vector<std::vector<std::vector<node_sai>>> intervals_fast(cache);

  // #pragma omp parallel
  // #pragma omp for
#pragma omp parallel for num_threads(4)
  for (int c = 1; c < 5; ++c) {
    fm6_set_intv(index, c, sai);

    std::vector<node_sai> e = {{sai, tags[0].size()}};
    intervals_fast[cache - 1].push_back(e);
    auto tmp_int = get_intervals(index, sai.x[0], sai.x[2], tags, adj_f, 6, {});
    intervals_fast[cache - 1][c - 1].insert(
        intervals_fast[cache - 1][c - 1].end(), tmp_int.begin(), tmp_int.end());
  }
  std::cerr << "Preparing cache, this will take time\n";
  for (int i = cache - 2; i >= 0; i--) {
    for (int j = 0; j < intervals_fast[i + 1].size(); j++) {

      auto e = ext_by_alph(index, tags, adj_f, labels_map_f,
                           intervals_fast[i + 1][j]);
      intervals_fast[i + 1].clear();
      for (auto ee : e) {
        intervals_fast[i].push_back(ee);
      }
    }
    std::cerr << intervals_fast[i].size() << "\n";
  }
  std::vector<std::vector<std::vector<unsigned int>>> intervals_f(
      (int)std::pow(4, cache));
  int k = 0;
  for (auto s : intervals_fast[0]) {
    std::vector<std::vector<unsigned int>> tmp_v;
    for (auto i : s) {
      std::vector<unsigned int> t = {i.sai.x[0], i.sai.x[2], i.curr_node};
      tmp_v.push_back(t);
    }

    intervals_f[k] = tmp_v;
    k++;
  }
  uintmat3_dump(intervals_f, int_out.c_str());
}

#endif // GRAPHINDEX_BWT_ROPE_H
