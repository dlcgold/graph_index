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
                               std::vector<std::vector<uint64_t>> &tags,
                               std::vector<std::vector<uint64_t>> &adj,
                               std::vector<uint64_t> labels_map,
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
ext_by_alph(const rld_t *index, std::vector<std::vector<uint64_t>> &tags,
            std::vector<std::vector<uint64_t>> &adj,
            std::vector<uint64_t> labels_map, std::vector<node_sai> int_s) {

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
  std::vector<std::pair<std::string, uint64_t>> labels;

  double t_start = realtime();
  rlcsa_t *rlc = rlc_init();

  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  fflush(stderr);
  std::vector<std::vector<uint64_t>> adj(ingfa->n_seg);
  std::vector<uint64_t> labels_map(ingfa->n_seg + 1);
  int insn = 0;
  uint64_t ml = 0;
  uint64_t al = 0;
  uint64_t ca = 0;
  bool fc = true;
  // std::cout << "nodes: " << ingfa->n_seg << "\n";
  std::cout << "Loading the graph \n";
  for (uint64_t i = 0; i < ingfa->n_seg; ++i) {
    // auto l = ingfa->seg[i].len;
    std::string seq_t = ingfa->seg[i].seq;
    transform(seq_t.begin(), seq_t.end(), seq_t.begin(), ::toupper);
    labels.push_back(std::make_pair(seq_t, i));
    if (seq_t < std::string("ACGTN")) {
      al += seq_t.size() + 1;
    }
    if (fc && al > INT32_MAX) {
      fc = false;
      for (char ch : seq_t) {
        if (ch == 'A') {
          ca++;
        } else {
          break;
        }
      }
    }
    char *segname = ingfa->seg[i].name;
    int32_t segid = i;
    uint32_t vid = segid << 1;
    ++vid;
    gfa_arc_t *inca_vid = gfa_arc_a(ingfa, vid);
    uint32_t n_inca_vid = gfa_arc_n(ingfa, vid);
    for (uint64_t j = 0; j < n_inca_vid; ++j) {
      adj[i].push_back(gfa_name2id(ingfa, ingfa->seg[inca_vid[j].w >> 1].name));
    }
    labels_map[i] = std::stoi(segname);
    if (labels_map[i] > ml)
      ml = labels_map[i];
    // labels_map[i] = std::string(segname);
  }

  std::string dummy_s = "ACGTN";
  if (al + 1 > INT32_MAX) {
    std::cerr << "Possible issues related to BWT index. Trying fioxing with " << ca + 2 << " A\n";
    dummy_s = std::string(ca + 2, 'A') + std::string("CGTN");
  } else {
    std::cerr << "Safe range of " << al << " bases\n";
  }
  std::vector<std::vector<uint64_t>> adj_f(ingfa->n_seg + 1);
  std::vector<uint64_t> labels_map_f(ingfa->n_seg + 1);
  std::cout << "Sorting the labels\n";
  labels.insert(labels.begin(), {dummy_s, ingfa->n_seg});
  labels_map[ingfa->n_seg] = ml + 1;
  std::stable_sort(labels.begin(), labels.end(),
                   [](const std::pair<std::string, uint64_t> &a,
                      const std::pair<std::string, uint64_t> &b) {
                     auto aa(a.first);
                     for (auto &x : aa)
                       x = (x != 'N' ? x : 'Z');
                     auto bb(b.first);
                     for (auto &x : bb)
                       x = (x != 'N' ? x : 'Z');
                     return aa < bb;
                   });
  // labels.push_back({"$ACGTN", ingfa->n_seg});
  // labels.insert(labels.begin(), {"$ACGTN", ingfa->n_seg});
  // labels_map[ingfa->n_seg] = ml + 1;

  std::cout << "Build the index\n";
  auto map_l = std::vector<uint64_t>(labels.size());
  auto j = 0;
  for (auto l : labels) {
    map_l[l.second] = j;
    j++;
    // std::cout << l.first << " " << l.second << "\n";
  }
  // for (unsigned int i = 0; i < map_l.size(); i++) {
  //   printf("%d\n", map_l[i]);
  // }

  for (uint64_t i = 0; i < adj.size(); i++) {
    for (uint64_t j = 0; j < adj[i].size(); j++) {
      adj_f[map_l[i]].push_back(map_l[adj[i][j]]);
    }
  }

  // printf("done\n");
  for (uint64_t i = 0; i < labels_map.size(); i++) {
    labels_map_f[map_l[i]] = labels_map[i];
  }

  // for (unsigned int i = 0; i < labels_map.size(); i++) {
  //   printf("%d\n", labels_map[i]);
  // }
  // printf("\n");
  // for (unsigned int i = 0; i < labels_map_f.size(); i++) {
  //   printf("%d\n", labels_map_f[i]);
  // }
  // int cc = 0;
  // for (auto v : adj_f) {
  //   std::cout << cc << ": ";
  //   for (auto n : v) {
  //     std::cout << n << " ";
  //   }
  //   std::cout << "\n";
  //   cc++;
  // }

  int start = true;
  for (auto ll : labels) {
    auto seq = ll.first.c_str();
    auto l = ll.first.size();
    // std::cout << seq << "\n";
    s = (uint8_t *)seq;
    if (seq[0] == '$') {
      start = false;
      s[0] = 0;
    }
    //  change encoding
    for (uint64_t j = 0; j < l; ++j)
      s[j] = fm6_i(s[j]);

    //++insn;
    if (buf.l + l > INT32_MAX) { // FIXME: hardcoded
      // NOTE: insert works only if string is "long enough" (at least ~13)
      // std::fprintf(stderr, "stringbuffer: ");
      // for (uint64_t i = 0; i < buf.l; i++) {
      //    std::fprintf(stderr, "%c", "$ACGTN"[buf.s[i]]);
      // }
      rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
      fprintf(stderr, "[M::%s] inserted %ld symbols\n", __func__, (long)buf.l);
      buf.l = 0;
      insn = 0;
    }
    kputsn((char *)s, l + 1, &buf);
    ++insn;
  }

  if (buf.l) {
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
    fprintf(stderr, "[M::%s] inserted %ld symbols\n", __func__, (long)buf.l);
  }
  // int dol = 0;
  // for (unsigned int i = 0; i < buf.l; i++) {
  //   printf("%c", "$ACGTN"[buf.s[i]]);
  // }
  // printf("\n");
  uintmat_dump(adj_f, graph_out.c_str());
  uintvec_dump(labels_map_f, labels_out.c_str());
  // rlc_print_bwt(rlc);

  rlc_dump(rlc, rope_out.c_str());
  int bwt_l = rlc->l;

  std::vector<std::vector<uint64_t>> tags(2);
  tags[0].resize(ingfa->n_seg + 1);
  tags[1].resize(ingfa->n_seg + 1);

  std::fprintf(stderr,
               "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
               __func__, realtime() - t_start, cputime());
  rlc_destroy(rlc);

  for (uint64_t off = 0; off < labels.size(); off++) {
    if (off == 0) {
      // tags[0][off] = ingfa->n_seg - 1;
      tags[0][off] = ingfa->n_seg;
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
  // for (unsigned int i = 0; i < bwt_l; i++) {
  //   printf("%c", "$ACGTN"[get_bwt_symb(index, i)]);
  // }
  std::vector<std::vector<std::vector<node_sai>>> intervals_fast(cache);

  // #pragma omp parallel
  // #pragma omp for
  // #pragma omp parallel for num_threads(4)
  for (int c = 1; c < 5; ++c) {
    fm6_set_intv(index, c, sai);

    std::vector<node_sai> e = {{sai, tags[0].size()}};
    intervals_fast[cache - 1].push_back(e);
    auto tmp_int = get_intervals(index, sai.x[0], sai.x[2], tags, adj_f, 6, {});
    intervals_fast[cache - 1][c - 1].insert(
        intervals_fast[cache - 1][c - 1].end(), tmp_int.begin(), tmp_int.end());
  }
  std::cerr << "\nPreparing cache, this will take time\n";
  for (int i = cache - 2; i >= 0; i--) {
    for (uint64_t j = 0; j < intervals_fast[i + 1].size(); j++) {

      auto e = ext_by_alph(index, tags, adj_f, labels_map_f,
                           intervals_fast[i + 1][j]);
      // intervals_fast[i + 1].clear();
      for (auto ee : e) {
        intervals_fast[i].push_back(ee);
      }
    }
    intervals_fast[i + 1].clear();
    std::cerr << intervals_fast[i].size() << "\n";
  }
  std::vector<std::vector<std::vector<uint64_t>>> intervals_f(
      (int)std::pow(4, cache));
  int k = 0;
  for (auto s : intervals_fast[0]) {
    std::vector<std::vector<uint64_t>> tmp_v;
    for (auto i : s) {
      std::vector<uint64_t> t = {i.sai.x[0], i.sai.x[2], i.curr_node};
      tmp_v.push_back(t);

      // printf("%ld \t %ld \t %ld\n", i.sai.x[0], i.sai.x[2], i.curr_node);
    }

    intervals_f[k] = tmp_v;
    k++;
  }
  uintmat3_dump(intervals_f, int_out.c_str());
}

#endif // GRAPHINDEX_BWT_ROPE_H
