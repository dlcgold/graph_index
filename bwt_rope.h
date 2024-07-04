//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_H
#define GRAPHINDEX_BWT_ROPE_H

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

#include <cstdint>

auto findDuplicates(std::vector<std::string> &v) {
  std::sort(v.begin(), v.end());
  std::vector<std::string> res;
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i] == v[i - 1] and (res.empty() or res.back() != v[i]))
      res.push_back(v[i]);
  return res;
}

void build_bwt_rope(const char *gfa_file, std::string out_prefix, int threads) {
  std::string robe_out = out_prefix + ".bwt";
  std::string tag_out = out_prefix + ".dollars";
  std::string graph_out = out_prefix + ".graph";
  std::string labels_out = out_prefix + ".labels";
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

    auto seq = seq_t.c_str();
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
  // std::cout << "nodes: " << labels.size() << "\n";
  int sl = 0;
  for (auto l : labels) {
    sl += l.first.size();
  }
  // std::cout << "nodes size labels " << sl << "\n";
  // std::cout << "nodes size buffer " << buf.l << "\n";
  if (insn) {
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
    // buf.l = 0;
  }
  // int dol = 0;
  //  std::fprintf(stderr, "stringbuffer: ");
  //  for (unsigned int i = 0; i < buf.l; i++) {
  //    std::fprintf(stderr, "%c", "$ACGTN"[buf.s[i]]);
  //  }
  // std::cout << "nodes size buffer dol " << dol << "\n";
  uintmat_dump(adj, graph_out.c_str());
  uintvec_dump(labels_map, labels_out.c_str());
  std::vector<std::vector<unsigned int>> tags(2);
  tags[0].resize(ingfa->n_seg);
  tags[1].resize(ingfa->n_seg);

  // rlc_print_bwt(rlc);
  rlc_dump(rlc, robe_out.c_str());
  std::fprintf(stderr,
               "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
               __func__, realtime() - t_start, cputime());
  rlc_destroy(rlc);

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
  // int c = 0;
  // for (auto l : labels) {
  //   std::cout << c << " " << l.first << " " << l.second << "\n";
  //   c++;
  // }
  int off = 0;
  for (auto l : labels) {
    if (l.second == 0) {
      tags[0][off] = ingfa->n_seg - 1;
      tags[1][off] = off;
    } else {
      tags[0][off] = l.second - 1;
      tags[1][off] = off;
    }
    off++;
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
}

#endif // GRAPHINDEX_BWT_ROPE_H
