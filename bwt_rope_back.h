//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_H
#define GRAPHINDEX_BWT_ROPE_H

#include "gfa.h"
#include "kseq.h"
#include "rlcsa.hpp"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "utils.h"

#include <cstdint>

void build_bwt_rope(const char *gfa_file, std::string out_prefix, int threads) {
  std::string robe_out = out_prefix + ".ser";
  std::string tag_out = out_prefix + ".tag";
  std::string graph_out = out_prefix + ".graph";
  gfa_t *ingfa;
  ingfa = gfa_read(gfa_file);
  std::vector<std::pair<std::string, unsigned int>> labels;

  int64_t m = (int64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1;
  double t_start = realtime();
  rlcsa_t *rlc = rlc_init();

  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  double ct, rt;
  fflush(stderr);
  bool start = true;
  std::vector<std::vector<unsigned int>> adj(ingfa->n_seg);
  int insn = 0;
  for (unsigned int i = 0; i < ingfa->n_seg; ++i) {
    auto l = ingfa->seg[i].len;
    std::string seq_t = ingfa->seg[i].seq;
    labels.push_back(std::make_pair(seq_t, i));

    auto seq = seq_t.c_str();
    s = (uint8_t *)seq;

    //  change encoding
    for (unsigned int j = 0; j < l; ++j)
      s[j] = fm6_i(s[j]);

    kputsn((char *)s, l + 1, &buf);
    ++insn;
    if (insn == INT_MAX) { // FIXME: hardcoded
      // NOTE: insert works only if string is "long enough" (at least ~13)

      fprintf(stderr, "buf.l %d ", buf.l);
      for (unsigned int i = 0; i < buf.l; i++) {
        fprintf(stderr, "%d ", buf.s[i]);
      }
      fprintf(stderr, "\n");
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
    for (int j = 0; j < n_inca_vid; ++j) {
      adj[i].push_back(gfa_name2id(ingfa, ingfa->seg[inca_vid[j].w >> 1].name));
    }
  }
  for (unsigned int i = 0; i < buf.l; i++) {
    fprintf(stderr, "%d ", buf.s[i]);
  }
  if (insn) {

    fprintf(stderr, "buf.l %d ", buf.l);
    for (unsigned int i = 0; i < buf.l; i++) {
      fprintf(stderr, "%d ", buf.s[i]);
    }
    fprintf(stderr, "\n");
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);
    buf.l = 0;
  }

  uintmat_dump(adj, graph_out.c_str());

  std::vector<std::vector<unsigned int>> tags(2);
  tags[0].resize(ingfa->n_seg);
  tags[1].resize(ingfa->n_seg);

  rlc_print_bwt(rlc);
  rlc_dump(rlc, robe_out.c_str());
  fprintf(stderr,
          "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
          __func__, realtime() - t_start, cputime());
  rlc_destroy(rlc);

  std::stable_sort(labels.begin(), labels.end());
  for (auto l : labels) {
    std::cout << l.first << " " << l.second << "\n";
  }
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

  uintmat_dump(tags, tag_out.c_str());

  for (auto tag : tags[0]) {
    std::cout << tag << " ";
  }
  std::cout << std::endl;
  for (auto tag : tags[1]) {
    std::cout << tag << " ";
  }
  std::cout << std::endl;

  rld_t *index = rld_restore(robe_out.c_str());

  rldintv_t sai;
  for (int c = 0; c < 6; ++c) {
    fm6_set_intv(index, c, sai);
    printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  }

  std::vector<std::vector<unsigned int>> tags2(2);
  tags2[0].resize(ingfa->n_seg, ingfa->n_seg);
  tags2[1].resize(ingfa->n_seg, ingfa->n_seg);

  for (auto tag : tags2[0]) {
    std::cout << tag << " ";
  }
  std::cout << std::endl;
  rldintv_t osai[6];
  int i = 0;
  unsigned int t = 0;

  for (unsigned int j = 0; j < ingfa->n_seg; ++j) {
    auto l = ingfa->seg[j].len + 2;
    std::string seq_t = ingfa->seg[j].seq;
    seq_t = "$" + seq_t + "$";
    // auto l = ingfa->seg[j].len + 1;
    // std::string seq_t = ingfa->seg[j].seq;
    // seq_t = seq_t + "$";
    std::cout << seq_t;
    auto ss = (uint8_t *)seq_t.c_str();
    ss[0] = 0;
    ss[l - 1] = 0;
    // change encoding
    for (i = 0; i < l; ++i)
      ss[i] = ss[i] < 128 ? seq_nt6_table[ss[i]] : 5;

    // init search
    i = l - 1;
    fm6_set_intv(index, ss[i], sai);
    --i;

    // backward extensions
    for (; i >= 0; --i) {
      rld_extend(index, &sai, osai, 1);
      sai = osai[ss[i]];
      if (sai.x[2] <= 0) {
        // errors += 1;
        break;
      }
    }

    printf("final: [%ld, %ld, %ld] \n", sai.x[0], sai.x[1], sai.x[2]);
    auto dol = rld_rank11(index, sai.x[0], 0);

    if (t == 0) {
      tags2[0][dol] = ingfa->n_seg - 1;
      tags2[1][ingfa->n_seg - 1] = dol;
    } else {
      while (tags2[0][dol] != ingfa->n_seg) {
        dol++;
      }
      tags2[0][dol] = t - 1;
      tags2[1][t - 1] = dol;
    }
    printf("dol = %d, tag[dol] =%d\n", dol, tags2[0][dol]);
    // printf("%d, %d\n", dol, t);
    t++;
  }
  for (auto tag : tags2[0]) {
    std::cout << tag << " ";
  }
  std::cout << std::endl;
  for (auto tag : tags2[1]) {
    std::cout << tag << " ";
  }
  std::cout << std::endl;

  gfa_destroy(ingfa);
  free(buf.s);
}

#endif // GRAPHINDEX_BWT_ROPE_H
