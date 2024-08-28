//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_H
#define GRAPHINDEX_BWT_ROPE_H

#include "gfa.h"
#include "kseq.h"
#include "mrope.h"
#include "rld0.h"
#include "rle.h"
#include <algorithm>
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
  // std::vector<std::string> labels;
  int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES,
      so = MR_SO_IO;
  // std::cerr << block_len << ". " << max_nodes << ", " << so << "\n";

  int64_t m = (int64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1;
  double t_start = realtime();
  mrope_t *mr = mr_init(max_nodes, block_len, so);
  if (threads > 0)
    mr_thr_min(mr, 10);
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  double ct, rt;
  // fprintf(stderr, "nodes: %d\n", ingfa->n_seg);
  fflush(stderr);
  bool start = true;
  // std::vector<unsigned int> T;
  std::vector<std::vector<unsigned int>> adj(ingfa->n_seg);
  for (unsigned int i = 0; i < ingfa->n_seg; ++i) {
    auto l = ingfa->seg[i].len;
    std::string seq_t = ingfa->seg[i].seq;
    labels.push_back(std::make_pair(seq_t, i));
    // labels.push_back(seq_t);
    std::reverse(seq_t.begin(), seq_t.end());

    // if (start) {
    //   seq_t = "A" + seq_t;
    //   l++;
    // }
    auto seq = seq_t.c_str();
    s = (uint8_t *)seq;

    //  change encoding
    for (unsigned int j = 0; j < l; ++j) {
      s[j] = fm6_i(s[j]);
      // if (start && j == 0) {
      //   s[j] = 0;
      //   start = false;
      // }
    }
    // // std::cout << std::endl;
    // for (unsigned int j = 0; j < l; ++j) {
    //   fprintf(stderr, "%d ", s[j]);
    // }
    // std::cout << std::endl;
    kputsn((char *)s, l + 1, &buf);

    // std::reverse(seq_t.begin(), seq_t.end());
    // for (unsigned int j = 0; j < l; ++j) {
    //   s[j] = fm6_comp(s[j]);
    //   // if (start && j == 0) {
    //   //   s[j] = 0;
    //   //   start = false;
    //   // }
    // }

    // kputsn((char *)s, l + 1, &buf);

    // if (l & 1)
    //   s[i] = fm6_comp(s[i]);
    // kputsn((char *)s, l + 1, &buf);

    mr_insert1(mr, (const uint8_t *)s);

    // if (buf.l >= m) {
    //   ct = cputime(), rt = realtime();

    //   // mr_insert_multi(mr, buf.l, (const uint8_t *)buf.s, threads);
    //   fprintf(stderr,
    //           "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
    //           __func__, (long)buf.l, realtime() - rt, cputime() - ct);
    //   buf.l = 0;
    // }

    char *segname = ingfa->seg[i].name;
    int32_t segid = i;
    uint32_t vid = segid << 1;
    ++vid;
    gfa_arc_t *inca_vid = gfa_arc_a(ingfa, vid);
    uint32_t n_inca_vid = gfa_arc_n(ingfa, vid);
    for (int j = 0; j < n_inca_vid; ++j) {
      adj[i].push_back(gfa_name2id(ingfa, ingfa->seg[inca_vid[j].w >> 1].name));
    }
    // break;
  }
  // for (unsigned int j = 0; j < buf.l; ++j) {
  //   fprintf(stderr, "%c", map[buf.s[j]]);
  // }
  // std::cout << std::endl;
  // int cc = 0;
  // for (auto v : adj) {
  //   std::cout << cc << ": ";
  //   for (auto n : v) {
  //     std::cout << n << " ";
  //   }
  //   std::cout << "\n";
  //   cc++;
  // }
  uintmat_dump(adj, graph_out.c_str());

  // for (auto tt : T) {
  //   std::cerr << tt;
  // }
  // std::cout << std::endl;
  // std::sort(labels.begin(), labels.end());
  // std::vector<unsigned int> tags(ingfa->n_seg);
  std::vector<std::vector<unsigned int>> tags(2);
  tags[0].resize(ingfa->n_seg);
  tags[1].resize(ingfa->n_seg);

  // dump index to stdout
  mritr_t itr;
  const uint8_t *block;
  rld_t *e = 0;
  rlditr_t di;
  e = rld_init(6, 3);
  rld_itr_init(e, &di, 0);
  mr_itr_first(mr, &itr, 1);
  int i = 0;
  // int sum = 0;
  while ((block = mr_itr_next_block(&itr)) != 0) {
    const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
    while (q < end) {
      int c = 0;
      int64_t l;
      rle_dec1(q, c, l);

      // fprintf(stderr, "%c, %d\n", map[c], l);
      //     sum += l;
      rld_enc(e, &di, l, c);
    }
    // std::cout << "\n";
  }
  // fprintf(stderr, "tot: %d ", sum);
  rld_enc_finish(e, &di);
  rld_dump(e, robe_out.c_str());
  fprintf(stderr,
          "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
          __func__, realtime() - t_start, cputime());
  mr_destroy(mr);

  std::stable_sort(labels.begin(), labels.end());

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
  for (auto t : tags[0]) {
    std::cout << t << "\n";
  }

  uintmat_dump(tags, tag_out.c_str());
  gfa_destroy(ingfa);
  free(buf.s);
}

#endif // GRAPHINDEX_BWT_ROPE_H
