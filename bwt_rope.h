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
      so = MR_SO_RLO;
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
  // int ii = 0;
  // for (auto &la : labels) {
  //   std::cerr << la.first << " " << la.second << "\n";
  //   tags[ii] = la.second;
  //   ii++;
  // }

  // uintvec_dump(tags, tag_out.c_str());
  // for (unsigned int i = 0; i < buf.l; i++) {
  //   fprintf(stderr, "%d ", buf.s[i]);
  // }
  //  fprintf(stderr, "%d\n", buf.s[0] == buf.s[buf.l - 1]);
  //  if (buf.l) {
  //    ct = cputime(), rt = realtime();
  //    mr_insert_multi(mr, buf.l, (const uint8_t *)buf.s, threads);
  //    fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU
  //    sec\n",
  //            __func__, (long)buf.l, realtime() - rt, cputime() - ct);
  //    buf.l = 0;
  //  }

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

  // rld_t *index = rld_restore(robe_out.c_str());

  // rldintv_t sai;
  // for (int c = 0; c < 6; ++c) {
  //   fm6_set_intv(index, c, sai);
  //   // printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  // }
  // rldintv_t osai[6];
  // i = 0;
  // unsigned int t = 0;
  // for (unsigned int j = 0; j < ingfa->n_seg; ++j) {
  //   auto l = ingfa->seg[j].len;
  //   std::string seq_t = ingfa->seg[j].seq;
  //   // std::cout << seq_t;
  //   auto ss = (uint8_t *)seq_t.c_str();

  //   // change encoding
  //   for (i = 0; i < l; ++i)
  //     ss[i] = ss[i] < 128 ? seq_nt6_table[ss[i]] : 5;

  //   // init search
  //   i = l - 1;
  //   fm6_set_intv(index, ss[i], sai);
  //   --i;

  //   // backward extensions
  //   for (; i >= 0; --i) {
  //     rld_extend(index, &sai, osai, 1);
  //     sai = osai[ss[i]];
  //     if (sai.x[2] <= 0) {
  //       // errors += 1;
  //       break;
  //     }
  //   }

  //   // printf("final: [%ld, %ld, %ld] \n", sai.x[0], sai.x[1], sai.x[2]);
  //   auto dol = rld_rank11(index, sai.x[0], 0);

  //   if (t == 0) {
  //     tags[0][dol] = ingfa->n_seg - 1;
  //     tags[1][ingfa->n_seg - 1] = dol;
  //   } else {
  //     tags[0][dol] = t - 1;
  //     tags[1][t - 1] = dol;
  //   }
  //   // printf("dol = %d, tag[dol] =%d\n", dol, tags[0][dol]);
  //   //  printf("dol = %d, tag[dol] =%d\n", dol, tags[dol]);
  //   //   printf("%d, %d\n", dol, t);
  //   t++;
  // }

  // for (auto tt : tags[0]) {
  //   std::cout << tt << " ";
  // }
  // std::cout << "\n";
  // for (auto tt : tags[1]) {
  //   std::cout << tt << " ";
  // }
  // for (auto tt : tags) {
  //   std::cout << tt << " ";
  // }
  std::sort(labels.begin(), labels.end());
  // for (auto l : labels) {
  //   std::cout << l.first << " " << l.second << "\n";
  // }
  int off = 0;
  for (auto l : labels) {
    // tags[0][off] = l.second;
    // tags[1][l.second] = off;

    if (l.second == 0) {
      tags[0][off] = ingfa->n_seg - 1;
      tags[1][ingfa->n_seg - 1] = off;
    } else {
      tags[0][off] = l.second - 1;
      tags[1][l.second - 1] = off;
      //   }
    }
    off++;
  }
  // for (auto tt : tags[0]) {
  //   std::cout << tt << " ";
  // }
  // std::cout << "\n";
  // for (auto tt : tags[1]) {
  //   std::cout << tt << " ";
  // }

  uintmat_dump(tags, tag_out.c_str());
  gfa_destroy(ingfa);
  free(buf.s);
}

#endif // GRAPHINDEX_BWT_ROPE_H
