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

std::string rc(const std::string &s) {
  std::string rs = "";
  for (char c : s) {
    switch (c) {
    case 'A':
      rs += 'T';
      break;
    case 'C':
      rs += 'G';
      break;
    case 'G':
      rs += 'C';
      break;
    case 'T':
      rs += 'A';
      break;
    case 'N':
      rs += 'N';
      break;
    default:
      break;
    }
  }
  std::reverse(rs.begin(), rs.end());
  return rs;
}

void build_bwt_rope(const char *gfa_file, std::string out_prefix, int threads) {
  std::string rope_out = out_prefix + ".bwt";
  std::string tag_out = out_prefix + ".dollars";
  std::string graph_out = out_prefix + ".graphi";
  std::string graph_out2 = out_prefix + ".grapho";
  std::string labels_out = out_prefix + ".labels";
  gfa_t *ingfa;
  ingfa = gfa_read(gfa_file);
  std::vector<std::pair<std::string, unsigned int>> labels(ingfa->n_seg * 2);

  double t_start = realtime();
  rlcsa_t *rlc = rlc_init();

  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  fflush(stderr);
  std::vector<std::vector<unsigned int>> adj(ingfa->n_seg * 2);
  std::vector<std::vector<unsigned int>> adj2(ingfa->n_seg * 2);
  std::vector<std::vector<unsigned int>> tags(5);
  tags[0].resize(ingfa->n_seg * 2);
  tags[1].resize(ingfa->n_seg * 2);
  tags[2].resize(ingfa->n_seg * 2);
  tags[3].resize(ingfa->n_seg * 2);
  tags[4].resize(ingfa->n_seg);

  std::vector<unsigned int> labels_map(ingfa->n_seg * 2);
  int insn = 0;
  unsigned int max_l = 0;
  unsigned long int avg_l = 0;
  // std::cout << "nodes: " << ingfa->n_seg << "\n";
  for (unsigned int i = 0; i < ingfa->n_seg; ++i) {

    auto l = ingfa->seg[i].len;
    if (l > max_l) {
      max_l = l;
    }
    avg_l += l;
    std::string seq_t = ingfa->seg[i].seq;
    // labels.push_back(std::make_pair(seq_t, i));
    labels[i] = std::make_pair(seq_t, i);
    labels[i + ingfa->n_seg] = std::make_pair(rc(seq_t), i + ingfa->n_seg);
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
      adj[i + ingfa->n_seg].push_back(
          gfa_name2id(ingfa, ingfa->seg[inca_vid[j].w >> 1].name));
    }

    vid = segid << 1;
    gfa_arc_t *outa_vid = gfa_arc_a(ingfa, vid);
    uint32_t n_outa_vid = gfa_arc_n(ingfa, vid);
    for (int j = 0; j < n_outa_vid; ++j) {
      adj2[i].push_back(
          gfa_name2id(ingfa, ingfa->seg[outa_vid[j].w >> 1].name));
      adj2[i + ingfa->n_seg].push_back(
          gfa_name2id(ingfa, ingfa->seg[outa_vid[j].w >> 1].name));
    }
    labels_map[i] = std::stoi(segname);
    labels_map[i + ingfa->n_seg] = std::stoi(segname);
  }

  std::vector<std::vector<unsigned int>> adj_f(ingfa->n_seg * 2 + 1);
  std::vector<std::vector<unsigned int>> adj_f2(ingfa->n_seg * 2);
  std::vector<unsigned int> labels_map_f(ingfa->n_seg * 2);

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
    // if (l.second >= ingfa->n_seg)
    //   continue;
    map_l[l.second] = j;
    j++;
  }

  for (unsigned int i = 0; i < adj.size(); i++) {
    for (unsigned int j = 0; j < adj[i].size(); j++) {
      adj_f[map_l[i]].push_back(map_l[adj[i][j]]);
    }
  }

  adj_f[ingfa->n_seg * 2].push_back((unsigned int)max_l);
  adj_f[ingfa->n_seg * 2].push_back((unsigned int)avg_l / ingfa->n_seg);
  for (unsigned int i = 0; i < adj2.size(); i++) {
    for (unsigned int j = 0; j < adj2[i].size(); j++) {
      adj_f2[map_l[i]].push_back(map_l[adj2[i][j]]);
    }
  }

  for (unsigned int i = 0; i < labels_map.size(); i++) {
    labels_map_f[map_l[i]] = labels_map[i];
  }

  // std::cout << "---------.\n";
  // for (auto ll : labels_map_f) {
  //   std::cout << ll << "\n";
  // }
  // std::cout << "---------.\n";
  // for (auto ll : map_l) {
  //   std::cout << ll << "\n";
  // }
  int i;
  for (auto ll : labels) {
    auto seq = ll.first.c_str();
    auto l = ll.first.size();
    s = (uint8_t *)seq;

    //  change encoding
    for (int j = 0; j < l; ++j)
      s[j] = fm6_i(s[j]);

    kputsn((char *)s, l + 1, &buf);

    // // Add reverse to buffer
    // for (i = 0; i < (l >> 1); ++i) {
    //   int tmp = s[l - 1 - i];
    //   tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
    //   s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    //   s[i] = tmp;
    // }
    // if (l & 1)
    //   s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    // kputsn((char *)s, l + 1, &buf);

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
  uintmat_dump(adj_f2, graph_out2.c_str());
  uintvec_dump(labels_map_f, labels_out.c_str());
  rlc_dump(rlc, rope_out.c_str());

  // rlc_print_bwt(rlc);

  std::fprintf(stderr,
               "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
               __func__, realtime() - t_start, cputime());
  rlc_destroy(rlc);
  // int c = 0;
  // for (auto &l : labels) {

  //   std::cout << c << " " << l.first << " " << l.second << "\n";
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
  //
  std::unordered_map<int, unsigned int> im;
  for (unsigned int i = 0; i < labels.size(); ++i) {
    im[labels[i].second] = i;
  }
  int rcc = 0;
  for (int off = 0; off < labels.size(); off++) {
    if (off == 0) {
      tags[0][off] = ingfa->n_seg * 2 - 1;
      tags[1][off] = off;
    } else {
      tags[0][off] = off - 1;
      tags[1][off] = off;
    }
    tags[2][off] = labels[off].first.size();
    tags[3][off] = (labels[off].second < ingfa->n_seg) ? 0 : 1;
    if (tags[3][off] == 0) {
      auto it = im.find(labels[off].second + ingfa->n_seg);
      if (it != im.end()) {
        auto pos = it->second;
        tags[4][rcc] = pos;
      }
      // for (int off2 = 0; off2 < labels.size(); off2++) {
      //   std::cout << off << " vs " << off2 << "\n";
      //   if (labels[off2].second == labels[off].second + ingfa->n_seg) {
      //     tags[4][rcc] = off2;
      //     break;
      //   }
      // }
      rcc++;
    }
  }

  // for (int i = 0; i < tags[0].size(); i++) {
  //   std::cout << tags[0][i] << "\t" << tags[1][i] << "\t" << tags[2][i] <<
  //   "\t"
  //             << tags[3][i] << "\n";
  // }
  // for (int i = 0; i < tags[4].size(); i++) {
  //   std::cout << tags[4][i] << "\n";
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

  uintmat_dump(tags, tag_out.c_str());

  // rld_t *index = rld_restore(rope_out.c_str());

  // rldintv_t sai;
  // for (int c = 0; c < 6; ++c) {
  //   fm6_set_intv(index, c, sai);
  //   printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  // }
  // rldintv_t osai[6];
  // i = 0;
  // unsigned int t = 0;
  // for (unsigned int j = 0; j < labels.size(); ++j) {
  //   // auto l = ingfa->seg[j].len + 2;
  //   // std::string seq_tt = ingfa->seg[j].seq;
  //   // std::string seq_t = "$" + seq_tt + "$";
  //   auto l = labels[j].first.size() + 1;
  //   std::string seq_tt = "$" + labels[j].first;
  //   std::string seq_t = seq_tt;
  //   std::cout << seq_t << "\n";
  //   auto q = (uint8_t *)seq_t.c_str();

  //   // change encoding
  //   for (i = 0; i < l; ++i)
  //     q[i] = q[i] < 128 ? seq_nt6_table[q[i]] : 5;
  //   q[0] = 0;
  //   // ss[l - 1] = 0;
  //   for (i = 0; i < l; ++i)
  //     printf("%d", q[i]);
  //   printf("\n");

  //   int i, c, ret;
  //   rldintv_t ik, ok[6];

  //   int x = 0;
  //   fm6_set_intv(index, q[x], ik);
  //   // int s = (j == 0) ? labels.size() - 1 : j - 1;
  //   // ik.x[0] = s;
  //   // ik.x[1] = s;
  //   // ik.x[2] = 1;
  //   printf("%d: [%ld, %ld, %ld]\n", q[x], ik.x[0], ik.x[1], ik.x[2]);

  //   ik.info = x + 1;
  //   for (i = x + 1; i < l; ++i) { // forward search
  //     c = fm6_comp(q[i]);
  //     rld_extend(index, &ik, ok, 0);
  //     printf("%d: [%ld, %ld, %ld]\n", c, ok[c].x[0], ok[c].x[1], ok[c].x[2]);
  //     if (ok[c].x[2] == 0)
  //       break; // cannot be extended
  //     ik = ok[c];

  //     ik.info = i + 1;
  //   }

  //   printf("----------\n");
  //   unsigned int t = 0;
  //   for (unsigned int j = 0; j < labels.size(); ++j) {
  //     if (labels[j].second >= ingfa->n_seg)
  //       continue;
  //     // auto l = ingfa->seg[j].len + 2;
  //     // std::string seq_tt = ingfa->seg[j].seq;
  //     // std::string seq_t = "$" + seq_tt + "$";
  //     auto l = labels[j].first.size() + 1;
  //     std::string seq_tt = "$" + labels[j].first;
  //     std::string seq_t = seq_tt;
  //     std::cout << seq_t << "\n";
  //     auto q = (uint8_t *)seq_t.c_str();

  //     // change encoding
  //     for (i = 0; i < l; ++i)
  //       q[i] = q[i] < 128 ? seq_nt6_table[q[i]] : 5;
  //     q[0] = 0;
  //     // ss[l - 1] = 0;
  //     for (i = 0; i < l; ++i)
  //       printf("%d", q[i]);
  //     printf("\n");
  //     int i, c, ret;
  //     rldintv_t ik, ok[6];

  //     int x = 0;
  //     fm6_set_intv(index, q[x], ik);
  //     // int s = (j == 0) ? labels.size() - 1 : j - 1;
  //     ik.x[0] = t;
  //     ik.x[1] = tags[4][t];
  //     ik.x[2] = 1;
  //     printf("%d: [%ld, %ld, %ld]\n", q[x], ik.x[0], ik.x[1], ik.x[2]);

  //     ik.info = x + 1;
  //     for (i = x + 1; i < l; ++i) { // forward search
  //       c = fm6_comp(q[i]);
  //       rld_extend(index, &ik, ok, 0);
  //       printf("%d: [%ld, %ld, %ld]\n", c, ok[c].x[0], ok[c].x[1],
  //       ok[c].x[2]); if (ok[c].x[2] == 0)
  //         break; // cannot be extended
  //       ik = ok[c];

  //       ik.info = i + 1;
  //     }
  //     t++;
  //   }
  //   // for (int k = 0; k < labels.size(); k++) {
  //   //   for (int s = 0; s < labels.size(); s++) {

  //   //     x = 0;
  //   //     fm6_set_intv(index, q[x], ik);

  //   //     ik.x[0] = k;
  //   //     ik.x[1] = s;
  //   //     ik.x[2] = 1;
  //   //     // printf("%d: [%ld, %ld, %ld]\n", q[x], ik.x[0], ik.x[1],
  //   ik.x[2]);

  //   //     ik.info = x + 1;
  //   //     for (i = x + 1; i < l; ++i) { // forward search
  //   //       c = fm6_comp(q[i]);
  //   //       rld_extend(index, &ik, ok, 0);

  //   //       if (ok[c].x[2] == 0)
  //   //         break; // cannot be extended
  //   //       ik = ok[c];

  //   //       ik.info = i + 1;
  //   //     }
  //   //     if (ok[c].x[2] != 0) {
  //   //       std::cout << "testing " << k << " and " << s << " for " <<
  //   seq_tt
  //   //                 << " at " << j << "\n";
  //   //       printf("%d: [%ld, %ld, %ld]\n", c, ok[c].x[0], ok[c].x[1],
  //   //              ok[c].x[2]);
  //   //       printf("........\n");
  //   //     }
  //   //   }
  //   // }
  //   printf("-----------------------------\n");
  //   // init search
  //  Forward search:
  // int end = 0;
  // fm6_set_intv(index, ss[end], sai);
  // sai.info = end + 1;
  // printf("%d: [%ld, %ld, %ld] start\n", ss[end], sai.x[0], sai.x[1],
  //        sai.x[2]);

  // // fm6_set_intv(index, ss[end], sai);
  // // printf("%d: [%ld, %ld, %ld]\n", ss[end], sai.x[0], sai.x[1],
  // sai.x[2]); while (sai.x[2] != 0 && end < l) {
  //   end++;
  //   rldintv_t osai[6];
  //   // rldintv_t sai_f;
  //   // sai_f.x[0] = sai.x[1];
  //   // sai_f.x[1] = sai.x[0];
  //   // sai_f.x[2] = sai.x[2];
  //   // rld_extend(index, &sai_f, osai, 1);
  //   // // sai = osai[ss[end] >= 1 && ss[end] <= 4 ? 5 - ss[end] :
  //   ss[end]];
  //   // sai_f = osai[fm6_comp(ss[end])];
  //   // sai.x[0] = sai_f.x[1];
  //   // sai.x[1] = sai_f.x[0];
  //   // sai.x[2] = sai_f.x[2];
  //   // printf("%d: [%ld, %ld, %ld]\n", fm6_comp(ss[end]), sai.x[0],
  //   sai.x[1],
  //   //        sai.x[2]);
  //   rld_extend(index, &sai, osai, 0);
  //   sai = osai[fm6_comp(ss[end])];
  //   sai.info = end + 1;
  //   // sai = osai[ss[end] >= 1 && ss[end] <= 4 ? 5 - ss[end] : ss[end]];
  //   printf("%d: [%ld, %ld, %ld]\n", fm6_comp(ss[end]), sai.x[0],
  //   sai.x[1],
  //          sai.x[2]);
  //   if (sai.x[2] == 0) {
  //     for (int c = 0; c < 6; ++c) {
  //       sai = osai[fm6_comp(c)];
  //       printf("\t%d: [%ld, %ld, %ld]\n", fm6_comp(c), sai.x[0],
  //       sai.x[1],
  //              sai.x[2]);
  //     }
  //   }
  // }

  // rld_extend(index, &sai, osai, 0);
  // for (int c = 0; c < 6; ++c) {
  //   sai = osai[fm6_comp(c)];
  //   printf("\t%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  // }

  // init search
  //
  // rldintv_t sai;
  // for (int c = 0; c < 6; ++c) {
  //   fm6_set_intv(index, c, sai);
  //   printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  // }
  // i = l - 1;
  // fm6_set_intv(index, ss[i], sai);
  // printf("%d: [%ld, %ld, %ld]\n", ss[i], sai.x[0], sai.x[1], sai.x[2]);
  // --i;

  // // backward extensions
  // for (; i >= 0; --i) {
  //   rld_extend(index, &sai, osai, 1);

  //   sai = osai[ss[i]];
  //   printf("%d: [%ld, %ld, %ld]\n", ss[i], sai.x[0], sai.x[1], sai.x[2]);
  //   if (sai.x[2] <= 0) {
  //     // errors += 1;
  //     break;
  //   }
  // }
  // rld_extend(index, &sai, osai, 1);
  // for (int c = 0; c < 6; ++c) {
  //   sai = osai[c];
  //   printf("\t%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  // }
  //
  // rld_extend(index, &sai, osai, 1);
  // sai = osai[0];
  // printf("[%ld, %ld, %ld]\n", sai.x[0], sai.x[1], sai.x[2]);
  //}

  // printf("-----------\n");

  // printf("\n");
  // for (unsigned int c = 0; c < labels.size(); ++c) {
  //   if (tags[3][c] == 1) {
  //     continue;
  //   }
  //   if (c == 0) {
  //     c = tags[0].size() - 1;
  //   }
  //   std::cout << labels[c].first << "\n";
  //   sai.x[0] = c;
  //   sai.x[1] = c;
  //   sai.x[2] = 1;
  //   printf("%d\n", c);
  //   printf("%d: [%ld, %ld, %ld] -> \n\n", c, sai.x[0], sai.x[1], sai.x[2]);
  //   rld_extend(index, &sai, osai, 0);
  //   for (int cc = 0; cc < 6; ++cc) {
  //     sai = osai[cc];
  //     printf("%d: [%ld, %ld, %ld]\n", cc, sai.x[0], sai.x[1], sai.x[2]);
  //   }
  //   printf("\n");
  // }

  gfa_destroy(ingfa);
  free(buf.s);
}

#endif // GRAPHINDEX_BWT_ROPE_H
