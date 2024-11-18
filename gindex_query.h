//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPHINDEX_BWT_ROPE_QUERY_H
#define GRAPHINDEX_BWT_ROPE_QUERY_H
#include "utils.h"

#include "common.h"
#include "kseq.h"
#include "rld0.h"
#include "utils.h"
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <filesystem>
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

void gindex_query_path(std::string index_pre, const char *query_file,
                       int threads) {
  double t_start = realtime();
  auto i_file = index_pre + ".bwt";
  auto t_file = index_pre + ".dollars";
  auto g_file = index_pre + ".graph";
  auto l_file = index_pre + ".labels";
  fprintf(stderr,
          "[M::%s] requested full path, cache will not be used. Expect slow "
          "execution time\n",
          __func__);

  rld_t *index = rld_restore(i_file.c_str());

  std::vector<std::vector<uint64_t>> tags = uintmat_load(t_file.c_str());

  // for (auto n : tags[0]) {
  //   std::cerr << n << " ";
  // }
  // std::cout << "\n-----------\n";
  // for (auto n : tags[1]) {
  //   std::cerr << n << " ";
  // }
  // std::cerr << "\n";

  std::vector<std::vector<uint64_t>> adj = uintmat_load(g_file.c_str());
  std::vector<uint64_t> labels_map = uintvec_load(l_file.c_str());
  // for (auto v : adj) {
  //   std::cout << cc << ": ";
  //   for (auto n : v) {
  //     std::cout << n << " ";
  //   }
  //   std::cout << "\n";
  //   cc++;
  // }
  std::fprintf(stderr, "[M::%s] index loaded in %.3f sec\n", __func__,
               realtime() - t_start);

  uint8_t *s;
  // q interval
  rldintv_t sai;
  for (int c = 0; c < 6; ++c) {
    fm6_set_intv(index, c, sai);
    // printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  }
  t_start = realtime();
  gzFile fp = gzopen(query_file, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint64_t r_c = 0;
  uint64_t al = 0;
  while ((l = kseq_read(ks)) >= 0) {
    fprintf(stderr, "[M::%s] analyzed %ld reads\r", __func__, r_c);
    al += l;
    s = (uint8_t *)ks->seq.s;
    // change encoding
    for (i = 0; i < l; ++i) {
      s[i] = fm6_i(s[i]);
    }
    // for (i = 0; i < l; ++i) {
    //   printf("%d ", s[i]);
    // }
    // std::cout << "\n";
    // rldintv_t osai[6];
    i = l - 1;
    fm6_set_intv(index, s[i], sai);
    if (i == 0) {

      for (uint64_t k = sai.x[0]; k < sai.x[0] + sai.x[2]; k++) {
        auto node_f = findnode(index, k, tags[0].size(), tags);
        // std::fprintf(stderr, "Match at node %d\n", node_f);
        std::fprintf(stderr, "%s\t%ld\n", ks->name.s, labels_map[node_f]);
      }

    } else {
      --i;
      auto s_b = sai.x[1];
      sai.x[1] = tags[0].size();
      // auto nodes = ext(index, s, i, l, sai, tags, adj, tags[0].size());
      ext_path(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(),
               ks->name.s);
      // std::fprintf(stderr, "finished");
      sai.x[1] = s_b;
    }
    r_c++;
    // std::fprintf(stderr, "\n-----------------\n");
  }
  std::fprintf(
      stderr,
      "[M::%s] all %ld of avg. length %ld queries analyzed in %.3f sec\n",
      __func__, r_c, al / r_c, realtime() - t_start);
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);
}

void gindex_query(std::string index_pre, const char *query_file, int threads) {
  double t_start = realtime();
  auto i_file = index_pre + ".bwt";
  auto t_file = index_pre + ".dollars";
  auto g_file = index_pre + ".graph";
  auto l_file = index_pre + ".labels";
  auto c_file = index_pre + ".cache";

  rld_t *index = rld_restore(i_file.c_str());

  // std::cerr << "size: " << index->l << "\n";
  std::vector<std::vector<uint64_t>> tags = uintmat_load(t_file.c_str());

  std::vector<std::vector<uint64_t>> adj = uintmat_load(g_file.c_str());
  std::vector<uint64_t> labels_map = uintvec_load(l_file.c_str());

  std::vector<std::vector<std::vector<uint64_t>>> c_int;
  std::vector<std::vector<node_sai>> cache;
  std::vector<std::string> suff;
  std::map<std::string, uint64_t> suff_map;
  std::vector<std::string> sfs;
  uint64_t ls = 0;
  std::fprintf(stderr, "[M::%s] index loaded in %.3f sec\n", __func__,
               realtime() - t_start);

  if (std::filesystem::exists(c_file)) {
    c_int = uintmat3_load(c_file.c_str());
    cache.resize(c_int.size());
    suff.resize(cache.size());
    ls = log4(cache.size());
    // std::cerr << "seeds of size " << cache.size() << " " << ls << "\n";
    t_start = realtime();
    fprintf(stderr, "[M::%s] seeds of size %ld having %ld cache intervals\n",
            __func__, ls, cache.size());
    for (uint64_t i = 0; i < cache.size(); ++i) {
      std::string suff_t(ls, ' ');
      int cur = i;

      for (uint64_t j = 0; j < ls; j++) {
        suff_t[j] = "ACGT"[cur % 4];
        cur /= 4;
      }
      // std::cout << suff_t << "\n";
      suff_map[suff_t] = i;
      sfs.push_back(suff_t);
      // std::cout << suff[i] << "\n";
    }
    //
    int k = 0;
    for (auto s : c_int) {
      std::vector<node_sai> tv;
      rldintv_t sai;
      // uint64_t curr_node;
      for (auto v : s) {
        // std::cerr << "suff: "  << sfs[k] << " ->";
        // fprintf(stderr, " cache %ld %ld %ld -> ", v[0], v[1], v[2]);
        // auto en =  get_end_nodes(index, v[0], v[1], tags, adj);
        // for(auto enn: en)  {fprintf(stderr, " %ld ", enn); if (enn==57905700)
        // exit(1);}fprintf(stderr, "\n");
        sai.info = 0;
        sai.x[0] = v[0];
        sai.x[1] = v[0];
        sai.x[2] = v[1];
        tv.push_back({sai, v[2], {v[2]}});
      }
      cache[k] = tv;
      k++;
    }
    std::fprintf(stderr, "[M::%s] cache loaded in %.3f sec\n", __func__,
                 realtime() - t_start);
  } else {
    std::fprintf(stderr, "[M::%s] cache not found\n", __func__);
  }

  /* uint64_t cc = 0;
     for (auto v : adj) {
        std::cout << cc << ": ";
        for (auto n : v) {
          std::cout << n << " ";
        }
        std::cout << "\n";
        cc++;
      }*/

  t_start = realtime();
  uint64_t al = 0;
  // auto int_tt = cache[suff_map[std::string("CTGTA")]];
  /*for(auto iii: int_tt){
 fprintf(stderr, "CTGTA %ld %ld %ld -> ", iii.sai.x[0], iii.sai.x[2],
 iii.curr_node); auto en =  get_end_nodes(index, iii.sai.x[0], iii.sai.x[2],
 tags, adj); for(auto enn: en)  {fprintf(stderr, " %ld ", enn); if
 (enn==57905700) ;}fprintf(stderr, "\n");
 }*/
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
  uint64_t r_c = 0;

  std::fflush(stdout);
  // #pragma omp parallel
  while ((l = kseq_read(ks)) >= 0) {
    // std::fprintf(stderr, "%d\r", r_c);
    fprintf(stderr, "[M::%s] analyzed %ld reads\r", __func__, r_c);
    s = (uint8_t *)ks->seq.s;
    al += l;
    std::string suf(ls, ' ');
    int sufc = 0;
    // change encoding
    for (i = 0; i < l; ++i) {
      s[i] = fm6_i(s[i]);
      if (ls != 0 && i >= l - ls) {
        suf[sufc] = map[s[i]];
        sufc++;
      }
    }

    // std::cerr << "suff: "  << suf << "\n";
    std::vector<node_sai> int_v;
    if (ls != 0) {
      int_v = cache[suff_map[suf]];
    }
    /*for(auto v: int_v){
    fprintf(stderr, "cache %ld %ld\n", v.sai.x[0], v.sai.x[2]);
    if(v.sai.x[2]==1){
    auto en =  get_end_nodes(index, v.sai.x[0], v.sai.x[2], tags, adj);
    for(auto enn: en){
    fprintf(stderr, "end for %ld %ld\n", v.sai.x[0], enn);
    }
    }
    }*/
    //   std::cout << "start size: " << int_v.size() << std::endl;
    if (ls != 0 && int_v.empty())
      continue;
    if (ls != 0)
      l -= (ls - 1);
    i = l - 1;
    fm6_set_intv(index, s[i], sai);
    ext(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(), ks->name.s,
        r_c, int_v);
    // for(auto is: interval_spectrum[0]) std::cerr << is << " ";
    r_c++;
  }
  std::fprintf(
      stderr,
      "[M::%s] all %ld of avg. length %ld queries analyzed in %.3f sec\n",
      __func__, r_c, al / r_c, realtime() - t_start);
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);
}

#endif // GRAPHINDEX_BWT_ROPE_QUERY_H
