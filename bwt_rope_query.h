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

unsigned long long int interval_size = 0;
unsigned long long int interval_size_avg = 1;

std::vector<std::vector<unsigned int>> interval_spectrum;

void ext(const rld_t *index, const uint8_t *s, int i, unsigned int l,
         rldintv_t &sai, std::vector<std::vector<unsigned int>> &tags,
         std::vector<std::vector<unsigned int>> &adj,
         std::vector<unsigned int> labels_map, unsigned int curr_node,
         char *read_name, unsigned int r_c, std::vector<node_sai> int_s) {

  interval_spectrum.push_back({});
  unsigned int tollerance = 0;

  uint8_t symb = i > 0 ? s[i] : 5;
  i--;

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
            unsigned int end_node = tags[0][rld_rank11(index, sai.x[0], 0)];
            if (end_node == tags[0].size() - 1) {
              std::fprintf(stdout, "@%s\t%d\n", read_name, labels_map[0]);
            } else {
              std::fprintf(stdout, "@%s\t%d\n", read_name,
                           labels_map[end_node + 1]);
            }
          } else {
            unsigned int node_f = findnode(index, k, tags[0].size(), tags);
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

  std::vector<std::vector<unsigned int>> adj = uintmat_load(g_file.c_str());
  std::vector<unsigned int> labels_map = uintvec_load(l_file.c_str());

  std::vector<std::vector<std::vector<uint64_t>>> c_int =
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
      sai.info = 0;
      sai.x[0] = v[0];
      sai.x[1] = v[0];
      sai.x[2] = v[1];
      tv.push_back({sai, v[2], {v[2]}});
    }
    cache[k] = tv;
    k++;
  }

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

    l -= (ls - 1);
    i = l - 1;

    ext(index, s, i, l, sai, tags, adj, labels_map, tags[0].size(), ks->name.s,
        r_c, int_v);

    r_c++;
  }

  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);
}

#endif // GRAPHINDEX_BWT_ROPE_QUERY_H
