//
// Created by dlcgold on 09/05/24.
//

#ifndef GRAPH_INDEX_UTILS_H
#define GRAPH_INDEX_UTILS_H

// #include "common.h"
#include <arpa/inet.h>
#include <cstdio>
#include <fstream>
#include <stdint.h>
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>

// #include "bwt_rope.h"

// #include "kseq.h"

char map[5] = {'$', 'A', 'C', 'G', 'T'};

/* #define max(a, b) \ */
/*   ({ \ */
/*     __typeof__(a) _a = (a); \ */
/*     __typeof__(b) _b = (b); \ */
/*     _a > _b ? _a : _b; \ */
/*   }) */

#define fm6(a) ((a) < 128 ? seq_nt6_table[(a)] : 5)
#define fm6_i(a) ((a) < 128 ? seq_nt6_table_i[(a)] : 5)

#define fm6_comp(a) ((a) >= 1 && (a) <= 4 ? 5 - (a) : (a))

#define fm6_set_intv(e, c, ik)                                                 \
  ((ik).x[0] = (e)->cnt[(int)(c)],                                             \
   (ik).x[2] = (e)->cnt[(int)(c) + 1] - (e)->cnt[(int)(c)],                    \
   (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

// LD: why do we need this?
static const unsigned int seq_nt6_table_i[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
    5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

void uintvec_dump(std::vector<unsigned int> vec, const char *file) {
  std::ofstream ofs;
  ofs.open(file, std::ios::out | std::ios::binary);

  uint32_t sz = htonl(vec.size());
  ofs.write((const char *)&sz, sizeof(uint32_t));
  for (uint32_t i = 0, end_i = vec.size(); i < end_i; ++i) {
    int32_t val = htonl(vec[i]);
    ofs.write((const char *)&val, sizeof(int32_t));
  }

  ofs.close();
}

std::vector<unsigned int> uintvec_load(const char *file) {

  std::ifstream ifs;
  ifs.open(file, std::ios::in | std::ios::binary);

  uint32_t sz = 0;
  ifs.read((char *)&sz, sizeof(uint32_t));
  sz = ntohl(sz);
  std::vector<unsigned int> vec(sz);

  for (uint32_t i = 0; i < sz; ++i) {
    int32_t val = 0;
    ifs.read((char *)&val, sizeof(int32_t));
    val = ntohl(val);
    vec[i] = val;
  }
  return vec;
}

void uintmat_dump(std::vector<std::vector<unsigned int>> mat,
                  const char *file) {
  std::ofstream ofs;
  ofs.open(file, std::ios::out | std::ios::binary);

  uint32_t sz = htonl(mat.size());
  ofs.write((const char *)&sz, sizeof(uint32_t));
  for (auto &vec : mat) {
    uint32_t szv = htonl(vec.size());
    ofs.write((const char *)&szv, sizeof(uint32_t));
    for (uint32_t i = 0, end_i = vec.size(); i < end_i; ++i) {
      int32_t val = htonl(vec[i]);
      ofs.write((const char *)&val, sizeof(int32_t));
    }
  }

  ofs.close();
}

std::vector<std::vector<unsigned int>> uintmat_load(const char *file) {
  std::ifstream ifs;
  ifs.open(file, std::ios::in | std::ios::binary);

  uint32_t sz = 0;
  ifs.read((char *)&sz, sizeof(uint32_t));
  sz = ntohl(sz);
  std::vector<std::vector<unsigned int>> mat(sz);
  for (uint32_t i = 0; i < sz; ++i) {

    uint32_t szv = 0;
    ifs.read((char *)&szv, sizeof(uint32_t));
    szv = ntohl(szv);
    std::vector<unsigned int> vec(szv);
    for (uint32_t j = 0; j < szv; ++j) {
      int32_t val = 0;
      ifs.read((char *)&val, sizeof(int32_t));
      val = ntohl(val);
      vec[j] = val;
    }
    mat[i] = vec;
  }
  return mat;
}

void uintmat3_dump(std::vector<std::vector<std::vector<unsigned int>>> mat,
                   const char *file) {
  std::ofstream ofs;
  ofs.open(file, std::ios::out | std::ios::binary);
  uint32_t sz = htonl(mat.size());
  ofs.write((const char *)&sz, sizeof(uint32_t));
  for (auto vec2 : mat) {
    uint32_t szv2 = htonl(vec2.size());
    ofs.write((const char *)&szv2, sizeof(uint32_t));
    for (auto &vec : vec2) {
      uint32_t szv = htonl(vec.size());
      ofs.write((const char *)&szv, sizeof(uint32_t));
      for (uint32_t i = 0, end_i = vec.size(); i < end_i; ++i) {
        int32_t val = htonl(vec[i]);
        // std::cout << vec[i] << " ";
        ofs.write((const char *)&val, sizeof(int32_t));
      }
    }
  }

  ofs.close();
}

std::vector<std::vector<std::vector<unsigned int>>>
uintmat3_load(const char *file) {
  std::ifstream ifs;
  ifs.open(file, std::ios::in | std::ios::binary);

  uint32_t sz = 0;
  ifs.read((char *)&sz, sizeof(uint32_t));
  sz = ntohl(sz);
  std::vector<std::vector<std::vector<unsigned int>>> mat(sz);
  for (uint32_t i = 0; i < sz; ++i) {

    uint32_t szv2 = 0;
    ifs.read((char *)&szv2, sizeof(uint32_t));
    szv2 = ntohl(szv2);
    std::vector<std::vector<unsigned int>> vec2(szv2);
    for (uint32_t k = 0; k < szv2; ++k) {
      uint32_t szv = 0;
      ifs.read((char *)&szv, sizeof(uint32_t));
      szv = ntohl(szv);
      std::vector<unsigned int> vec(szv);
      for (uint32_t j = 0; j < szv; ++j) {
        int32_t val = 0;
        ifs.read((char *)&val, sizeof(int32_t));
        val = ntohl(val);

        vec[j] = val;
      }
      vec2[k] = vec;
    }
    mat[i] = vec2;
  }
  return mat;
}

/* static inline uint kputsn(const char *p, uint64_t l, kstring_t *s) { */
/*   if (s->l + l + 1 >= s->m) { */
/*     char *tmp; */
/*     s->m = s->l + l + 2; */
/*     kroundup32(s->m); */
/*     if ((tmp = (char *)realloc(s->s, s->m))) */
/*       s->s = tmp; */
/*     else */
/*       return EOF; */
/*   } */
/*   memcpy(s->s + s->l, p, l); */
/*   s->l += l; */
/*   s->s[s->l] = 0; */
/*   return l; */
/* } */

/* static inline double cputime() { */
/*   struct rusage r; */
/*   getrusage(RUSAGE_SELF, &r); */
/*   return (double)r.ru_utime.tv_sec + (double)r.ru_stime.tv_sec + */
/*          1e-6 * (double)(r.ru_utime.tv_usec + r.ru_stime.tv_usec); */
/* } */

/* static inline double realtime() { */
/*   struct timeval tp; */
/*   struct timezone tzp; */
/*   gettimeofday(&tp, &tzp); */
/*   return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6; */
/* } */

#endif // GRAPH_INDEX_UTILS_H
