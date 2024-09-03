#ifndef COMMON_H_
#define COMMON_H_

struct node_sai {
  rldintv_t sai;
  unsigned int curr_node;
  std::vector<unsigned int> path;

  int end() const { return sai.x[0] + sai.x[2]; }
  bool operator<(const node_sai &other) const {
    return sai.x[0] < other.sai.x[0];
  }
  bool operator==(const node_sai &other) const {
    return sai.x[0] == other.sai.x[0] && sai.x[2] == other.sai.x[2];
  }
};

#endif // COMMON_H_
