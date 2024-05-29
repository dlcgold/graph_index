#include <getopt.h>
#include <iostream>
#include <string>

#include "bwt_rope.h"

void printHelp() {
  std::cout << "Usage: gindex [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
}

int main(int argc, char **argv) {

  if (argc == 1) {
    printHelp();
    exit(EXIT_SUCCESS);
  }
  std::string gfa_file = "";
  std::string out_prefix = "";
  int threads = 1;
  while (true) {
    static struct option long_options[] = {
        {"input", required_argument, nullptr, 'i'},
        {"output", required_argument, nullptr, 'o'},
        {"threads", required_argument, nullptr, 't'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    int c = getopt_long(argc, argv, "i:o:t:", long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 'i':
      gfa_file = optarg;
      break;
    case 'o':
      out_prefix = optarg;
      break;
    case 't':
      threads = atoi(optarg);
      break;
    case 'h':
      printHelp();
      exit(EXIT_SUCCESS);
    default:
      printHelp();
      exit(EXIT_FAILURE);
    }
  }
  build_bwt_rope(gfa_file.c_str(), out_prefix, threads);
  return 0;
}
