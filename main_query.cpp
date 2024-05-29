#include <getopt.h>
#include <iostream>
#include <string>

#include "bwt_rope_query.h"

void printHelp() {
  std::cout << "Usage: gindex [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
}

int main(int argc, char **argv) {

  if (argc == 1) {
    printHelp();
    exit(EXIT_SUCCESS);
  }
  std::string index_file = "";
  std::string query_file = "";
  int threads = 1;
  while (true) {
    static struct option long_options[] = {
        {"input_prefix", required_argument, nullptr, 'i'},
        {"query", required_argument, nullptr, 'q'},
        {"threads", required_argument, nullptr, 't'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    int c = getopt_long(argc, argv, "i:q:t:", long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 'i':
      index_file = optarg;
      break;
    case 'q':
      query_file = optarg;
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
  query_bwt_rope(index_file, query_file.c_str(), threads);
  return 0;
}
