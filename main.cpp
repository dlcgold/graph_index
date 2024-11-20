#include "argparse/argparse.hpp"
#include "gindex_index.h"
#include "gindex_query.h"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("gindex");

  argparse::ArgumentParser index_command("index");

  std::string gfa_file = "";
  std::string out_prefix = "";
  int threads = 1;
  int cache = 0;

  index_command.add_description(
      "Index a given GFA file to perform exact pattern-matching");
  index_command.add_argument("-i", "--input")
      .help("input graph in GFA format")
      .store_into(gfa_file);
  index_command.add_argument("-t", "--thread")
      .default_value(1)
      .help("number of threads (default 1)")
      .store_into(threads);
  index_command.add_argument("-c", "--cache")
      .default_value(0)
      .help("lenght of kmer to preprocess (default 1)")
      .store_into(cache);
  index_command.add_argument("-o", "--output")
      .help("output prefix")
      .store_into(out_prefix);

  program.add_subparser(index_command);

  argparse::ArgumentParser query_command("query");

  std::string index_file = "";
  std::string query_file = "";
  int threadsq = 1;
  bool full = false;
  bool nocache = false;
  query_command.add_description(
      "Query a given index, queries in FASTA/FASTQ format");
  query_command.add_argument("-i", "--input")
      .help("input prefix for the index")
      .store_into(index_file);
  query_command.add_argument("-q", "--query")
      .help("query FASTA/FASTQ file")
      .store_into(query_file);
  query_command.add_argument("-t", "--thread")
      .default_value(1)
      .help("number of threads (default 1)")
      .store_into(threadsq);
  query_command.add_argument("-f", "--full")
      .default_value(false)
      .implicit_value(true)
      .help("compute all nodes path for a match (warning: cache will not be "
            "used, very slow)")
      .store_into(full);
  query_command.add_argument("-n", "--nocache")
      .default_value(false)
      .implicit_value(true)
      .help("ignore cache")
      .store_into(nocache);
  program.add_subparser(query_command);

  if (argc == 1) {
    std::cerr << program;
    return 1;
  }
  try {
    program.parse_args(argc, argv);
  } catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  if (program.is_subcommand_used("index")) {
    gindex_index(gfa_file.c_str(), out_prefix, threads, cache);
  } else {
    if (full == true) {
      gindex_query_path(index_file, query_file.c_str(), threadsq);
    } else {
      gindex_query(index_file, query_file.c_str(), threadsq, nocache);
    }
  }
  //
  return 0;
}
