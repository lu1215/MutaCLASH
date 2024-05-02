#ifndef _OUTPUT_PRINTER_H_
#define _OUTPUT_PRINTER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <ctime>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <unistd.h>
#include <iomanip>

#include "bwt.h"
#include "bwt_search.h"
#include "loader.h"
#include "bio_alphabet.h"
#include "minimizer_sort.h"
#include "bwt_frag_merger.h"
#include "score.h"
#include "dimer_folding.h"
#include "read_merge.h"

struct DuplexType {
  int island_ID1, island_ID2;
  bool strand_1, strand_2;
  DuplexInfo info;
  
  void operator=(const DuplexType &d) {
    island_ID1 = d.island_ID1; island_ID2 = d.island_ID2;
    strand_1 = d.strand_1; strand_2 = d.strand_2;
    info = d.info;
    return;
  }
};


class OutputPrinter {
 public:
  void OutputPrinterCLAN(
    char const* const* ref_header, char const* const* ref_seq,
    char const* const* ref_header_aux, char const* const* ref_seq_aux,
    char const* const* header, char const* const* seq,          // the header and the sequence of the read set
    const int begin, const int size, const int upper,    // the begin, the number of sequences, and the upper bound of sequences to map in this batch
    std::vector< std::vector<FragType> > &result,      // the output 
    std::ofstream &out_fh               // the output file handle
  );
  
  void ParseCIGAR(const char* cigar, int &align_len, int &num_gap);
  
  void OutputPrinterBLAST(
    char const* const* ref_header, char const* const* ref_seq,
    char const* const* header, char const* const* seq,          // the header and the sequence of the read set
    const int begin, const int size, const int upper,   // the begin and the number of sequences to map in this batch
    std::vector< std::vector<FragType> > &result,       // the output 
    const std::string &db_name, const long long int db_size,   // the name of the database
    std::ofstream &out_fh               // the output file handle
  );
  
  void OutputPrinterSAM(
    char const* const* ref_header, char const* const* ref_seq,
    char const* const* header, char const* const* seq,
    const int begin, const int size, const int upper,
    std::vector< std::vector<FragType> > &result,
    const long long int db_size,
    std::ofstream &out_fh
  );
  
  void OutputEncoded(
    char const* const* ref_seq, const int num_refs,
    const int begin, const int size, const int upper,
    std::vector< std::vector<FragType> > &result, 
    std::ofstream &out_fh  
  );
  
  void ReadEncoded(
    const int max_entry, const std::string &file_name,
    std::vector< std::vector<FragType> > &result
  );
  
  void OutputAnnotationTable(
    const bool print_simplex, char const* const* ref_header, 
    std::vector<GPosType> &islands, 
    std::vector<DuplexType> &duplex_annot,
    std::ofstream &out_fh
  );
  
  void OutputPredictedStructure();
};

#endif
