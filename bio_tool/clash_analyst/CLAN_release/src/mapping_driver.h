#ifndef _MAPPING_DRIVER_H_
#define _MAPPING_DRIVER_H_

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
#include "parameters.h"

class MappingDriver {
 public:
  void ReadMapBatch(
    char **ref_header, char **ref_seq, int num_refs,
    char **ref_header_aux, char **ref_seq_aux, int num_refs_aux,
    char **header, char **seq, int num_reads,          // the header and the sequence of the read set
    const int begin, const int size, const int upper,  // the begin, the number of sequences, and the upper bound of sequences to map in this batch
    BWT &bwt, BWT &bwt_aux,                            // the indexed BWT for the main and auxiliary databases
    std::vector< std::vector<FragType> > &result,      // the output 
    const CLANSearchParam &param_set
  );
  
  
};

#endif


