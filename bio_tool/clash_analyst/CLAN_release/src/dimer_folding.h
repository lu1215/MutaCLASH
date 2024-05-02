#ifndef _DIMER_FOLDING_H_
#define _DIMER_FOLDING_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <list>
#include <queue>
#include <stack>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <unordered_map>
#include <ctype.h>

#include "free_energy.h"

struct RNAStackType {
  int begin, end;
  float energy;
  
  void operator=(const RNAStackType &s) {
    begin = s.begin; end = s.end;
    energy = s.energy;
    return;
  }
};

struct DuplexInfo {
  float coverage;
  float dimer_energy;
  std::vector<int> dimer_pair;
  RNAStackType stack_info;
  
  void operator=(const DuplexInfo &d) {
    coverage = d.coverage;
    dimer_energy = d.dimer_energy;
    dimer_pair = d.dimer_pair;
    stack_info = d.stack_info;
    return;
  }
};


class DimerFolding  {
 public:
  
  DimerFolding()  {}
  ~DimerFolding() {}
  
  void FoldDimerHelixEnergy(
    std::string &seq,     // the sequence to be folded 
    const int terminal,   // the connection terminal of the two RNA strands, no intramolecular pair is predicted
    DuplexInfo &dimer_info,
    const int max_bulge_size = 10
  );
  
  
 private:
  bool FormatRNASeq(std::string &seq, std::string &formatted_seq)  {
    for(int i = 0; i < seq.length(); ++ i) {
      switch (toupper(seq[i]))  {
        case 'A':
          formatted_seq[i] = 'A';
          continue;
        case 'C':
          formatted_seq[i] = 'C';
          continue;
        case 'G':
          formatted_seq[i] = 'G';
          continue;
        case 'T':
          formatted_seq[i] = 'U';
          continue; 
        case 'U':
          formatted_seq[i] = 'U';
          continue; 
        default:
          std::cerr << "Warning: CLAN: DimerFolding::FormatRNASeq: unrecognized character for RNA sequence: " << seq[i] << "." << std::endl;
          return false;
      }
    }
    return true;
  }
};


#endif
