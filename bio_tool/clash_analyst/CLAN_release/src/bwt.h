#ifndef _BWT_H_
#define _BWT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <queue>
#include <stack>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <tuple>
#include <unordered_map>

#include "divsufsort64.h"
#include "bio_alphabet.h"
#include "concatenator.h"

// the gap between two recorded occurrences in the FM-index
#define BWT_FM_GAP 64

/* Definition from 64-bit libdivsufsort */
// typedef uint8_t sauchar_t;
typedef sauchar_t BWTCHAR;
// typedef int64_t saidx_t;
typedef saidx64_t BWTIDX;
// typedef int32_t saint_t
typedef saint_t BWTINT;

class BWT {
 public:
  explicit BWT() : 
    build_success_(false),
    ref_success_(false),
    bwt_(NULL),
    position_(NULL),
    acc_freq_(NULL)  
  {}
  ~BWT()  {}
  void Purge(void)  {
    if(build_success_)  {
      delete [] bwt_;
      delete [] str_;
      delete [] position_;
      delete [] acc_freq_;
      for(auto it = fm_index_.begin(); it != fm_index_.end(); ++ it)
        delete [] it->second;
    }
    if(ref_success_)  {
      delete [] acc_len_;
      for(int i = 0; i < num_seqs_; ++ i) {
        delete [] header_[i];
      }
      delete [] header_;
    }
    return;
  }
  void PrintBWT(void) {std::cout << bwt_ << std::endl;}
  
  // return the size
  inline BWTIDX GetSize() {return size_;}
  
  // constructing the BWT
  void Construct(const BioAlphabet &alphabet, const char *text);
  void ConstructFromIndex(
      const BioAlphabet &alphabet, const char *text, 
      const std::string &file_prefix
  );
  // constructing the BWT without read ID anf position infomation
  void ConstructNoPos(const BioAlphabet &alphabet, const char *text);
  // construct the reference sequence set information for location mapping
  void ConstructNoPosFromIndex(
      const BioAlphabet &alphabet, const char *text, 
      const std::string &file_prefix
  );
  // retrieve the reference position given a corresponding BWT position
  // returns the absolute location
  BWTIDX BWTToStrPos(const BWTIDX bwt_pos);
  void StrPosToRefPos(const BWTIDX str_pos, std::string &ref_header, BWTIDX &ref_pos);
  void StrPosToRefPos(const BWTIDX str_pos, int &ref_header_idx, BWTIDX &ref_pos);
  
  void ConstructLocationInfo(const int n, char **header, char **seq);
  // update the current "range" while the preceeding character is "c"
  std::pair<BWTIDX, BWTIDX> UpdateRange(const char c, const std::pair<BWTIDX, BWTIDX> &range);
  void UpdateRange(
      const char c, const BWTIDX old_begin, const BWTIDX old_end, 
      BWTIDX &new_begin, BWTIDX &new_end
  );
  // counting the occurrence of character "c" at positions in the BWT that are less than "pos" 
  inline BWTIDX CountOccurrence(const char c, const BWTIDX pos) {
    BWTIDX occ = fm_index_[c][(BWTIDX) pos / BWT_FM_GAP];
    for(BWTIDX i = pos - (pos % BWT_FM_GAP); i < pos; ++ i) occ += bwt_[i] == c ? 1 : 0;
    return occ;
  }
  // counting the occurrence of chars that are lexicographically less than "c"
  // at positions in the BWT that are less than "pos"
  BWTIDX CountLexicoLess(const char c, const BWTIDX pos);
  
  // dumping the BWT index to disk
  void WriteBWTIndex(const std::string &file_name);
  // dumping the positions to disk
  void WritePosition(const std::string &file_name);
  // dumping the FM index to disk
  void WriteFMIndex(const std::string &file_name);
  // loading the BWT index from disk
  void LoadBWTIndex(const std::string &file_name);
  // loading the position index from disk
  void LoadPosition(const std::string &file_name);
  // loading the FM index from the disk
  void LoadFMIndex(const std::string &file_name);
  
  void WriteIndex(const std::string &file_prefix);
  void WriteIndexNoPos(const std::string &file_prefix);
  void PrintSubstr(const BWTIDX begin, const BWTIDX len)    {
    for(BWTIDX i = 0; i < len; ++ i)    std::cerr << str_[begin + i];
    std::cerr << std::endl;
    return;
  }
  
  friend class BWTSearch;

 protected:
  bool build_success_;
  bool ref_success_;
  BioAlphabet alphabet_;  // the alphebet in use
  BWTCHAR *bwt_;          // the BWT string
  BWTCHAR *str_;          // the raw string
  //BWTIDX *tmp_position_;  // the tmp holder for the positions generated
  BWTIDX *position_;      // the gapped positions
  int num_seqs_;          // the number of sequences used to build the BWT  
  BWTCHAR **header_;       // the headers of the original string
  BWTIDX *acc_len_;       // the accumulated length of the sequence set
  // the FM-index
  std::unordered_map<BWTCHAR, BWTIDX*> fm_index_;
  // the frequency of the characters     
  std::unordered_map<BWTCHAR, BWTIDX> char_freq_;
  // the accumulated frequency, which is equivalent to index in the sorted column
  BWTIDX *acc_freq_;   
  BWTIDX size_;           // the size of the text 
  
  // counting the frequency of the characters in the alphabet, excluding "delim"
  // sorts out the BWT, and also builds the FM-index
  void BuildIndex(const BWTCHAR* text);
  void BuildIndexNoPos(const BWTCHAR* text);
  
  
};

#endif
