#ifndef _BWT_FRAG_MERGER_
#define _BWT_FRAG_MERGER_

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

#include "bwt.h"
#include "bwt_search.h"

#ifndef DIFF_ARRAY_SIZE
#define DIFF_ARRAY_SIZE 4096
#endif

struct TargetType {
  int sid;
  BWTIDX begin, end;
  bool strand;        // 'true' stands for positive strand, 'false' for negative strand
  bool db;            // 'true' stands for main DB, 'false' for auxiliary DB
  char *cigar;        // cigar string for the alignment
  int num_mismatch;
  
  void operator=(const TargetType &b)  {
    sid = b.sid; begin = b.begin; end = b.end;
    strand = b.strand; db = b.db; num_mismatch = b.num_mismatch;
    cigar = new char [CIGAR_ARRAY_SIZE]; strcpy(cigar, b.cigar);
    return;
  }
};

struct FragType {
  BWTINT q_begin, q_end;
  std::vector<TargetType> targets; 
  int sol_id;
  int score;          // best score for the alignment
  bool db;            // 'true' represents the main DB; 'false' represents the auxiliary DB
  void operator=(const FragType &b)  {
    q_begin = b.q_begin; q_end = b.q_end;
    sol_id = b.sol_id; score = b.score; db = b.db;
    targets.resize(b.targets.size());
    for(int i = 0; i < b.targets.size(); ++ i) {
      targets[i] = b.targets[i];
    }
    return;
  }
};

class BWTFragMerger {
 public:
  explicit BWTFragMerger(void)  {}
  ~BWTFragMerger()  {}
  
  //void AlignToFragType(BWT &bwt, std::vector<AlignType> &all_pos, std::vector<FragType> &all_frag, const int num_hits);

  /* obsolete  
  void MergeFragments(
    std::vector<FragType> &all_frag, std::vector<FragType> &merged_frag,
    int mismatch = 2, int gap = 2
  );
  */

  void MergeFragments(
      BWT &bwt, const int query_len,
      std::vector<MatchType> &frags, std::vector<MatchType> &merged_frags,
      const int gap = 5
  );

  // frag_cost is the additional cost for adding a new fragment to the matching
  int FindBestMatching(
      const std::string &query_seq, BWT &bwt, BWT &bwt_aux, const int num_fragments,    
      std::vector<MatchType> &all_frag,
      std::vector<FragType> &selected,
      const bool ensure_partition, const bool use_insert_penalty,
      const int frag_cost = 10, const int max_overlap = 5, const float min_map_portion = 0.6
  );
  
  int BoundaryCost(const int begin, const int end, const int length, const int len_forgiven = 2)    {
    // computing the boundary cost, the boundary is consdered as the shorter boundary between the two directions
    // we assume that the other part could be missing 
    int n1 = begin - len_forgiven; int n2 = length - end - len_forgiven;
    n1 = n1 >= 0 ? n1 : 0; n2 = n2 >= 0 ? n2 : 0;
    int n = n1 < n2 ? n1 :  n2;
    return -(n * n);
  }
};

#endif
