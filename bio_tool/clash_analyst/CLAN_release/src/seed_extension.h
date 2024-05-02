#ifndef _SEED_EXTENSION_H_
#define _SEED_EXTENSION_H_

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
#include <sstream>

#include "score.h"

class SeedExtension {
 public:

  template <typename S1IndexType, typename S2IndexType>
  int UngappedExt(
    const char *s1, S1IndexType &l1, S1IndexType &r1, S1IndexType len1,  
    const char *s2, S2IndexType &l2, S2IndexType &r2, S2IndexType len2,
    int &num_mismatch,                  // number of mismatches recorded, used for computing sequence identity
    const bool check_special_char = false,  // if meet special character, terminate extension
    const char special_char = (char) 1,     // the special character, ignored when check_special_char set to false
    const int lower_bound = -5              // the lower bound for ungapped extension
  ) {
    S2IndexType tl2 = l2, tr2 = r2;
    //std::cerr << "extending: " << l1 << " " << r1 << " " << l2 << "  " << r2 << std::endl;
    //std::cerr << s1 + l1 << std::endl;
    //std::cerr << s2 + l2 << std::endl;
    assert(lower_bound < 0);
    // extend to left
    int curr_score = 0, opt_score = 0, curr_mismatch = 0, opt_mismatch = 0;
    S1IndexType opt_l1 = l1; S2IndexType opt_l2 = l2;
    while(l1 > 0 && l2 > 0 && curr_score > lower_bound) {    
      -- l1; -- l2;
      // check if special character presents
      if(check_special_char && (s1[l1] == special_char || s2[l2] == special_char)) {
        //std::cerr << "special char found:  " << s1[l1] << " " << s2[l2] << "  " << special_char << std::endl;
        break;
      }
      curr_score += s1[l1] == s2[l2] ? NT_MATCH : NT_MISMATCH;
      curr_mismatch += s1[l1] == s2[l2] ? 0 : 1;
      if(curr_score >= opt_score)  {  // if extension is able to be made, record it even with the same score
        opt_score = curr_score; opt_mismatch = curr_mismatch, opt_l1 = l1; opt_l2 = l2;
      }
      //if(tl2 == 23 && tr2 == 54)
      //  std::cerr << l1 << " " << l2 << "  " << s1[l1] << "  " << s2[l2] << "  " << curr_score << " " << opt_score << std::endl;
    }   
    int score_gain = opt_score;
    num_mismatch += opt_mismatch;
    
    // extend to right
    curr_score = 0, opt_score = 0, curr_mismatch = 0, opt_mismatch = 0;
    S1IndexType opt_r1 = r1; S2IndexType opt_r2 = r2;
    //std::cerr << r1 << "  " << len1 << "  " << r2 << "  " << len2  << std::endl;
    while(r1 < len1 - 1 && r2 < len2 - 1 && curr_score > lower_bound) {
      ++ r1; ++ r2;
      // check if special character presents
      if(check_special_char && (s1[r1] == special_char || s2[r2] == special_char)) {
        //std::cerr << "special char found:  " << s1[r1] << " " << s2[r2] << "  " << special_char << std::endl;
        break;
      }
      curr_score += s1[r1] == s2[r2] ? NT_MATCH : NT_MISMATCH;
      curr_mismatch += s1[r1] == s2[r2] ? 0 : 1;
      if(curr_score >= opt_score)  {  // if extension is able to be made, record it even with the same score
        opt_score = curr_score; opt_mismatch = curr_mismatch; opt_r1 = r1; opt_r2 = r2;
      }
      //if(tl2 == 0 && tr2 == 9)
      //  std::cerr << r1 << " " << r2 << "  " << s1[r1] << "  " << s2[r2] << "  " << curr_score << " " << opt_score << std::endl;
    }  
    // update interval
    l1 = opt_l1; r1 = opt_r1; l2 = opt_l2; r2 = opt_r2; 
    num_mismatch += opt_mismatch;
    return score_gain + opt_score;
  }
  
  template <typename S1IndexType, typename S2IndexType>
  int GappedExt(
    const char *s1, S1IndexType &l1, S1IndexType &r1, S1IndexType len1,  
    const char *s2, S2IndexType &l2, S2IndexType &r2, S2IndexType len2,
    char *cigar,                 // the cigar string that is going to be updated with gaps
    int &num_mismatch,           // number of mismatches recorded, used for computing sequence identity
    const bool check_special_char = false,  // if meet special character, terminate extension
    const char special_char = (char) 1,     // the special character, ignored when check_special_char set to false
    const int max_gap_size = 2,         // the maximum size of the gap
    const int lower_bound = -5          // the lower bound for ungapped extension
  ) {
    if(max_gap_size <= 0) return 0;
    assert(lower_bound < 0);
    // record the original terminals
    S1IndexType ol1 = l1, or1 = r1;
    S2IndexType ol2 = l2, or2 = r2;
    // determin the gap penalty
    int gap_penalty;
    // extend to left
    int curr_score = 0, opt_score = 0, curr_mismatch = 0, opt_mismatch = 0;
    S1IndexType opt_l1 = l1; S2IndexType opt_l2 = l2;
    int lsize = 0; int lext = 0; char loperation = 'M';
    for(int i = 1; i <= max_gap_size; ++ i) {
      gap_penalty = NT_GAP_OPEN + i * NT_GAP_EXT;
      // open gap in S1
      // reset the boundaries and score
      l1 = ol1 - i; l2 = ol2; curr_score = 0; curr_mismatch = 0;
      while(l1 > 0 && l2 > 0 && curr_score > lower_bound) {    
        -- l1; -- l2;
        if(check_special_char && (s1[l1] == special_char || s2[l2] == special_char)) {
          //std::cerr << "special char found:  " << s1[l1] << " " << s2[l2] << "  " << special_char << std::endl;
          break;
        }
        curr_score += s1[l1] == s2[l2] ? NT_MATCH : NT_MISMATCH;
        curr_mismatch += s1[l1] == s2[l2] ? 0 : 1;
        if(curr_score + gap_penalty >= opt_score)  {  // if extension is able to be made, record it even with the same score
          opt_score = curr_score + gap_penalty; opt_mismatch = curr_mismatch; 
          opt_l1 = l1; opt_l2 = l2;
          lsize = i; loperation = 'D'; lext = ol2 - l2;
        }
        //if(tl2 == 23 && tr2 == 54)
        //  std::cerr << l1 << " " << l2 << "  " << s1[l1] << "  " << s2[l2] << "  " << curr_score << " " << opt_score << std::endl;
      }   
      // open gap in S2
      l1 = ol1; l2 = ol2 - i; curr_score = 0; curr_mismatch = 0;
      while(l1 > 0 && l2 > 0 && curr_score > lower_bound) {    
        -- l1; -- l2;
        if(check_special_char && (s1[l1] == special_char || s2[l2] == special_char)) {
          //std::cerr << "special char found:  " << s1[l1] << " " << s2[l2] << "  " << special_char << std::endl;
          break;
        }
        curr_score += s1[l1] == s2[l2] ? NT_MATCH : NT_MISMATCH;
        curr_mismatch += s1[l1] == s2[l2] ? 0 : 1;
        if(curr_score + gap_penalty >= opt_score)  {  // if extension is able to be made, record it even with the same score
          opt_score = curr_score + gap_penalty; opt_mismatch = curr_mismatch;
          opt_l1 = l1; opt_l2 = l2;
          lsize = i; loperation = 'I'; lext = ol1 - l1;
        }
      }
    }
    int score_gain = opt_score;
    num_mismatch += opt_mismatch;
    //std::cerr << "score_gain:  " << score_gain << " location: " << ol1 << " " << ol2 << " " << opt_l1 << "  " << opt_l2 << std::endl;
    // extend to right
    curr_score = 0; opt_score = 0; curr_mismatch = 0; opt_mismatch = 0;
    S1IndexType opt_r1 = r1; S2IndexType opt_r2 = r2;
    int rsize = 0; int rext = 0; char roperation = 'M';
    for(int i = 1; i <= max_gap_size; ++ i) {
      gap_penalty = NT_GAP_OPEN + i * NT_GAP_EXT;
      // open gap in S1
      // reset the boundaries and score
      r1 = or1 + i; r2 = or2; curr_score = 0; curr_mismatch = 0;
      while(r1 < len1 - 1 && r2 < len2 - 1 && curr_score > lower_bound) {
        ++ r1; ++ r2;
        if(check_special_char && (s1[r1] == special_char || s2[r2] == special_char)) {
          //std::cerr << "special char found:  " << s1[r1] << " " << s2[r2] << "  " << special_char << std::endl;
          break;
        }
        curr_score += s1[r1] == s2[r2] ? NT_MATCH : NT_MISMATCH;
        curr_mismatch += s1[r1] == s2[r2] ? 0 : 1;
        if(curr_score + gap_penalty >= opt_score)  {  // if extension is able to be made, record it even with the same score
          opt_score = curr_score + gap_penalty; opt_mismatch = curr_mismatch;
          opt_r1 = r1; opt_r2 = r2;
          rsize = i; roperation = 'D'; rext = r2 - or2;
        }
        //std::cerr << r1 << " " << r2 << "  " << s1[r1] << "  " << s2[r2] << "  " << curr_score << " " << opt_score << std::endl;
      }    
      // open gap in S2
      r1 = or1; r2 = or2 + i; curr_score = 0, curr_mismatch = 0;
      while(r1 < len1 - 1 && r2 < len2 - 1 && curr_score > lower_bound) {
        ++ r1; ++ r2;
        if(check_special_char && (s1[r1] == special_char || s2[r2] == special_char)) {
          //std::cerr << "special char found:  " << s1[r1] << " " << s2[r2] << "  " << special_char << std::endl;
          break;
        }
        curr_score += s1[r1] == s2[r2] ? NT_MATCH : NT_MISMATCH;
        curr_mismatch += s1[r1] == s2[r2] ? 0 : 1;
        if(curr_score + gap_penalty >= opt_score)  {  // if extension is able to be made, record it even with the same score
          opt_score = curr_score + gap_penalty; opt_mismatch = curr_mismatch; 
          opt_r1 = r1; opt_r2 = r2;
          rsize = i; roperation = 'I'; rext = r1 - or1;
        }
        //std::cerr << r1 << " " << r2 << "  " << s1[r1] << "  " << s2[r2] << "  " << curr_score << " " << opt_score << std::endl;
      } 
    }
    score_gain += opt_score;
    num_mismatch += opt_mismatch;
    l1 = opt_l1; r1 = opt_r1; l2 = opt_l2; r2 = opt_r2;
    // update the CIGAR string
    //std::stringstream sl, sr;
    // DEBUG
    if(lsize > 0)  {
      //std::cerr << "L before CIGAR: " << cigar << std::endl;
      //std::cerr << "L components: " << lext << "  " << lsize << " " << loperation << std:: endl;
      char *ncigar = new char [100];
      sprintf(ncigar, "%d%c%d%c%s", lext, 'M', lsize, loperation, cigar);  
      strcpy(cigar, ncigar);
      //std::cerr << "L after CIGAR: " << cigar << std::endl;
      delete [] ncigar;
    }
    if(rsize > 0) {
      //std::cerr << "R before CIGAR: " << cigar << std::endl;
      //std::cerr << "R components: " << rext << "  " << rsize << " " << roperation << std:: endl;
      char *ncigar = new char [100];
      sprintf(ncigar, "%s%d%c%d%c", cigar, rsize, roperation, rext, 'M');
      strcpy(cigar, ncigar);
      //std::cerr << "R after CIGAR: " << cigar << std::endl;
      delete [] ncigar;
    }
    return score_gain;
  }
  
};


#endif
