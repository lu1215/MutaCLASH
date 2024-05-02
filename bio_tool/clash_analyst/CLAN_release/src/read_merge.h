#ifndef _READ_MERGE_H_
#define _READ_MERGE_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <list>
#include <queue>
#include <stack>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <tuple>
#include <unordered_map>

#include "bwt_frag_merger.h"
#include "parameters.h"

//TODO: consider strand information

#ifndef INTPREC
#define INTPREC 1000
#endif

struct GPosType {
  int sid;
  BWTIDX begin, end;
  bool strand;
  float coverage;
  
  void operator=(const GPosType &gp)  {
    sid = gp.sid;
    begin = gp.begin;
    end = gp.end;
    strand = gp.strand;
    coverage = gp.coverage;
  }
};

// TODO: change interface to coverage rather than DuplexInfo
// Move DuplexInfo definition to dimer_folding.h

class ReadMerge {
 public:
 
  ReadMerge() {}
  ~ReadMerge() {}
  
  void MergeReads(
    const int num_ref, char const* const* ref_header, char const* const* ref_seq,
    const int ave_read_len,
    std::vector< std::vector<FragType> > &mapping,
    std::vector<GPosType> &islands,
    std::vector< std::unordered_map<int, float> > &duplex,
    const CLANAnnotateParam &param_set
  );
  
  void PrintIslands(char const* const* ref_header, std::vector<GPosType> &islands) {
    for(int i = 0; i < islands.size(); ++ i) {
      std::cout << ref_header[islands[i].sid] << "  " <<
       (long long int) islands[i].begin << " " << (long long int) islands[i].end << 
       std::endl;
    }
    return;
  }
  
  void PrintDuplex(char const* const* ref_header, std::vector<GPosType> &islands, std::vector< std::unordered_map<int, int> > &duplex) {
    for(int i = 0; i < duplex.size(); ++ i) {
      for(auto it = duplex[i].begin(); it != duplex[i].end(); ++ it) {
        std::cout << ref_header[islands[i].sid] << " " << ref_header[islands[it->first].sid] << " " << it->second << std::endl;
      }
    }
    return;
  }
  
  
 private:
  void Init(const int num_ref, char const* const* ref_seq)  {
    num_ref_ = num_ref;
    coverage_pos_.resize(num_ref_);
    coverage_rev_.resize(num_ref_);
    island_hash_pos_.resize(num_ref_);
    island_hash_pos_.resize(num_ref_);
    ref_len_.resize(num_ref_);
    for(int i = 0; i < num_ref_; ++ i) {
      ref_len_[i] = std::string(ref_seq[i]).length();
      if(ref_len_[i] >  0)  {
        coverage_pos_[i].resize(ref_len_[i], 0);
        coverage_rev_[i].resize(ref_len_[i], 0);
        //memset(coverage_pos_[i], 0, sizeof(int) * ref_len_[i]); 
        //memset(coverage_rev_[i], 0, sizeof(int) * ref_len_[i]);
      }
    }
    //coverage_allocated_ = true;
    //std::cerr << num_ref_ << std::endl;
    return;
  }
  
  /*
  void ReleaseCoverage()  {
    std::cerr << num_ref_ << std::endl;
    if(!coverage_allocated_)  return;
    for(int i = 0; i < num_ref_; ++ i) {
      //std::cerr << "length:  " <<  i << "  " << ref_len_[i] << std::endl;
      if(ref_len_[i] > 0) {  
        delete [] coverage_pos_[i];
        delete [] coverage_rev_[i];
      }
      //std::cerr << "done" << std::endl;
    }
    //std::cerr << "finish second dimension" << std::endl;
    delete [] coverage_pos_;
    delete [] coverage_rev_;
    //std::cerr << "finish first dimension" << std::endl;
    coverage_allocated_ = false;
    return;
  }
  */
  
  // distribute the reads against the reference genomes
  void DistributeReads(std::vector< std::vector<FragType> > &mapping, const bool strand_specific = false);
  
  // define the islands from the coverage information
  void DetectIslands(
      std::vector<GPosType> &islands, const int ave_read_len,
      const int min_island_len = 20, const int max_island_len = 200, 
      const int max_gap_size = 10, const bool strand_specific = false
  );

  void FindDuplex(
    std::vector< std::vector<FragType> > &mapping, 
    std::vector<GPosType> &islands, const int ave_read_len,
    std::vector< std::unordered_map<int, float> > &duplex,
    const bool strand_specific = false, const int min_coverage = 3
  );
  
  void PartitionIslands(
      const GPosType &candidate, const int ave_read_len, 
      const int min_island_len, const int max_island_len, std::vector<GPosType> &islands
  );
  // if we ingore sites with coverage less than min_coverage, what is the longest island we get
  BWTIDX CheckMaxIslandLen(const std::vector<int> coverage, const int min_coverage, const BWTIDX begin, const BWTIDX end);
  
  void IncorporateIslands(
    const GPosType &candidate,
    std::vector<GPosType> &islands, const int ave_read_len,
    const int min_island_len = 20, const int max_island_len = 200,
    const int window_size = 10
  );

  int num_ref_;
  std::vector<BWTIDX> ref_len_;
  bool coverage_allocated_;
  std::vector< std::vector<int> > coverage_pos_;
  std::vector< std::vector<int> > coverage_rev_;
  std::vector< std::unordered_map<BWTIDX, int> > island_hash_pos_;
  std::vector< std::unordered_map<BWTIDX, int> > island_hash_rev_;
  
};



#endif
