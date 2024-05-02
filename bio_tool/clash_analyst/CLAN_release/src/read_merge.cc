#include "read_merge.h"

using namespace std;

void ReadMerge::MergeReads(
  const int num_ref, char const* const* ref_header, char const* const* ref_seq,
  const int ave_read_len,
  std::vector< std::vector<FragType> > &mapping,
  std::vector<GPosType> &islands,
  std::vector< std::unordered_map<int, float> > &duplex,
  const CLANAnnotateParam &param_set
) {
  Init(num_ref, ref_seq);
  //cerr << "Init array done" << endl;
  DistributeReads(mapping, param_set.strand_specific);
  //cerr << "Distribute reads done" << endl; 
  DetectIslands(
    islands, ave_read_len, param_set.min_island_len, param_set.max_island_len, 
    param_set.max_island_gap, param_set.strand_specific
  );
  //cerr << "Detect islands done" << endl;
  FindDuplex(mapping, islands, ave_read_len, duplex, param_set.strand_specific, param_set.min_coverage);
  //cerr << "Find duplex done" << endl;
  //ReleaseCoverage();
  //cerr << "Release memory done" << endl;
  return;
}

void ReadMerge::DistributeReads(std::vector< std::vector<FragType> > &mapping, const bool strand_specific) {

  for(int i = 0; i < mapping.size(); ++ i) {
    if(mapping[i].size() <= 0) continue;
    float cov_sol = INTPREC * 2 / mapping[i].size();    // multiplied by 2 because the solutions come in pair as both RNA arms
    for(int j = 0; j < mapping[i].size(); ++ j) {
      if(mapping[i][j].targets.size() <= 0) continue;
      float cov_mapping = cov_sol / mapping[i][j].targets.size();
      for(int k = 0; k < mapping[i][j].targets.size(); ++ k) {
        //cerr << "locations: " << mapping[i][j].targets[k].begin << "  " << mapping[i][j].targets[k].end << endl;
        for(int l = mapping[i][j].targets[k].begin; l <= mapping[i][j].targets[k].end; ++ l) {
          //cerr << "cov: " << coverage_pos_[mapping[i][j].targets[k].sid][l] << "  " << coverage_rev_[mapping[i][j].targets[k].sid][l] << endl;
          //cerr << l << "  " << ref_len_[mapping[i][j].targets[k].sid] << endl; 
          //if(mapping[i][j].targets[k].sid == 775) cerr << "Reads distributed: " << (int) cov_mapping << endl;
          
          if(!strand_specific || mapping[i][j].targets[k].strand)  {
            coverage_pos_[mapping[i][j].targets[k].sid][l] += (int) cov_mapping;
            //cerr << "coverage:  " << coverage_pos_[mapping[i][j].targets[k].sid][l] << endl;
          } else  {
            coverage_rev_[mapping[i][j].targets[k].sid][l] += (int) cov_mapping;
          }
        }        
      }
    }
  }
  return;
}

BWTIDX ReadMerge::CheckMaxIslandLen(
  const std::vector<int> coverage, const int min_coverage, const BWTIDX begin, const BWTIDX end
) {
  BWTIDX curr_len = 0, max_len = 0;
  bool is_open = false;
  BWTIDX i;
  for(i = begin; i <= end; ++ i) {
    if(is_open && coverage[i] >= min_coverage)  {
      ++ curr_len;
    } else if(!is_open && coverage[i] >= min_coverage) {
      curr_len = 1; is_open = true;
    } else if(is_open && coverage[i] < min_coverage) {
      is_open = false; 
      if(curr_len > max_len)  max_len = curr_len;
    } // if island is not open and the coverage is less than threashold, do nothing
  }
  if(is_open && curr_len > max_len) max_len = curr_len; // record the last one
  return max_len;
}

void ReadMerge::PartitionIslands(
  const GPosType &candidate, const int ave_read_len, 
  const int min_island_len, const int max_island_len, std::vector<GPosType> &islands
) {
  // sort the coverages
  BWTIDX len = candidate.end - candidate.begin + 1;
  vector<int> cov_holder(len, 0);
  for(BWTIDX i = 0; i < len; ++ i) {
    cov_holder[i] = candidate.strand ? 
      coverage_pos_[candidate.sid][candidate.begin + i] : coverage_rev_[candidate.sid][candidate.begin + i];
  }
  sort(cov_holder.begin(), cov_holder.end());
  // using binary search to find the minimum coverage that partition the entire island into smaller ones;
  // such that no smaller island is longer than max size allowed
  BWTIDX left = 0, right = len - 1, middle;
  BWTIDX curr_max_len;
  while(left < right) {
    middle = floor((left + right) / 2);
    if(candidate.strand)  {
      curr_max_len = CheckMaxIslandLen(coverage_pos_[candidate.sid], cov_holder[middle], candidate.begin, candidate.end);
    } else  {
      curr_max_len = CheckMaxIslandLen(coverage_rev_[candidate.sid], cov_holder[middle], candidate.begin, candidate.end);
    }
    if(curr_max_len > max_island_len) {
      // need to use even higher coverage threshold
      left = middle + 1;
    } else if(curr_max_len < max_island_len) {
      // check if we can use lower coverage threshold to preserve more information
      right = middle - 1;
    } else  {
      break;
    }
  }
  // now cut the island and record the partitioned islands into the output
  int cut_threshold = curr_max_len <= max_island_len ? cov_holder[middle] : cov_holder[middle + 1];
  GPosType curr_island = candidate;
  BWTIDX curr_len = 0; BWTIDX total_cov = 0;
  bool is_open = false;
  BWTIDX i;
  for(i = candidate.begin; i <= candidate.end; ++ i) {
    if(is_open && 
        (( candidate.strand && coverage_pos_[candidate.sid][i] >= cut_threshold) || 
         (!candidate.strand && coverage_rev_[candidate.sid][i] >= cut_threshold))
    )  {
      ++ curr_len;
      total_cov += candidate.strand ? coverage_pos_[candidate.sid][i] : coverage_rev_[candidate.sid][i];
    } else if(!is_open && 
        (( candidate.strand && coverage_pos_[candidate.sid][i] >= cut_threshold) || 
         (!candidate.strand && coverage_rev_[candidate.sid][i] >= cut_threshold))
    ) {
      curr_len = 1; is_open = true;
      total_cov = candidate.strand ? coverage_pos_[candidate.sid][i] : coverage_rev_[candidate.sid][i];
    } else if(is_open && 
        (( candidate.strand && coverage_pos_[candidate.sid][i] < cut_threshold) || 
         (!candidate.strand && coverage_rev_[candidate.sid][i] < cut_threshold))
    ) {
      is_open = false; 
      if(curr_len > min_island_len)  {
        curr_island.end = i - 1;
        curr_island.begin = curr_island.end - curr_len + 1;
        curr_island.coverage = total_cov / (INTPREC * ave_read_len);
        islands.push_back(curr_island);
      }
    } // if island is not open and the coverage is less than threashold, do nothing
  }
  if(curr_len < min_island_len) {
    curr_island.end = i - 1;
    curr_island.begin = curr_island.end - curr_len + 1;
    curr_island.coverage = total_cov / (INTPREC * ave_read_len);
    islands.push_back(curr_island);
  }
  return;
}


void ReadMerge::IncorporateIslands(
  const GPosType &candidate,
  std::vector<GPosType> &islands, const int ave_read_len,
  const int min_island_len, const int max_island_len,
  const int window_size
) {
  //cerr << "IncorporateIsland called" << endl;
  BWTIDX len = candidate.end - candidate.begin + 1;
  //cerr << "check length:  " << len << " " << min_island_len << endl;
  if(len < min_island_len)  return; // case 1: the island dos not satisify minimum length requirement
  else if(len >= min_island_len & len <= max_island_len)  { // case 2: the island has the right size, record it
    GPosType is = candidate;
    is.coverage /= (INTPREC * ave_read_len);
    //cerr << "island coverage: " << is.coverage << " " << ave_read_len << endl;
    islands.push_back(is);
    for(BWTIDX k = is.begin; k <= is.end; ++ k) {
      island_hash_pos_[is.sid][k] = islands.size() - 1;
    }
  } else if(len > max_island_len)  {  // case 3: the island is too long, needs partitioning
    PartitionIslands(candidate, ave_read_len, min_island_len, max_island_len, islands);  
  }
  return;
}

void ReadMerge::DetectIslands(
    std::vector<GPosType> &islands, const int ave_read_len,
    const int min_island_len, const int max_island_len, 
    const int max_gap_size, const bool strand_specific
) {
  for(int i = 0; i < num_ref_; ++ i) {
    // take care of positive strand or non strand-specific cases
    GPosType current_pos; current_pos.sid = i; current_pos.begin = -1; current_pos.strand = true; current_pos.coverage = 0;
    int uncovered_len = 0;
    bool island_open = false;
    for(BWTIDX j = 0; j < ref_len_[i]; ++ j) {
      //cerr << "pos: " << j << "  " << coverage_rev_[i][j] << endl;
      if(coverage_pos_[i][j] > 0) {
        if(island_open) { current_pos.end = j; current_pos.coverage += coverage_pos_[i][j]; }
        else  { island_open = true; current_pos.begin = current_pos.end = j; current_pos.coverage = coverage_pos_[i][j]; uncovered_len = 0; }
      } else  {
        // only need to consider if island is open
        if(island_open)  {
          ++ uncovered_len;
          if(uncovered_len > max_gap_size)  {
            //cerr << "Called here" << endl;
            IncorporateIslands(current_pos, islands, ave_read_len, min_island_len, max_island_len);
            // terminate the current island
            island_open = false;
          }
        }
      }
    }
    // if the current island remains open
    if(island_open && current_pos.end - current_pos.begin > min_island_len)  {
      IncorporateIslands(current_pos, islands, ave_read_len, min_island_len, max_island_len);
    }
  }
  // quit if we do not need to consider strand specific case
  if(!strand_specific)  return;
  
  // otherwise redo the same on negative strand
  for(int i = 0; i < num_ref_; ++ i) {
    // take care of positive strand or non strand-specific cases
    GPosType current_pos; current_pos.sid = i; current_pos.begin = -1; current_pos.strand = false;
    int uncovered_len = 0;
    bool island_open = false;
    for(BWTIDX j = 0; j < ref_len_[i]; ++ j) {
      //cerr << "rev: " << j << "  " << coverage_rev_[i][j] << endl;
      if(coverage_rev_[i][j] > 0) {
        if(island_open) { current_pos.end = j; current_pos.coverage += coverage_rev_[i][j]; }
        else  { island_open = true; current_pos.begin = current_pos.end = j; current_pos.coverage = coverage_rev_[i][j]; uncovered_len = 0; }
      } else  {
        // only need to consider if island is open
        if(island_open)  {
          ++ uncovered_len;
          if(uncovered_len > max_gap_size)  {
            IncorporateIslands(current_pos, islands, ave_read_len, min_island_len, max_island_len);
            island_open = false;
          }
        }
      }
    }
    // if the current island remains open
    if(island_open && current_pos.end - current_pos.begin > min_island_len)  {
      IncorporateIslands(current_pos, islands, ave_read_len, min_island_len, max_island_len);
    }
  }
  return;
}

void ReadMerge::FindDuplex(
  std::vector< std::vector<FragType> > &mapping, 
  std::vector<GPosType> &islands, const int ave_read_len,
  std::vector< std::unordered_map<int, float> > &duplex,
  const bool strand_specific, const int min_coverage
) {
  duplex.resize(islands.size());
  for(int i = 0; i < mapping.size(); ++ i) {
    if(mapping[i].size() <= 0) continue;
    float cov_sol = INTPREC * 2 / mapping[i].size();
    //cerr << "cov_sol:  " << cov_sol << "  " << mapping[i].size() << endl;
    // the solutions come in pair
    for(int j = 0; j < mapping[i].size(); ++ j) {
      // check if the mapping involves two strands; if not, it is not a duplex read and we should disgard it
      // TODO: deal with mapping to a single region
      if(j >= mapping[i].size() - 1 || mapping[i][j].sol_id != mapping[i][j + 1].sol_id)  {
        continue;
      }
      // otherwise the mapping forms a pair
      int n = mapping[i][j].targets.size() < mapping[i][j + 1].targets.size() ? mapping[i][j].targets.size() : mapping[i][j + 1].targets.size();
      //cerr << "n: " << n << " " << mapping[i][j].targets.size() << "  " << mapping[i][j + 1].targets.size() << endl;
      if(n <= 0) continue;
      float cov_duplex = cov_sol / n;
      //cerr << "cov_duplex:  " << cov_duplex << "  " << n << endl;
      //for(auto it = island_hash_pos_.begin(); it != island_hash_pos_.end(); ++ it) {
      //  cerr << "In hash: " << it->first << " " << it->second << endl;
      //}
      
      int island_x = -1, island_y = -1;
      for(int k = 0; k < mapping[i][j].targets.size(); ++ k) {
        if(!strand_specific || mapping[i][j].targets[k].strand)  {
          if(island_hash_pos_[mapping[i][j].targets[k].sid].find(mapping[i][j].targets[k].begin) 
              == island_hash_pos_[mapping[i][j].targets[k].sid].end()) continue;
          island_x = island_hash_pos_[mapping[i][j].targets[k].sid][mapping[i][j].targets[k].begin];
        } else  {
          if(island_hash_rev_[mapping[i][j].targets[k].sid].find(mapping[i][j].targets[k].begin) 
              == island_hash_rev_[mapping[i][j].targets[k].sid].end()) continue;
          island_x = island_hash_rev_[mapping[i][j].targets[k].sid][mapping[i][j].targets[k].begin];
        }
        for(int l = 0; l < mapping[i][j + 1].targets.size(); ++ l) {
          if(!strand_specific || mapping[i][j + 1].targets[l].strand)  {
            if(island_hash_pos_[mapping[i][j + 1].targets[l].sid].find(mapping[i][j + 1].targets[l].begin) 
              == island_hash_pos_[mapping[i][j + 1].targets[l].sid].end()) continue;
            island_y = island_hash_pos_[mapping[i][j + 1].targets[l].sid][mapping[i][j + 1].targets[l].begin];
          } else  {
            if(island_hash_rev_[mapping[i][j + 1].targets[l].sid].find(mapping[i][j + 1].targets[l].begin) 
              == island_hash_rev_[mapping[i][j + 1].targets[l].sid].end()) continue;
            island_y = island_hash_rev_[mapping[i][j + 1].targets[l].sid][mapping[i][j + 1].targets[l].begin];
          }
          
          //cerr << "islands located" << endl;
          
          if(island_x < 0 || island_y < 0) continue;
          
          int i1 = island_x < island_y ? island_x : island_y;
          int i2 = island_x < island_y ? island_y : island_x;
          
          
          if(duplex[i1].find(i2) == duplex[i1].end())  {
            duplex[i1][i2] = 0;
          }
          duplex[i1][i2] += (int) cov_duplex;
                 
        }  
      }     
      // remember to further increase index j
      ++ j;
    }
  }
  for(int i = 0; i < duplex.size(); ++ i) {
    unordered_map<int, int> to_delete;
    for(auto it = duplex[i].begin(); it != duplex[i].end(); ++ it) {
      //cerr << "duplex coverage: " << it->second << "  " << min_coverage << endl;  
      it->second /= INTPREC;    
      //cerr << "duplex coverage: " << it->second << "  " << min_coverage << endl;  
      if(it->second < min_coverage) to_delete[it->first] = 1;
    }
    for(auto it = to_delete.begin(); it != to_delete.end(); ++ it) {
      duplex[i].erase(it->first);
    }
  }
  return;
}


