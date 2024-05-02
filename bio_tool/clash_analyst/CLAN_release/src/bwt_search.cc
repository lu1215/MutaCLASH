#include "bwt_search.h"

using namespace std;

// sorting the AlignType
bool SortAlignTypeBySize(const AlignType &a, const AlignType &b)  {
  int a_dist = a.q_end - a.q_begin;
  int b_dist = b.q_end - b.q_begin;
  if(a_dist > b_dist || (a_dist == b_dist && a.q_begin < b.q_begin))  {
    return true;
  }  
  return false;
} 

void BWTSearch::SearchExact(BWT &bwt, const char *str, AlignType &pos) {
  
  int len = strlen(str);
  if(len <= 0)  {
    cerr << "Searching empty sequence. An empty list will be returned." << endl;
    return;
  }
  pair<BWTIDX, BWTIDX> range;
  range.first = 0, range.second = bwt.GetSize();
  int i;
  for(i = strlen(str) - 1; i >= 0; -- i) {
    // indicating that there is no solution
    if(range.first >= range.second)  break;
    range = bwt.UpdateRange(str[i], range);   
    //cerr << "BWT Update:  " << i << "  " << range.first << " " << range.second << endl; 
  }
  // record the positions
  pos.q_begin = i + 1; pos.q_end = strlen(str) - 1; pos.cost = 0;
  pos.bwt_begin = range.first; pos.bwt_end = range.second;
  return;
}

void BWTSearch::SearchAllSubRegions(
    BWT &bwt, const int min_len, const int max_hits, const int flank_len,
    const char *str, std::vector<AlignType> &all_pos
) {
  //cout << min_len << "  " << max_hits << "  " << max_sub_hits << "  " << flank_len << endl;
  //cerr << "sequence:    " << str << endl;
  int len = strlen(str);
  if(len <= 0)  {
    cerr << "BWTSearch::Search: Searching empty sequence. An empty list will be returned." << endl; return;
  }
  BWTIDX bwt_size = bwt.GetSize();
  // reverse the string
  string s = str; s = string(s.rbegin(), s.rend());
  //cout << "Query sequence:  " << s << endl;
  // construct a vector of ranges
  vector<AlignType> cand_pos; cand_pos.resize(len);
  vector<queue<AlignType> > record_range; record_range.resize(len);
  for(int i = 0; i < cand_pos.size(); ++ i) {
    //cout << "=======================" << endl;
    cand_pos[i].bwt_begin = 0; cand_pos[i].bwt_end = bwt_size;
    cand_pos[i].q_begin = cand_pos[i].q_end = len - i - 1; 
    ++ cand_pos[i].q_begin; // this is because we need to subtract 1 for search of the first character
    //cerr << "i: " << i << " q_end: " << cand_pos[i].q_end << endl;

    //cout << i << "  " << s[i] << endl;
    for(int j = 0; j <= i; ++ j) {
      if(cand_pos[j].bwt_begin >= cand_pos[j].bwt_end) continue;
      BWTIDX x, y; 
      bwt.UpdateRange(s[i], cand_pos[j].bwt_begin, cand_pos[j].bwt_end, x, y);
      //if(cand_pos[j].q_end == 9)  {
        //cerr << "pushed: " << i << "  " << j  << "  " << s[i] << " query location:  " << cand_pos[j].q_begin << "  " << cand_pos[j].q_end << " " << x << " " << y << endl;
      //}
      if(x < y) cand_pos[j].q_begin --; // if the range is valid, extend the matched range
      cand_pos[j].bwt_begin = x; cand_pos[j].bwt_end = y;
      if(cand_pos[j].bwt_begin < cand_pos[j].bwt_end && cand_pos[j].q_end - cand_pos[j].q_begin + 1 >= min_len) {
        //cerr << len - i << "  " << len - j << " pushed... " << (long long int) cand_pos[j].bwt_begin << "  " << (long long int) cand_pos[j].bwt_end << endl;
        // check if the queue is full        
        if(record_range[j].size() >= flank_len)  {
            //cerr << "popped: " << record_range[j].front().q_begin << " " << record_range[j].front().q_end << endl;
            record_range[j].pop();  
        }
        record_range[j].push(cand_pos[j]);   
      }
    }
  }
  // record the final identified ranges
  for(int i = 0; i < cand_pos.size(); ++ i) {
    if(cand_pos[i].bwt_begin < cand_pos[i].bwt_end && cand_pos[i].q_end - cand_pos[i].q_begin + 1 >= min_len) {
      if(record_range[i].size() >= flank_len)  record_range[i].pop();  
      record_range[i].push(cand_pos[i]); 
      //cerr << "pushed: " << cand_pos[i].q_begin << "  " << cand_pos[i].q_end << " " << (long long int) cand_pos[i].bwt_begin << " " << (long long int) cand_pos[i].bwt_end << endl;
    }
  }
  
  // finally record all ranges
  for(auto it = record_range.begin(); it != record_range.end(); ++ it) {
    //cerr << "identified seeds:  " << it->front().q_begin << "  " << it->front().q_end << endl;  
    if(it->empty() || (long long int) it->back().bwt_end - (long long int) it->back().bwt_begin + 1 > max_hits) {  continue; }
    while(!it->empty()) {
      //cerr << "identified seeds:  " << it->front().q_begin << "  " << it->front().q_end << endl;  
      if((long long int) it->front().bwt_end - (long long int) it->front().bwt_begin + 1 <= max_hits) {
        //cerr << "identified seeds:  " << it->front().q_begin << "  " << it->front().q_end << endl; 
        all_pos.push_back(it->front()); 
      }
      it->pop();
    }
  }
  // DEBUG check all searched sequences
  /*
  cout << "size of fragments: " << all_pos.size() << endl;
  string sp = str;
  for(int i = 0; i < all_pos.size(); ++ i) {
    cerr << "HIT ENTRY: " << all_pos[i].q_begin << "    " << all_pos[i].q_end << "  " << all_pos[i].bwt_begin << "    " << all_pos[i].bwt_end << endl;
    for(int j = all_pos[i].q_begin; j <= all_pos[i].q_end; ++ j)   {
        cerr << sp[j];
    }
    cerr << endl;
    bwt.PrintSubstr(bwt.BWTToStrPos(all_pos[i].bwt_begin), all_pos[i].q_end - all_pos[i].q_begin + 1);
  }
  */
  return;
}

void BWTSearch::ExtendAllSubRegions(
    BWT &bwt, const char *str, 
    std::vector<AlignType> &all_pos, std::vector<MatchType> &extended_frag
) {
  extern char DELIM;
  int num_pos = 0;
  for(auto it = all_pos.begin(); it != all_pos.end(); ++ it) {
    int diff = (long long int) (it->bwt_end) - (long long int) (it->bwt_begin);
    //cerr << "!!!!!  " << diff << " " << it->q_begin << " " << it->q_end << endl;
    //cerr << str + it->q_begin << endl;
    num_pos += diff;
  }
  SeedExtension seed_extend;
  vector<MatchType> all_match(2 * num_pos);
  
  int q_len = strlen(str);
  int idx = 0;
  //cout << "Total hits to extend:  " << num_pos << endl;
  unordered_map<BWTIDX, bool> ref_begin;
  std::stringstream buffer;
  for(auto it = all_pos.begin(); it != all_pos.end(); ++ it) {
    //cerr << "region: " << it->q_begin << " " << it->q_end << endl;
    int len = it->q_end - it->q_begin + 1;
    int idx_chunk_begin = idx; int max_score = 0;
    for(BWTIDX i = it->bwt_begin; i < it->bwt_end; ++ i) {
      //cout << "died at " << i << endl;
      
      // record the original hit before extension
      all_match[idx].str_pos = bwt.BWTToStrPos(i);    
      all_match[idx].q_begin = it->q_begin; all_match[idx].q_end = it->q_end;
      all_match[idx].str_len = len;
      all_match[idx].score = len * NT_MATCH;  // assume all seeds that correspond to perfect matches
      all_match[idx].num_mismatch = 0;
      MatchType record = all_match[idx];
      extended_frag.push_back(record);
      extended_frag[extended_frag.size() - 1].cigar = new char [CIGAR_ARRAY_SIZE];
      sprintf(extended_frag[extended_frag.size() - 1].cigar, "%d%c", all_match[idx].str_len, 'M');
      all_match[idx].cigar = new char [CIGAR_ARRAY_SIZE];
      sprintf(all_match[idx].cigar, "%d%c", all_match[idx].str_len, 'M');
      
      //cerr << "check cigar:   " << extended_frag[extended_frag.size() - 1].cigar << endl;
      
      // attempt to extend the seed ungapped
      //cerr << "Before ungapped extension:  " << all_match[idx].q_begin << "  " << all_match[idx].q_end << " " << all_match[idx].score << "  " << all_match[idx].cigar << endl;
      // ungapped extension of fragments
      BWTIDX str_end = all_match[idx].str_pos + (BWTIDX) len - 1;
      //if(all_match[idx].str_pos < 0) {
      //  cerr << "Before ExtendSeeds: " << all_match[idx].str_pos << endl;
      //  exit(1);
      //}
      // TODO: should use the reference string for the extension rather than using the entire concat sequence
      all_match[idx].score += seed_extend.UngappedExt<BWTIDX, BWTINT>(
          (const char*) bwt.str_, all_match[idx].str_pos, str_end, bwt.size_,
          str, all_match[idx].q_begin, all_match[idx].q_end, (BWTINT) q_len,
          all_match[idx].num_mismatch, true, DELIM
      );
      
      //if(all_match[idx].str_pos < 0) {
      //  cerr << "After ExtendSeeds: " << all_match[idx].str_pos << endl;
      //  exit(1);
      //}
      all_match[idx].str_len = str_end - all_match[idx].str_pos + 1;
      sprintf(all_match[idx].cigar, "%d%c", all_match[idx].str_len, 'M');
      max_score = all_match[idx].score > max_score ? all_match[idx].score : max_score;
      
      //cerr << "After ungapped extension:  " << all_match[idx].q_begin << "  " << all_match[idx].q_end << " " << all_match[idx].score << "  " << all_match[idx].cigar << endl;
      //cerr << "Ungapped Extended:  " << all_match[idx].q_begin << "  " << all_match[idx].q_end << " " << all_match[idx].score << endl;
      ++ idx;
    }
    // only record ungapped mappings that have the highest score
    for(int ii = idx_chunk_begin; ii < idx; ++ ii) {
      //cerr << "Gapped max score:  " << all_match[ii].score << " " << max_score << " " << (bool) (ref_begin.find(all_match[ii].str_pos) == ref_begin.end()) << endl;
      if(all_match[ii].score >= max_score )  {
        // performs gapped alignment extension
        BWTIDX str_end = all_match[ii].str_pos + all_match[ii].str_len - 1;
        // DEBUG: try skipping ungapped alignment
        // TODO: check this GappedExt function
        all_match[ii].score += seed_extend.GappedExt<BWTIDX, BWTINT>(
          (const char*) bwt.str_, all_match[ii].str_pos, str_end, bwt.size_,
          str, all_match[ii].q_begin, all_match[ii].q_end, (BWTINT) q_len,
          all_match[ii].cigar, all_match[ii].num_mismatch, true, DELIM
        );
        all_match[ii].str_len = str_end - all_match[ii].str_pos + 1;
        //ref_begin[all_match[ii].str_pos] = true;
        //cerr << "pushed:  " << all_match[ii].q_begin << "   " << all_match[ii].q_end << "   " << all_match[idx].score << "  " << endl;
        MatchType record = all_match[ii];
        extended_frag.push_back(record);
        extended_frag[extended_frag.size() - 1].cigar = new char [CIGAR_ARRAY_SIZE];
        sprintf(extended_frag[extended_frag.size() - 1].cigar, "%s", all_match[ii].cigar);
      }
    }  
  }
  
  for(int ii = 0; ii < idx; ++ ii) {
    delete [] all_match[ii].cigar;
  }
  return;
} 

void BWTSearch::Search(
    BWT &bwt, BWT &rev_bwt, 
    const char *str, AlignType &pos
) {
  int len = strlen(str);
  if(len <= 0)  {
    cerr << "BWTSearch::Search: Searching empty sequence. An empty list will be returned." << endl; return;
  }
  int *bound = new int [len];
  CalLowerBound(rev_bwt, str, bound);
  SearchInExact(bwt, rev_bwt, bound, str, pos);  
  AlignType prev_pos = pos;
  if(pos.q_begin >= 0)  {
      string rev_seq = str;
      rev_seq = string(rev_seq.rbegin(), rev_seq.rend());
      CalLowerBound(bwt, rev_seq.c_str(), bound);   
      SearchInExact(rev_bwt, bwt, bound, rev_seq.c_str(), pos);
  }
  int prev_span = prev_pos.q_end - prev_pos.q_begin + 1;
  int span = pos.q_end - pos.q_begin + 1;
  if(prev_span > span || (prev_span == span && prev_pos.cost < pos.cost)) { 
    pos = prev_pos;
  } else  {
    // recall that "pos" was the search result of the reverse BWT
    // need to flip the starting and ending positions for the query
    pos.q_end = len - pos.q_begin - 1; pos.q_begin = 0;
    // TODO: notice that the BWT indexes are for the reverse BWT
    // TODO: and should be converted to exact reference location before returning
  }
  delete [] bound;
  return;
}


void BWTSearch::Enqueue(
    std::priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> &candidate, 
    ExtInfo &phase_info
) {
  if(candidate.size() < max_queue_size || phase_info.cost <= candidate.top().cost)
    candidate.push(phase_info);
  return;
}

void BWTSearch::SearchInExact(
    BWT &bwt, BWT &rev_bwt, const int *bound,
    const char *str, AlignType &pos
) {
  
  int len = strlen(str);
  int lower_cost = cost;
  // construct initial stack 
  priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> candidate;
  ExtInfo init;
  init.range.first = 0; init.range.second = bwt.GetSize(); init.q_pos = len - 1; init.cost = 0;
  init.cost_bound = bound[len - 1];
  Enqueue(candidate, init);
  pos.bwt_begin = 0; pos.bwt_end = bwt.GetSize(); 
  pos.q_begin = pos.q_end = len - 1; pos.cost = cost;
  // extend each element in the stack with fast mode
  while(!candidate.empty()) {
    // record the best extension so par
    if(candidate.top().q_pos < pos.q_begin || 
        (candidate.top().q_pos == pos.q_begin && candidate.top().cost < pos.cost)
    ) { 
      pos.q_begin = candidate.top().q_pos + 1; pos.cost = candidate.top().cost;
      pos.bwt_begin = candidate.top().range.first; pos.bwt_end = candidate.top().range.second;
    }
       
    if(candidate.top().q_pos < 0)  {
      if(candidate.top().cost < lower_cost) lower_cost = candidate.top().cost;
      candidate.pop();
    } else 
      ExtendInExactFast(bwt, str, candidate, bound, lower_cost);
  }
  // if the entire sequence cannot be aligned, report the longest alignment within 
  // the allowed number of errors
  return;
}

void BWTSearch::ExtendInExactFast(
    BWT &bwt, const char *str,
    std::priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> &candidate, 
    const int *bound, int lower_cost
) {

  if(candidate.empty()) return;
  // get the top extension
  ExtInfo current = candidate.top();
  candidate.pop(); 
  if(current.q_pos < 0 || current.cost + bound[current.q_pos] > lower_cost)  return;
  // perform extension
  ExtInfo next;
  // tries to find if the next is a match
  next.range = bwt.UpdateRange(str[current.q_pos], current.range);
  if(next.range.first < next.range.second)  {
    next.q_pos = current.q_pos - 1; next.cost = current.cost;
    next.cost_bound = next.cost;
    if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
    if(next.q_pos > 0)  {
      pair<BWTIDX, BWTIDX> peek_ahead = bwt.UpdateRange(str[current.q_pos - 1], next.range);
      if(peek_ahead.first < peek_ahead.second && next.cost_bound <= lower_cost) {
        Enqueue(candidate, next); return;
      }
    } else {
      // can't peek ahead, directly insert
      if(next.cost_bound <= lower_cost) {
        Enqueue(candidate, next); return;
      }
    }
  }
  // otherwise try all possible mismatch/insertion/deletions
  for(int i = 0; i < bwt.alphabet_.alphabet_size_; ++ i) {
    BWTCHAR c = bwt.alphabet_.GetInvCharMap(i);
    // refine the range
    next.range = bwt.UpdateRange(c, current.range);
    if(next.range.first >= next.range.second) continue;    
    // the match/mismatch case; subtract the index
    if(c == str[current.q_pos]) next.cost = current.cost;
    else next.cost = current.cost + m_cost;
    next.q_pos = current.q_pos - 1;
    next.cost_bound = next.cost;
    if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
    if(next.cost_bound <= lower_cost) {
      Enqueue(candidate, next);
    }   
    // the insertion case (any character "c" is inserted into the sequence)
    // add gap cost, keep index
    next.cost = current.cost + g_cost;
    next.q_pos = current.q_pos;
    next.cost_bound = next.cost;
    if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
    if(next.cost_bound <= lower_cost) {
      Enqueue(candidate, next);
    }
  }
  // the deletion case (skip current character, add gap cost)
  next = current;
  -- next.q_pos; next.cost += g_cost;
  next.cost_bound = next.cost;
  if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
  if(next.cost_bound <= lower_cost) {
    Enqueue(candidate, next);
  }
  return;
}

void BWTSearch::CalLowerBound(BWT &rev_bwt, const char *str, int *bound) {
  int cost = 0;
  pair<BWTIDX, BWTIDX> range;
  range.first = 0; range.second = rev_bwt.GetSize();
  for(int i = 0; i < strlen(str); ++ i) {
    range = rev_bwt.UpdateRange(str[i], range);
    if(range.first >= range.second)  {
      range.first = 0; range.second = rev_bwt.GetSize();
      cost += m_cost;
    }
    bound[i] = cost;
  }
  return;
}

// given a read of interest and minimum overlap, find all intervals (in both fw and re FM-indexes) 
// corresponding to the prefix of the reads that perfectly overlap with the given read
void BWTSearch::SearchBeginIntervals(const char* seq, const int min_len, IvInfo &search_info) {
  int n = strlen(seq);
  // if the total length of the read is less than the minimum overlap, return
  if(n < min_len) return;
  // initialize the search for the overlapping region
  pair<BWTIDX, BWTIDX> fw_range, re_range, fw_range_terminal, re_range_terminal;
  fw_range.first = re_range.first = 0; 
  fw_range.second = re_range.second = search_info.bwtF_->GetSize();
  int i;
  for(i = 0; i < min_len - 1; ++ i) {
    //cout << "search sequence: " << &seq[n - i - 1] << endl;
    char c = seq[n - i - 1];
    BWTIDX occbegin = search_info.bwtF_->CountOccurrence(c, fw_range.first);
    BWTIDX occend = search_info.bwtF_->CountOccurrence(c, fw_range.second);
    //cout << "initial range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
    //cout << "search occurrence:  " << c << " " << occbegin << "  " << occend << endl;
    if(occbegin >= occend)  break;
    //cout << "lexicographical smaller: " << search_info.bwtF_->CountLexicoLess(c, fw_range.first) << "  " << search_info.bwtF_->CountLexicoLess(c, fw_range.second) << endl;
    re_range.first = re_range.first
        + search_info.bwtF_->CountLexicoLess(c, fw_range.second)
        - search_info.bwtF_->CountLexicoLess(c, fw_range.first); 
    re_range.second = re_range.first + occend - occbegin;
    int c_id = search_info.bwtF_->alphabet_.GetCharMap(c);
    fw_range.first = search_info.bwtF_->acc_freq_[c_id + 1] + occbegin;
    fw_range.second = search_info.bwtF_->acc_freq_[c_id + 1] + occend;
    //cout << "updated range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
  }
  // if no sequence overlap for the given minumum overlap length, return
  if(i < min_len - 1) return;
  // otherwise also try to search the delimitor to detect begin intervals
  // recall that i is less than n - 1 because we do not want the read itself
  // or other wise the read is contained by other reads
  for(i = min_len - 1; i < n - 1; ++ i) {
    // extend the sequence
    //cout << "search sequence: " << &seq[n - i - 1] << " " << n << " " << i << endl;
    //cout << "initial range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
    char c = seq[n - i - 1];
    BWTIDX occbegin = search_info.bwtF_->CountOccurrence(c, fw_range.first);
    BWTIDX occend = search_info.bwtF_->CountOccurrence(c, fw_range.second);
    //cout << "search occurrence:  " << c << " " << occbegin << "  " << occend << endl;
    //cout << "lexicographical smaller: " << search_info.bwtF_->CountLexicoLess(c, fw_range.first) << "  " << search_info.bwtF_->CountLexicoLess(c, fw_range.second) << endl;
    if(occbegin >= occend)  break;
    re_range.first = re_range.first
        + search_info.bwtF_->CountLexicoLess(c, fw_range.second)
        - search_info.bwtF_->CountLexicoLess(c, fw_range.first);
    re_range.second = re_range.first + occend - occbegin;
    int c_id = search_info.bwtF_->alphabet_.GetCharMap(c);
    fw_range.first = search_info.bwtF_->acc_freq_[c_id + 1] + occbegin;
    fw_range.second = search_info.bwtF_->acc_freq_[c_id + 1] + occend;
    // also search for the delimitor (with the updated fw_range and re_range)
    //cout << "initial terminal range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
    occbegin = search_info.bwtF_->CountOccurrence(DELIM, fw_range.first);
    occend = search_info.bwtF_->CountOccurrence(DELIM, fw_range.second);
    //cout << "search terminal occurrence:  " << occbegin << "  " << occend << endl;
    if(occbegin >= occend)  continue;
    // note that delimitor is the lexico-smallest char, 
    // no need to add CountLexicoLess nor acc_freq
    re_range_terminal.first = re_range.first; 
    re_range_terminal.second = re_range.first + occend - occbegin;
    fw_range_terminal.first = occbegin;
    fw_range_terminal.second = occend;
    // record such interval
    search_info.intervals_.PushA(fw_range_terminal.first, fw_range_terminal.second);
    search_info.intervals_.PushB(re_range_terminal.first, re_range_terminal.second);
    search_info.intervals_.PushLen(i + 1);
    
    //cout << "overlap recorded !!!" << endl;
  }
  return;
}


void BWTSearch::FindIrreducible(
    IvInfo &search_info, std::vector<BWTIDX> &ir_positions, std::vector<int> &ir_overlap
) {
  
  search_info.intervals_.Reverse();  
  // check boundary conditions
  if(!search_info.intervals_.Check()) {
    cout << "Warning: corrupted intervals, no irreducible edges can be detected." << endl;
    return;
  }
  if(search_info.intervals_.GetSize() <= 0)  return;
  // the stack contains all intervals that ends with the same sequences
  stack<IvSet> candidates;
  candidates.push(search_info.intervals_);
  // recursively check each intervals
  bool used_char[256];  
  while(!candidates.empty()) {
    //cout << "============ handling each interval group ===============  " << candidates.size() << endl;
    IvSet current = candidates.top(); candidates.pop();
    //for(int h = 0; h < current.GetSize(); ++ h) {
    //  cout << "search_info " << h << ": " << current.len_[h] << "  " << current.ivA_[h].first << " " << current.ivA_[h].second << "  " << current.ivB_[h].first << " " << current.ivB_[h].second << endl;
    //}
    // check all presented characters in the alphabet
    // note that for the first time we do not check "$"-extension because it would
    // mean that the extension read is contained
    int i, j, k, n = current.GetSize();
    memset(used_char, 0, 256);
    for(i = 0; i < n; ++ i) { // for each interval
      // for each character in the alphabet (look at the reverse BWT)
      for(k = current.ivB_[i].first; k < current.ivB_[i].second; ++ k) { 
        //cout << "Interval info: " << k << ": " <<  current.len_[i] << "  " << current.ivA_[i].first << " " << current.ivA_[i].second << "  " << current.ivB_[i].first << " " << current.ivB_[i].second << endl;
        char c = (char) search_info.bwtR_->bwt_[k]; // the kth char in the reverse BWT string
        // check if the character has been tested
        if(c == DELIM || used_char[c])  continue;
        used_char[c] = true;
        IvSet next;
        // try to update the intervals by appending such character
        for(j = i; j < n; ++ j) {
          // try to append the character
          pair<BWTIDX, BWTIDX> r1 = search_info.bwtR_->UpdateRange(c, current.ivB_[j]);
          // try to check if the read ends after appending the character
          pair<BWTIDX, BWTIDX> r2 = search_info.bwtR_->UpdateRange(DELIM, r1); 
          //cout << "phase: " << j << ": " << r1.first << " " << r1.second << " " << r2.first << "  " << r2.second << endl; 
          // if the first read leads to a termination, we found an irreducible read
          // in cases where multiple reads end at the same time, take the first position
          // terminate current loop
          if(j == i && r1.second - r1.first == r2.second - r2.first)  {
            //cout << "irreducible read found!!!" << endl;
            ir_positions.push_back(r2.first); 
            ir_overlap.push_back(current.len_[j]);
            break;
          }
          // for other reads that do not terminate, add to interval set next
          if(r1.second - r1.first >= 1 && r2.second - r2.first <= 0)  {
            // we only need to record intervals at the reverse BWT
            //cout << "read extension recorded!!!" << endl;
            next.PushA(-1, -1); next.PushB(r1.first, r1.second); next.PushLen(current.len_[j]);
          }
        }
        if(next.GetSize() > 0)  candidates.push(next);
      }
    }
  }
  return;
}

bool BWTSearch::IsContainedRead(const char* seq, BWT &bwt, AlignType &pos)  {
  // constructing the extended string
  int n = strlen(seq);
  char *extended_seq = new char [n + 3];
  extended_seq[0] = DELIM;
  strcpy(&extended_seq[1], seq);
  extended_seq[n + 1] = DELIM; extended_seq[n + 2] = '\0';
  // searching the entire sequence with and without delimitor against the BWT
  SearchExact(bwt, seq, pos);
  int r1 = pos.bwt_end - pos.bwt_begin;
  SearchExact(bwt, extended_seq, pos);
  int r2 = pos.bwt_end - pos.bwt_begin;
  delete [] extended_seq;
  if(r2 == r1)  return false;
  return true;
}

