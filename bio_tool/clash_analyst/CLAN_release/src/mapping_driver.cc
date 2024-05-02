# include "mapping_driver.h"

using namespace std;

void MappingDriver::ReadMapBatch(
    char **ref_header, char **ref_seq, int num_refs,
    char **ref_header_aux, char **ref_seq_aux, int num_refs_aux,
    char **header, char **seq, int num_reads,          // the header and the sequence of the read set
    const int begin, const int size, const int upper,  // the begin, the number of sequences, and the upper bound of sequences to map in this batch
    BWT &bwt, BWT &bwt_aux,                            // the indexed BWT
    std::vector< std::vector<FragType> > &result,      // the output 
    const CLANSearchParam &param_set
) {
  assert(size > 0);
  //cerr << "begin mapping" << endl;
////////////////////////////////////////
  //cerr << "num refs: " << num_refs << endl;
  //BWTIDX *ref_len = new BWTIDX [num_refs];
  //for(int k = 0; k < num_refs; ++ k) {
  //  string s = ref_seq[k];
  //  ref_len[k] = (BWTIDX) s.length();
  //}
  //delete [] ref_len;
  //return;
////////////////////////////////////////
  
  
  for(int i = 0; i < size && begin + i < upper; ++ i) {
    //cerr << "Working on " << header[begin + i] << endl;
    //cerr << "seq  " << seq[begin + i] << endl;
    BWTSearch searcher;
    BWTFragMerger merger;

    vector<AlignType> all_hits, all_hits_aux;
    vector<MatchType> merged_frag, extended_frag, merged_frag_aux, extended_frag_aux;
    string s = seq[begin + i];
      
    // search against the main database
    searcher.SearchAllSubRegions(bwt, param_set.min_frag_len, param_set.max_hits, param_set.flank_len, s.c_str(), all_hits); 
    searcher.ExtendAllSubRegions(bwt, s.c_str(), all_hits, extended_frag);
    
    // DEBUG
    //exit(0);
    //for(int i = 0; i < all_hits.size(); ++ i) {
    //  cerr << "HIT ENTRY OUT: " << all_hits[i].q_begin << "    " << all_hits[i].q_end << endl;
    //}
    merger.MergeFragments(bwt, s.length(), extended_frag, merged_frag);
    for(auto it = merged_frag.begin(); it != merged_frag.end(); ++ it) {
      //cerr << "check fragments before merge:    " << it->q_begin << "   " << it->q_end << endl; 
      it->strand = true;
      it->db = true;
    }   
    
    
    
    // search against the auxiliary database, if necessary
    bool ensure_partition = false;
    if(num_refs_aux > 0 && ref_header_aux != NULL && ref_seq_aux != NULL)  {
      ensure_partition = true;
      searcher.SearchAllSubRegions(bwt_aux, param_set.min_frag_len, param_set.max_hits, param_set.flank_len, s.c_str(), all_hits_aux); 
      searcher.ExtendAllSubRegions(bwt_aux, s.c_str(), all_hits_aux, extended_frag_aux);
      merger.MergeFragments(bwt, s.length(), extended_frag_aux, merged_frag_aux);
      for(auto it = merged_frag_aux.begin(); it != merged_frag_aux.end(); ++ it) {
        it->strand = true;
        it->db = false;
      } 
    }
    merged_frag.insert(merged_frag.end(), merged_frag_aux.begin(), merged_frag_aux.end());
    
    int covered = merger.FindBestMatching(
        s, bwt, bwt_aux, param_set.num_frag, merged_frag, result[i], ensure_partition, 
        param_set.use_insert_penalty, param_set.frag_penalty, param_set.max_overlap, param_set.min_map_portion
    );
    
    // TODO: recover the reverse mapping mode
    //continue;
      
    //cerr << s << endl;
    //cerr << "reverse search begins" << endl;
    // if the data is not strand-specific, map to both strands 
    /*
    if(!param_set.strand_specific)  {
      string s_rv = string(s.rbegin(), s.rend());
      // search reverse direction
      for(int j = 0; j < s_rv.length(); ++ j) {
        switch(s_rv[j]) {
          case 'A': s_rv[j] = 'T'; break;
          case 'C': s_rv[j] = 'G'; break;
          case 'G': s_rv[j] = 'C'; break;
          case 'N': s_rv[j] = 'N'; break;
          case 'T': s_rv[j] = 'A'; break;
          default: break;
        } 
      }
      //cerr << s_rv << endl;
      searcher.SearchAllSubRegions(bwt, param_set.min_frag_len, param_set.max_hits, param_set.flank_len, s_rv.c_str(), all_hits_rv);
      searcher.ExtendAllSubRegions(bwt, s_rv.c_str(), all_hits_rv, extended_frag_rv);
      merger.MergeFragments(bwt, s_rv.length(), extended_frag_rv, merged_frag_rv);
      // add the newly detected fragments into the candidate pool
      for(auto it = merged_frag_rv.begin(); it != merged_frag_rv.end(); ++ it) {
        it->strand = false;
        int b = s.length() - it->q_end - 1;
        int e = s.length() - it->q_begin - 1;
        it->q_begin = b; it->q_end = e;
        merged_frag.push_back(*it);
      } 
    }
    //cerr << "before best matching" << endl;
    //continue;
    int covered = merger.FindBestMatching(
        bwt, param_set.num_frag, merged_frag, result[i], 
        param_set.use_insert_penalty, param_set.frag_penalty, param_set.max_overlap
    );
    */
    
  }
  
  
  return;
}


