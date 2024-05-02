#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

struct  CLANSearchParam {
  int num_threads;
  int max_hits;
  int min_frag_len;
  int flank_len;
  int num_frag;
  int frag_penalty;
  int max_overlap;
  float min_map_portion;
  bool strand_specific;
  bool use_insert_penalty;
  char f_reference[1000]; 
  char f_index[1000]; 
  char F_reference[1000];
  char F_index[1000];
  char f_input[1000]; 
  char f_output[1000]; 
};

struct CLANIndexParam {
  char f_reference[1000];
  char f_index[1000];
};

struct CLANOutputParam  {
  int mode;
  char f_reference[1000]; 
  char F_reference[1000]; 
  char f_reads[1000];
  char f_input[1000]; 
  char f_output[1000]; 
};

struct CLANAnnotateParam {
  char f_input[1000];
  char f_output[1000];
  char f_reference[1000]; 
  char f_reads[1000];
  bool strand_specific;
  bool print_simplex;
  int extend_len;
  int num_threads;
  int min_island_len;
  int max_island_len;
  int max_island_gap;
  int min_coverage;
  float max_dimer_dG;
  int min_seed_len;
};

#endif
