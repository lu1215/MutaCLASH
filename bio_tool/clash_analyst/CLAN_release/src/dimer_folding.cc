#include "dimer_folding.h"

using namespace std;

void DimerFolding::FoldDimerHelixEnergy(
  std::string &seq,     // the sequence to be folded 
  const int terminal,   // the connection terminal of the two RNA strands, no intramolecular pair is predicted
  DuplexInfo &dimer_info,
  const int max_bulge_size
) {
  assert(max_bulge_size <= 15); // because we only have turner energy model parameter setup up to 30

  // convert 'T' to 'U' and check sequence validity
  string rna_seq(seq.length(), 'N');
  if(!FormatRNASeq(seq, rna_seq)) return;
  
  int n = rna_seq.length();
  // use dynamic programming to predict dimer MFE
  int **mbp = new int* [n];
  int **opt_arrayP = new int* [n];
  int **opt_arrayQ = new int* [n];
  for(int i = 0; i < n; ++ i) {
    mbp[i] = new int [n];
    opt_arrayP[i] = new int [n];
    opt_arrayQ[i] = new int [n];
    for(int j = 0; j < n; ++ j) {
      mbp[i][j] = 999;
      opt_arrayP[i][j] = opt_arrayQ[i][j] = -1;
    }
  }
  
  // mbp holds the best score when i and j pairs
  dimer_info.dimer_energy = 0;
  int all_opt_p = -1, all_opt_q = -1;
  for(int i = 1; i < n; ++ i) {
    for(int j = 0; j < n - i; ++ j) {
      int p = j, q = j + i; // the index of the DP table that is being processed
      if(q < terminal) continue;  // this pair is intramolecular, skip
      string bp = rna_seq.substr(p, 1) + rna_seq.substr(q, 1);
      // if it cannot form base pair, continue
      if(free_energy::basepair_map.find(bp) == free_energy::basepair_map.end()) continue;
      // if it can form a base pair
      int cmfe = 0; // a single pair does not have energy (nearest neighbor model, but need to overwrite the MAX energy value 999)
      int opt_p = -1, opt_q = -1;
      // now search if lower energy can be reached
      
      // check stack
      if(mbp[p + 1][q - 1] < 999)  {
        string prev_bp = rna_seq.substr(p + 1, 1) + rna_seq.substr(q - 1, 1);
        if(free_energy::basepair_map.find(prev_bp) != free_energy::basepair_map.end())  {
          int cfe = mbp[p + 1][q - 1] 
            + free_energy::stack[free_energy::basepair_map[prev_bp]][free_energy::basepair_map[bp]];
          //if(cfe < -10000) cerr << "!!! stack " << p << " " << q << " " << p + 1 << " " << q - 1 << endl;
          if(cfe <  cmfe) {
            cmfe = cfe; opt_p = p + 1; opt_q = q - 1;
          }
        }
      }
      
      // check internal loop and bulge loop 
      for(int r = p + 1; r < q - 1 && r - p <= max_bulge_size; ++ r) {
        for(int s = q - 1; s > r && q - s <= max_bulge_size; -- s) {
          // r and s is the next closing pair
          if(mbp[r][s] >= 999) continue; 
          if(r - p == 1 || q - s == 1 || (r - p + q - s) <= 5)  { // bulge case
            int cfe = mbp[r][s] + free_energy::bulge[r - p + q - s - 2];
            //if(cfe < -10000) {
            //  cerr << "!!! bulge loop " << p << " " << q << " " << r << " " << s << endl;
            //  cerr << " prev sol: " << mbp[r][s] << endl;
            //  cerr << " added energy: " << free_energy::bulge[r - p + q - s - 2] << " " << r - p + q - s - 2 << endl;
            //}
            if(cfe < cmfe) {
              cmfe = cfe; opt_p = r; opt_q = s;
            }
          } else  { // internal loop case
            int cfe = mbp[r][s] + free_energy::internal[r - p + q - s - 2]
              + free_energy::internal_asymmetric * abs((r - p) - (q - s));
            // special closure energy
            if(bp == "AU" || bp == "UA" || bp == "GU" || bp == "UG")  cfe += free_energy::internal_closure_AUGU;
            //if(cfe < -10000) cerr << "!!! internal loop " << p << " " << q << " " << r << " " << s << endl;
            if(cfe < cmfe) {
              cmfe = cfe; opt_p = r; opt_q = s;
            }
          }
        }
      }
      
      // update the info
      mbp[p][q] = cmfe;
      opt_arrayP[p][q] = opt_p;
      opt_arrayQ[p][q] = opt_q;
      
      if(cmfe < dimer_info.dimer_energy)  {
        dimer_info.dimer_energy = cmfe; all_opt_p = opt_p; all_opt_q = opt_q;
      }
      
    }
  }
  
  /*
  for(int i = 0; i < n; ++ i) {
    for(int j = 0; j < n; ++ j) {
      cerr << mbp[i][j] << "\t";
    }
    cerr << endl;
  }
  */
  
  // traceback and find the optimal solution
  dimer_info.dimer_pair.resize(rna_seq.length());
  dimer_info.stack_info.energy = 0; dimer_info.stack_info.begin = dimer_info.stack_info.end = -1;
  int stack_p_begin = all_opt_p, stack_q_end = all_opt_q;
  for(int i = 0; i < dimer_info.dimer_pair.size(); ++ i) dimer_info.dimer_pair[i] = -1;
  while(all_opt_p >= 0 && all_opt_q >= 0) {
    dimer_info.dimer_pair[all_opt_p] = all_opt_q;
    dimer_info.dimer_pair[all_opt_q] = all_opt_p;
    int cp = opt_arrayP[all_opt_p][all_opt_q];
    int cq = opt_arrayQ[all_opt_p][all_opt_q];
    // check if is a perfect stack
    //cerr << all_opt_p << "\t" << all_opt_q << "\t" << cp << "\t" << cq << endl;
    if(cp != all_opt_p + 1 || cq != all_opt_q - 1)  {
      //cerr << "termination of a stak" << endl;
      float dG = mbp[stack_p_begin][stack_q_end] - mbp[all_opt_p][all_opt_q];
      if(dG < dimer_info.stack_info.energy) {
        //cerr << "updated information: " << dG << "  " << dimer_info.stack_info.energy << endl;
        dimer_info.stack_info.energy = dG; 
        dimer_info.stack_info.begin = stack_p_begin; dimer_info.stack_info.end = all_opt_p;
      }
      stack_p_begin = cp; stack_q_end = cq;
    }
    // update the index   
    all_opt_p = cp; all_opt_q = cq;
  }
  // check the last stack
  
  dimer_info.dimer_energy /= SCALE_TO_INT;
  dimer_info.stack_info.energy /= SCALE_TO_INT;
  
  // collect memory
  for(int i = 0; i < n; ++ i) {
    delete [] mbp[i]; delete [] opt_arrayP[i]; delete [] opt_arrayQ[i];
  }
  delete [] mbp;
  delete [] opt_arrayP; delete [] opt_arrayQ;
  
  return;
}
