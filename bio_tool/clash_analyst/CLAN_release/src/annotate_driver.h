#ifndef _ANNOTATE_DRIVER_H_
#define _ANNOTATE_DRIVER_H_

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

#include "parameters.h"
#include "read_merge.h"
#include "free_energy.h"
#include "output_printer.h"
#include "misc.h"

#define MAX_SEQ_BUFFER_SIZE 1000

class AnnotateDriver  {
 public:
 
  AnnotateDriver()  {}
  ~AnnotateDriver() {}
 
  void AnnotateMapping(
    char const* const* ref_header, char const* const* ref_seq, const int num_ref, 
    char const* const* header, char const* const* seq, const int num_read,  
    const int ave_read_len,
    std::vector< std::vector<FragType> > *mapping,
    const CLANAnnotateParam &param_set
  )  {
    std::ofstream out_fh;
    out_fh.open(param_set.f_output, std::ios::out);
    if(!out_fh.is_open())  {
      std::cerr << "CLAN:clan_annotate: unable to write to " << param_set.f_output << "; Abort." << std::endl;
      exit(1);
    }
    //std::cerr << "annotation begins!!!" << std::endl;
    // get length of the refs
    std::vector<BWTIDX> ref_len; ref_len.resize(num_ref);
    for(int i = 0; i < num_ref; ++ i) {
      ref_len[i] = std::string(ref_seq[i]).length();
    }
  
    // annotate the mapping
    std::vector<GPosType> *islands = new std::vector<GPosType>;
    std::vector< std::unordered_map<int, float> > *duplex = new std::vector< std::unordered_map<int, float> >;
    
    ReadMerge *read_merger = new ReadMerge;
    read_merger->MergeReads(num_ref, ref_header, ref_seq, ave_read_len, *mapping, *islands, *duplex, param_set);
    delete read_merger;
    
    if(param_set.print_simplex)  {
      std::vector<DuplexType> foo; 
      OutputPrinter output_printer;
      output_printer.OutputAnnotationTable(param_set.print_simplex, ref_header, *islands, foo, out_fh);
      out_fh.close();
      delete islands;
      delete duplex;
      return;
    }
    
    //std::cerr << "Duplex detection finished" << std::endl;
    int max_island_len = 0;
    for(int i = 0; i < islands->size(); ++ i) {
      if(max_island_len < (*islands)[i].end - (*islands)[i].begin + 1) {
        max_island_len = (*islands)[i].end - (*islands)[i].begin + 1;
      }
    }
    std::cerr << "Longest island: " << max_island_len << std::endl;
    std::cerr << "Num. of island: " << islands->size() << std::endl;
    
    
    // prepare sequence for folding
    std::vector<std::string> *island_seq = new std::vector<std::string>;
    island_seq->resize(islands->size());
    for(int i = 0; i < islands->size(); ++ i) {
      BWTIDX begin = (*islands)[i].begin - param_set.extend_len;
      BWTIDX end = (*islands)[i].end + param_set.extend_len;
      begin = begin < 0 ? 0 : begin;
      end = end < ref_len[(*islands)[i].sid] ? end : ref_len[(*islands)[i].sid] - 1;
      //std::cerr << begin << "\t" << end << "\t" << std::string(ref_seq[islands[i].sid]).length() << "\t" << i << "\t" << islands.size() << std::endl;
      BWTIDX len = end - begin + 1;
      char *tstr = new char [len + 1];
      memcpy(tstr, ref_seq[(*islands)[i].sid] + begin, sizeof(char) * len);
      tstr[len] = '\0';
      (*island_seq)[i] = std::string(tstr); 
      delete [] tstr;
    }

    // perform dimer folding
    int num_duplex = 0;
    for(int i = 0; i < duplex->size(); ++ i) {
      for(auto it = (*duplex)[i].begin(); it != (*duplex)[i].end(); ++ it) {
        ++ num_duplex;  
      }
    }
    
    std::cerr << "Island seqeunce preparation done: " << std::endl;
    std::cerr << "Num duplex: " << num_duplex << std::endl;
    
    std::vector<DuplexType> *duplex_annot = new std::vector<DuplexType>; 
    //duplex_annot->resize(num_duplex);
    DimerFolding dimer_folder;
    int idx = 0;

    
    #pragma omp parallel num_threads(param_set.num_threads)
    {
      #pragma omp for
      for(int i = 0; i < duplex->size(); ++ i) {
        for(auto it = (*duplex)[i].begin(); it != (*duplex)[i].end(); ++ it) {
          DuplexType di; di.island_ID1 = i; di.island_ID2 = it->first;
          di.strand_1 = (*islands)[i].strand; di.strand_2 = (*islands)[it->first].strand;
          if(param_set.strand_specific)  {
            std::string dimer_concat = (*island_seq)[i] + (*island_seq)[it->first];
            
            //std::cerr << "s1: " << (*island_seq)[i] << "  " << std::endl;
            //std::cerr << "s2: " << (*island_seq)[it->first] << "  " << std::endl;
            //std::cerr << "s1+s2: " << dimer_concat << std::endl;
            
            dimer_folder.FoldDimerHelixEnergy(
              dimer_concat, (*island_seq)[i].length(), 
              di.info
            );      
          } else  {
            bool sr_first, sr_second;
            float min_mfe = 999999;
            for(int m = 0; m < 4; ++ m) {
              sr_first = (bool) floor(m / 2);
              sr_second = (bool) (m % 2 == 0);
              std::string s1 = sr_first ? (*island_seq)[i] : misc::RevComplement((*island_seq)[i]);
              std::string s2 = sr_second ? (*island_seq)[it->first] : misc::RevComplement((*island_seq)[it->first]);
              
              std::string dimer_concat = s1 + s2;
              
              //std::cerr << "s1: " << (*island_seq)[i] << "  " << sr_first << std::endl;
              //std::cerr << "s2: " << (*island_seq)[it->first] << "  " << sr_second << std::endl;
              //std::cerr << "s1+s2: " << dimer_concat << std::endl;
              
              DuplexInfo dinfo; // note that the coverage is not initailized
              dimer_folder.FoldDimerHelixEnergy(
                dimer_concat, s1.length(), dinfo
              );
              if(dinfo.dimer_energy < min_mfe) {
                min_mfe = dinfo.dimer_energy;
                di.info = dinfo;
                // if the island strand is positive, and we take the positive strand, then the resulting strand is positive 
                // ... positive, negative, negative
                // ... negative, positive, negative
                // ... negative, negative, positive
                di.strand_1 = ((*islands)[i].strand == sr_first);
                di.strand_2 = ((*islands)[it->first].strand == sr_second);
              }  
            }
          }
          di.info.coverage = it->second;
          #pragma omp critical
          {
            if(di.info.stack_info.end - di.info.stack_info.begin + 1 >= param_set.min_seed_len &&
               di.info.dimer_energy <= param_set.max_dimer_dG
            )  {
              duplex_annot->push_back(di);
            }
          }
        }
      }
    }
    //std::cerr << "Folding finished" << std::endl;
    
    
    OutputPrinter output_printer;
    output_printer.OutputAnnotationTable(param_set.print_simplex, ref_header, *islands, *duplex_annot, out_fh);
    out_fh.close();
    
    delete island_seq;
    delete islands;
    delete duplex;
    delete duplex_annot;
    //std::cerr << "Output finished" << std::endl;
    
    /*
    std::string s1 = "AAACAACAGAACGAGCACTGGACTTGGAGCCAGAAGTCTTGGGCTCAAGCCC";
    std::string s2 = "TGGGGCCTCTGCTCAAGGCGGGAATCTCGGGTGTCTACACAGAGTCA";
    std::string s_concat = s1 + s2;
    
    DimerFolding dimer_folder;
    DuplexInfo duplex_fold;
    dimer_folder.FoldDimerHelixEnergy(s_concat, s1.length(), duplex_fold);
    std::cerr << duplex_fold.dimer_energy << std::endl;
    std::cerr << duplex_fold.stack_info.energy << "  " << duplex_fold.stack_info.begin << " " << duplex_fold.stack_info.end << std::endl;
    
    for(int k = 0; k < duplex_fold.dimer_pair.size(); ++ k) {
      if(duplex_fold.dimer_pair[k] >= 0 && k < duplex_fold.dimer_pair[k])  
        std::cerr << "(";
      else if(duplex_fold.dimer_pair[k] >= 0 && k > duplex_fold.dimer_pair[k])
        std::cerr << ")";
      else
        std::cerr << ".";
    }
    std::cerr << std::endl;
    for(int k = 0; k < duplex_fold.dimer_pair.size(); ++ k) {
      if(k >= duplex_fold.stack_info.begin && k <= duplex_fold.stack_info.end)  
        std::cerr << "(";
      else if(duplex_fold.dimer_pair[k] >= duplex_fold.stack_info.begin && duplex_fold.dimer_pair[k] <= duplex_fold.stack_info.end)
        std::cerr << ")";
      else
        std::cerr << ".";
    }
    std::cerr << std::endl;
    */
    return;
  } 
  
  void CleanMappingCIGAR(std::vector< std::vector<FragType> > &mapping)  {
    for(int i = 0; i < mapping.size(); ++ i) {
      for(int j = 0; j < mapping[i].size(); ++ j) {
        for(int k = 0; k < mapping[i][j].targets.size(); ++ k) {
          delete [] mapping[i][j].targets[k].cigar;
        }
      }
    }
    return;
  }
  
};

#endif

