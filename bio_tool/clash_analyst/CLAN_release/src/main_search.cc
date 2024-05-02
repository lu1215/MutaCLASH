#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <ctime>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <omp.h>

#include "bwt.h"
#include "bwt_search.h"
#include "loader.h"
#include "bio_alphabet.h"
#include "bwt_frag_merger.h"
#include "mapping_driver.h"
#include "output_printer.h"
#include "parameters.h"
#include "misc.h"

#ifndef NUM_FRAGS
  #define NUM_FRAGS 2
#endif

#ifndef MAP_CHUNK_SIZE
  #define MAP_CHUNK_SIZE 1000
#endif

using namespace std;


//int num_threads = 1;
//int max_hits = 20;
//int min_frag_len = 10;
//int flank_len = 4;
//int num_frag = 2;
//int frag_penalty = 10;
//int max_overlap = 5;
//bool strand_specific = false;
//bool use_insert_penalty = true;
//char *f_reference = new char [1000]; 
//char *f_index = new char [1000]; 
//char *f_input = new char [1000]; 
//char *f_output = new char [1000]; 

static CLANSearchParam param_set;
int main(int argc, char **argv)  {

  // setting default parameters
  param_set.num_threads = 1;
  param_set.max_hits = 20;
  param_set.min_frag_len = 10;
  param_set.flank_len = 4;
  param_set.num_frag = 2;
  param_set.frag_penalty = 3;
  param_set.max_overlap = 2;
  param_set.min_map_portion = 0.6;
  param_set.strand_specific = true;
  param_set.use_insert_penalty = true;
  // the following directories correspond to the first database (mandatory)
  param_set.f_reference[0] = '\0';
  param_set.f_index[0] = '\0';
  // the following directories correspond to the second database (optional)
  param_set.F_reference[0] = '\0';
  param_set.F_index[0] = '\0';
  param_set.f_input[0] = '\0';
  param_set.f_output[0] = '\0';

  int copt;	
	extern char *optarg;
  extern int optind;
  while ((copt=getopt(argc,argv,"r:o:f:F:d:D:t:m:n:k:l:p:v:c:seh")) != EOF)	{
    switch(copt) {
      case 'r':
        sscanf(optarg, "%s", param_set.f_input);
        continue;
      case 'o':
        sscanf(optarg, "%s", param_set.f_output);
        continue;
      case 'f':
        sscanf(optarg, "%s", param_set.f_reference);
        continue;
      case 'd':
        sscanf(optarg, "%s", param_set.f_index);
        continue;
      case 'F':
        sscanf(optarg, "%s", param_set.F_reference);
        continue;
      case 'D':
        sscanf(optarg, "%s", param_set.F_index);
        continue;
      case 't':
        sscanf(optarg, "%d", &param_set.num_threads);
        continue;
      case 'm':
        sscanf(optarg, "%d", &param_set.max_hits);
        continue;
      case 'k':
        sscanf(optarg, "%d", &param_set.flank_len);
        continue;
      case 'l':
        sscanf(optarg, "%d", &param_set.min_frag_len);
        continue;
      case 'p':
        sscanf(optarg, "%d", &param_set.frag_penalty);
        continue;
      case 'v':
        sscanf(optarg, "%d", &param_set.max_overlap);
        continue;
      case 'c':
        sscanf(optarg, "%f", &param_set.min_map_portion);
        continue;
      case 's':
        param_set.strand_specific = false;
        continue;
      case 'e':
        param_set.use_insert_penalty = false;
        continue;
			case 'h':
			default:
        cout << "==========================================================" << endl;
        cout << "\tCLAN: the CrossLinked reads ANalaysis tool" << endl;
        cout << "==========================================================" << endl;
        cout << endl;
        cout << "usage: clan_search -r [READ_FILE] -o [OUTPUT_FILE] -f [REFERENCE_FILE] -d [INDEX_PREFIX]" << endl;
        cout << endl;
        cout << "\tr: the file containing the reads, in FASTA format (mandatory)" << endl;
        cout << "\to: the file for writing the temporary mapping results (mandatory)" << endl;
        cout << "\tf: the reference database (e.g. the human genome, mandatory)" << endl;
        cout << "\td: the prefix of the indexes (prefix of the database index constructed by \"clan_index\", mandatory)" << endl;
        cout << "\tF: the second reference database (e.g. if the reference can be partitioned into two databases such as 3UTR and miRNA, optional)" << endl;
        cout << "\tD: the prefix of the indexes (prefix of the second database index constructed by \"clan_index\", optional)" << endl;
        cout << "\ts: map to both strands of the reference databases (default FALSE)" << endl;
        cout << "\te: disable insert penalty (penalize insert sequence between two arms, default TRUE)" << endl;
        cout << "\tt: number of threads to use (optional, default 1)" << endl;
        cout << "\tm: number of maximum hits for each maximal fragment (optional, default 20)" << endl;
        cout << "\tk: length of flanking bases when recording non-maximal fragments (optional, default 4)" << endl;
        cout << "\tl: minimum length for each fragment (optional, default 10)" << endl;
        cout << "\tp: penalty for introducing one extra strand (optional, default 3)" << endl;
        cout << "\tv: maximum overlap allowed between mapped regions (optional, default 2)" << endl;
        cout << "\tc: minimum portion of the read that needs to be mapped (optional, default 0.6)" << endl;
        cout << "\th: print this help message" << endl << endl;
        exit(0);
		}
		optind--;
	}	
  
  if(strlen(param_set.f_reference) <= 0 || strlen(param_set.f_index) <= 0 || 
      strlen(param_set.f_input) <= 0 || strlen(param_set.f_output) <= 0)  
  {
    cerr << "Mandatory argument missing; please type \"clan_search -h\" to view the help information." << endl;
    cerr << "Abort." << endl;
    exit(1);
  }
  if((strlen(param_set.F_reference) > 0 && strlen(param_set.F_index) == 0) || (strlen(param_set.F_reference) == 0 && strlen(param_set.F_index) > 0))  
  {
    cerr << "Information of auxiliary database is partially set (eithe sequence or index is missing);" << endl;
    cerr << "please type \"clan_search -h\" to view the help information." << endl;
    cerr << "Abort." << endl;
    exit(1);
  }
  
  ofstream out_fh;
  out_fh.open(param_set.f_output, ios::out);
  if(!out_fh.is_open())  {
    cerr << "clan_search: unable to write to " << param_set.f_output << "; Abort." << endl;
    exit(1);
  }

  double start_time = misc::MyTime();
  double check_time;
  BioAlphabet dna_alphabet(DNA);
  Loader loader;

  // construct index for the first database  
  int num_refs = loader.CountFastaNumSeqs(param_set.f_reference);
  char **ref_header = new char* [num_refs];
  char **ref_seq = new char* [num_refs];
  num_refs = loader.LoadFasta(dna_alphabet, param_set.f_reference, ref_header, ref_seq);
  string concat_seq; 
  Concatenator concat_obj(ref_seq, num_refs, concat_seq);
  BWT bwt;
  bwt.ConstructFromIndex(
      dna_alphabet, concat_seq.c_str(), param_set.f_index
  );
  //cout << "Done constructing object" << endl;
  bwt.ConstructLocationInfo(num_refs, ref_header, ref_seq);
  cerr << "CLAN: Finish loading main sequence database." << endl; 
  
  // construct index for the second database
  int num_refs_aux = 0;
  char **ref_header_aux = NULL;
  char **ref_seq_aux = NULL;
  string concat_seq_aux;
  BWT bwt_aux;
  if(strlen(param_set.F_reference) > 0 && strlen(param_set.F_index) > 0)  {
    num_refs_aux = loader.CountFastaNumSeqs(param_set.F_reference);
    ref_header_aux = new char* [num_refs_aux];
    ref_seq_aux = new char* [num_refs_aux];
    num_refs_aux = loader.LoadFasta(dna_alphabet, param_set.F_reference, ref_header_aux, ref_seq_aux);
    Concatenator concat_obj(ref_seq_aux, num_refs_aux, concat_seq_aux);
    bwt_aux.ConstructFromIndex(
        dna_alphabet, concat_seq_aux.c_str(), param_set.F_index
    );
    //cout << "Done constructing object" << endl;
    bwt_aux.ConstructLocationInfo(num_refs_aux, ref_header_aux, ref_seq_aux);
    cerr << "CLAN: Finish loading auxiliary sequence database" << endl;
  }  
      
  // load sequece data
  int num_reads;
  if(loader.IsFASTA(param_set.f_input)) {
    num_reads = loader.CountFastaNumSeqs(param_set.f_input);
  } else if(loader.IsFASTQ(param_set.f_input)) {
    num_reads = loader.CountFastqNumSeqs(param_set.f_input);
  } else  {
    cerr << "CLAN: unrecognized reads format, only FASTA or FASTQ format is suported. Abort." << endl;
    exit(1);
  }
  
  char **header = new char* [num_reads];
  char **seq = new char* [num_reads];
  if(loader.IsFASTA(param_set.f_input)) {
    num_reads = loader.LoadFasta(dna_alphabet, param_set.f_input, header, seq);
  } else  {
    num_reads = loader.LoadFastq(dna_alphabet, param_set.f_input, header, seq);
  }
  
  check_time = misc::MyTime();
  //PrintElapsed(start_time, check_time, "Loading index");
  start_time = misc::MyTime();

  BWTIDX num_hits = 0;
  // perform the search
  
  vector<vector<FragType> > selected_frag;
  vector<int> seq_order;
  
  MappingDriver map_driver;
  OutputPrinter out_printer;

  int num_chunks = 0;
  for(int i = 0; i < num_reads; i += MAP_CHUNK_SIZE) {
    ++ num_chunks;
  }
  out_fh.write((char*) &num_chunks, sizeof(int));  
  
  
  #pragma omp parallel num_threads(param_set.num_threads)
  {
    #pragma omp for
    for(int i = 0; i < num_reads; i += MAP_CHUNK_SIZE) {
      
      vector<vector<FragType> > *mapping_result = new vector<vector<FragType> >; 
      mapping_result->resize(MAP_CHUNK_SIZE);
      // clearing up the vector
      for(int j = 0; j < MAP_CHUNK_SIZE; ++ j) {
        (*mapping_result)[j].clear();
      }
      //cerr << "Entering mapping" << endl;
      map_driver.ReadMapBatch(
        ref_header, ref_seq, num_refs,
        ref_header_aux, ref_seq_aux, num_refs_aux, 
        header, seq, num_reads,
        i, MAP_CHUNK_SIZE, num_reads,
        bwt, bwt_aux, *mapping_result, param_set
      );
      # pragma omp critical
      {
        //cerr << "Preparing for output" << endl;
        out_printer.OutputEncoded(
          ref_seq, num_refs,
          i, MAP_CHUNK_SIZE, num_reads,
          *mapping_result, out_fh
        );
        // collect memory from mapping_result
        for(int k = 0; k < (*mapping_result).size(); ++ k)   {
          for(int l = 0; l < (*mapping_result)[k].size(); ++ l)   {
            for(int m = 0; m < (*mapping_result)[k][l].targets.size(); ++ m)   {
              delete [] (*mapping_result)[k][l].targets[m].cigar;
            }
          }
        }
        //cerr << "Output done" << endl;
        delete mapping_result;
      }
    }
  }
  out_fh.close();
  check_time = misc::MyTime();
  //PrintElapsed(start_time, check_time, "Writing all results");
  
  //cerr << "Preparing for input" << endl; 
  //vector<vector<FragType> > check_load;
  //out_printer.ReadEncoded(
  //  ref_seq, num_refs,
  //  num_reads, string(param_set.f_output), check_load
  //);
  //cerr << "Input done" << endl;
  
  for(int i = 0; i < num_refs; ++ i) {
    delete [] ref_header[i]; delete [] ref_seq[i];
  }
  delete [] ref_header; delete [] ref_seq;
  for(int i = 0; i < num_refs_aux; ++ i) {
    delete [] ref_header_aux[i]; delete [] ref_seq_aux[i];
  }
  delete [] ref_header_aux; delete [] ref_seq_aux;
  for(int i = 0; i < num_reads; ++ i) {
    delete [] header[i]; delete [] seq[i];
  }
  delete [] header; delete [] seq;
  
  return 0; 
}
