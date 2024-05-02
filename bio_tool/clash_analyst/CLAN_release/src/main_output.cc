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

using namespace std;

static CLANOutputParam param_set;

int main(int argc, char **argv)  {
  
  param_set.mode = 1;
  param_set.f_reference[0] = '\0';
  param_set.F_reference[0] = '\0';
  param_set.f_reads[0] = '\0';
  param_set.f_input[0] = '\0';
  param_set.f_output[0] = '\0';


  int copt;	
	extern char *optarg;
  extern int optind;
  while ((copt=getopt(argc,argv,"i:m:f:F:r:o:h")) != EOF)	{
    switch(copt) {
      case 'i':
        sscanf(optarg, "%s", param_set.f_input);
        continue;
      case 'o':
        sscanf(optarg, "%s", param_set.f_output);
        continue;
      case 'f':
        sscanf(optarg, "%s", param_set.f_reference);
        continue;
      case 'F':
        sscanf(optarg, "%s", param_set.F_reference);
        continue;
      case 'r':
        sscanf(optarg, "%s", param_set.f_reads);
        continue;
      case 'm':
        sscanf(optarg, "%d", &param_set.mode);
        continue;
			case 'h':
			default:
        cout << "==========================================================" << endl;
        cout << "\tCLAN: the CrossLinked reads ANalaysis tool" << endl;
        cout << "==========================================================" << endl;
        cout << endl;
        cout << "usage: clan_output -i [CLAN_FILE] -o [OUTPUT_FILE] -f [REFERENCE_FILE] -r [READ_FILE]" << endl;
        cout << endl;
        cout << "\ti: the .clan file output by clan_search (mandatory)" << endl;
        cout << "\to: the file for writing the mapping results (mandatory)" << endl;
        cout << "\tf: the main reference sequence database in FASTA format (e.g. the human genome, mandatory)" << endl;
        cout << "\tF: the auxiliary reference sequence database in FASTA format (optional)" << endl;
        cout << "\tr: the reads in FASTA format (madatory)" << endl; 
        //cout << "\tm: output mode (optional, 1: CLAN format; 2: BLAST tab format, 3: SAM format. default: 1)" << endl;
        cout << "\th: print this help message" << endl << endl;
        exit(0);
		}
		optind--;
	}	
  
  if(strlen(param_set.f_reference) <= 0 || strlen(param_set.f_input) <= 0 
      || strlen(param_set.f_output) <= 0 || strlen(param_set.f_reads) <= 0 
      || (param_set.mode != 1 && param_set.mode != 2 && param_set.mode != 3)
  )  {
    cerr << "Mandatory argument missing; please type \"clan_output -h\" to view the help information. Abort." << endl;
    exit(1);
  }
  ofstream out_fh;
  out_fh.open(param_set.f_output, ios::out);
  if(!out_fh.is_open())  {
    cerr << "CLAN:clan_output: unable to write to " << param_set.f_output << "; Abort." << endl;
    exit(1);
  }
  
  BioAlphabet dna_alphabet(DNA);
  Loader loader;
  
  int num_refs = loader.CountFastaNumSeqs(param_set.f_reference);
  char **ref_header = new char* [num_refs];
  char **ref_seq = new char* [num_refs];
  num_refs = loader.LoadFasta(dna_alphabet, param_set.f_reference, ref_header, ref_seq);
  
  int num_refs_aux = 0;
  char **ref_header_aux = NULL;
  char **ref_seq_aux = NULL;
  if(strlen(param_set.F_reference) > 0)  {
    num_refs_aux = loader.CountFastaNumSeqs(param_set.F_reference);
    ref_header_aux = new char* [num_refs_aux];
    ref_seq_aux = new char* [num_refs_aux];
    num_refs_aux = loader.LoadFasta(dna_alphabet, param_set.F_reference, ref_header_aux, ref_seq_aux);
  }
  
  int num_reads;
  if(loader.IsFASTA(param_set.f_reads)) {
    num_reads = loader.CountFastaNumSeqs(param_set.f_reads);
  } else if(loader.IsFASTQ(param_set.f_reads)) {
    num_reads = loader.CountFastqNumSeqs(param_set.f_reads);
  } else  {
    cerr << "CLAN: unrecognized reads format, only FASTA or FASTQ format is suported. Abort." << endl;
    exit(1);
  }
  char **header = new char* [num_reads];
  char **seq = new char* [num_reads];
  if(loader.IsFASTA(param_set.f_reads)) {
    num_reads = loader.LoadFasta(dna_alphabet, param_set.f_reads, header, seq);
  } else  {
    num_reads = loader.LoadFastq(dna_alphabet, param_set.f_reads, header, seq);
  }
  
  std::vector< std::vector<FragType> > mapping_results;
  OutputPrinter output_printer;
  output_printer.ReadEncoded(num_reads, string(param_set.f_input), mapping_results);
  //cerr << "Reads loaded!!!" << endl;
  //cerr << "num refs:  " << num_refs << endl;
  //cerr << ref_header[1619] << endl;  
  if(param_set.mode == 1)  {
    output_printer.OutputPrinterCLAN(
      ref_header, ref_seq, ref_header_aux, ref_seq_aux,
      header, seq,
      0, num_reads, num_reads,
      mapping_results,
      out_fh
    );
  } else if(param_set.mode == 2) {
    string concat_seq; 
    Concatenator concat_obj(ref_seq, num_refs, concat_seq);
    output_printer.OutputPrinterBLAST(
      ref_header, ref_seq,
      header, seq,
      0, num_reads, num_reads,
      mapping_results,
      (string) param_set.f_reference, concat_seq.length(),
      out_fh
    );
  } else if(param_set.mode == 3)  {
    string concat_seq; 
    Concatenator concat_obj(ref_seq, num_refs, concat_seq);
    output_printer.OutputPrinterSAM(
      ref_header, ref_seq,
      header, seq,
      0, num_reads, num_reads,
      mapping_results,
      concat_seq.length(),
      out_fh
    );
  }
  
  return 0;
} 

