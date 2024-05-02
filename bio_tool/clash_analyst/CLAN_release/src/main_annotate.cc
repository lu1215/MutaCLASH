#include <iostream>
#include <string>
#include <list>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bwt.h"
#include "bwt_search.h"
#include "loader.h"
#include "bio_alphabet.h"
#include "minimizer_sort.h"
#include "parameters.h"
#include "annotate_driver.h"

using namespace std;

static CLANAnnotateParam param_set;

int main(int argc, char **argv) {

  param_set.strand_specific = false;
  param_set.print_simplex = false;
  param_set.extend_len = 10;
  param_set.min_island_len = 20;
  param_set.max_island_len = 200;
  param_set.max_island_gap = 10;
  param_set.min_coverage = 0;
  param_set.max_dimer_dG = -10.0;
  param_set.min_seed_len = 5;
  param_set.num_threads = 1;
  param_set.f_reference[0] = '\0';
  param_set.f_reads[0] = '\0';
  param_set.f_input[0] = '\0';
  param_set.f_output[0] = '\0';
  
  
  int copt;	
	extern char *optarg;
  extern int optind;
  while ((copt=getopt(argc,argv,"i:f:r:o:t:l:x:g:c:e:k:sph")) != EOF)	{
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
      case 'r':  
        sscanf(optarg, "%s", param_set.f_reads);
        continue;
      case 't':  
        sscanf(optarg, "%d", &param_set.num_threads);
        continue;
      case 'l':  
        sscanf(optarg, "%d", &param_set.min_island_len);
        continue;
      case 'x':  
        sscanf(optarg, "%d", &param_set.max_island_len);
        continue;
      case 'g':  
        sscanf(optarg, "%d", &param_set.max_island_gap);
        continue;
      case 'c':  
        sscanf(optarg, "%d", &param_set.min_coverage);
        continue;
      case 'e':  
        sscanf(optarg, "%f", &param_set.max_dimer_dG);
        continue;
      case 'k':  
        sscanf(optarg, "%d", &param_set.min_seed_len);
        continue;
      case 's':
        param_set.strand_specific = true;
        continue;
      case 'p':
        param_set.print_simplex = true;
        continue;
			case 'h':
			default:
			  cout << endl;
        cout << "==========================================================" << endl;
        cout << "\tCLAN: the CrossLinked reads ANalaysis tool" << endl;
        cout << "==========================================================" << endl;
        cout << endl;
        cout << "usage: clan_annotate -i [CLAN_FILE] -o [OUTPUT_FILE] -f [REFERENCE_FILE] -r [READ_FILE]" << endl;
        cout << endl;
        cout << "\ti: the .clan file output by clan_search (mandatory)" << endl;
        cout << "\to: the file for writing the annotation results (mandatory)" << endl;
        cout << "\tf: the reference sequence in FASTA format (e.g. the human genome, mandatory)" << endl;
        cout << "\tr: the reads in FASTA format (madatory)" << endl; 
        cout << "\ts: enable strand-specific mapping (ignore mapping to reverse strand, default FALSE)" << endl;
        cout << "\tp: only print the islands (ignore any duplex information, default FALSE)" << endl;
        cout << "\tl: mininum length of an island (optional, default 20)" << endl;
        cout << "\tx: maximum length of an island (optional, default 200)" << endl;
        cout << "\tg: maximum length of uncovered bases within an island (optional, default 10)" << endl;
        cout << "\tc: minimum number of supporting reads for duplex (optional, default 0, i.e. output all)" << endl;
        cout << "\te: maximun dimer free energy of the duplex (optional, default -10.0)" << endl;
        cout << "\tk: minimum length of the perfect stack in duplex (optional, default 5)" << endl;
        cout << "\tt: number of threads to use (optional, default 1)" << endl;
        cout << "\th: print this help message" << endl << endl;
        exit(0);
		}
		optind--;
	}	
	
	
	if(strlen(param_set.f_reference) <= 0 || strlen(param_set.f_input) <= 0 || strlen(param_set.f_output) <= 0)  {
    cerr << "Mandatory argument missing; please type \"clan_annotate -h\" to view the help information. Abort." << endl;
    exit(1);
  }
  ofstream out_fh;
  out_fh.open(param_set.f_output, ios::out);
  if(!out_fh.is_open())  {
    cerr << "CLAN:clan_annotate: unable to write to " << param_set.f_output << "; Abort." << endl;
    exit(1);
  }
  out_fh.close();
  
  
  // load the reference sequence and reads
  BioAlphabet dna_alphabet(DNA);
  Loader loader;
  
  
  int num_refs = loader.CountFastaNumSeqs(param_set.f_reference);
  char **ref_header = new char* [num_refs];
  char **ref_seq = new char* [num_refs];
  num_refs = loader.LoadFasta(dna_alphabet, param_set.f_reference, ref_header, ref_seq);
  
  cerr << "Finish loading reference" << endl;
  
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
  
  cerr << "Finish loading reads" << endl;
  //cerr << seq[23] << endl;
  //exit(1);
  
  
  // read in the CLAN mapping output
  std::vector< std::vector<FragType> > *mapping_results = new std::vector< std::vector<FragType> >;
  OutputPrinter output_printer;
  output_printer.ReadEncoded(num_reads, string(param_set.f_input), *mapping_results);
  cerr << "Finish loading input" << endl;

  // checked here
  
  // estimate the total read length by checking the first 1000 reads
  int n = num_reads > 1000 ? 1000 : num_reads;
  int total_len = 0;
  for(int i = 0; i < n; ++ i) {
    total_len += strlen(seq[i]);
  } 
  int ave_len = total_len / n;
  
  // perform the annotation
  AnnotateDriver *annotate_driver = new AnnotateDriver;
  annotate_driver->AnnotateMapping(
    ref_header, ref_seq, num_refs, 
    header, seq, num_reads,
    ave_len, mapping_results, param_set
  );
  annotate_driver->CleanMappingCIGAR(*mapping_results);
  
  delete mapping_results;
  delete annotate_driver;
  
  // collect memory
  for(int i = 0; i < num_refs; ++ i) {
    delete [] ref_header[i]; delete [] ref_seq[i];
  }
  delete [] ref_header; delete [] ref_seq;
  for(int i = 0; i < num_reads; ++ i) {
    delete [] header[i]; delete [] seq[i];
  }
  delete [] header; delete [] seq;
  
  
  return 0;
}
