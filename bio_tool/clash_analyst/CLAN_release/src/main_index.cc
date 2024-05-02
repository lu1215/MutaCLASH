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

using namespace std;

/*
void PrintTime()  {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  std::cout << "currentDateTime()=" << buf << std::endl;
  return;
}
*/

bool check_index_valid(char *f_index)  {
  struct stat s_f_index;
  if( stat(f_index, &s_f_index) == 0 )  {
    if( s_f_index.st_mode & S_IFDIR ) {
      cerr << "clan_index: Error: The output index file prefix should not be a directory, please reset output index prefix. Abort" << endl;
      cerr << "\tExample: you should input \"/home/YourDirectory/fileID\" instead of \"/home/YourDirectory\"." << endl;
      return false;
    }
  }
  return true;
}

static CLANIndexParam param_set;

int main(int argc, char **argv)  {
  
  param_set.f_reference[0] = '\0';
  param_set.f_index[0] = '\0';

  int copt;	
	extern char *optarg;
  extern int optind;
  while ((copt=getopt(argc,argv,"f:d:h")) != EOF)	{
    switch(copt) {
      case 'f':
        sscanf(optarg, "%s", param_set.f_reference);
        continue;
      case 'd':
        sscanf(optarg, "%s", param_set.f_index);
        continue;
			case 'h':
			default:
        cout << "==========================================================" << endl;
        cout << "\tCLAN: the CrossLinked reads ANalaysis tool" << endl;
        cout << "==========================================================" << endl;
        cout << endl;
        cout << "usage: clan_index -f [REFERENCE_FILE] -d [INDEX_PREFIX]" << endl;
        cout << endl;
        cout << "\tf: the reference genome(s)/sequence database (mandatory)" << endl;
        cout << "\td: the prefix of the indexes (mandatory)" << endl;
        cout << "\th: print this help message" << endl << endl;
        exit(0);
		}
		optind--;
	}	

  // check argument setting
  if(strlen(param_set.f_reference) <= 0 || strlen(param_set.f_index) <= 0)  {
    cout << "Mandatory argument missing; please type \"clan_index -h\" to view the help information." << endl;
    cout << "Abort." << endl;
    exit(0);
  }

  // check if the specified index prefix is valid
  if(!check_index_valid(param_set.f_index))  {
    exit(1);
  }
  ofstream out_fh;
  string prefix = param_set.f_index;
  string idx = prefix + string(".bwt");
  out_fh.open(idx.c_str(), ios::out);
  if(!out_fh.is_open())  {
    cout << "clan_index: Error: unable to write index files" << "; Abort." << endl;
    exit(1);
  }
  out_fh.close();

  BioAlphabet dna_alphabet(DNA);
  Loader fasta_loader;
  
  int num_refs = fasta_loader.CountFastaNumSeqs(param_set.f_reference);
  char **ref_header = new char* [num_refs];
  char **ref_seq = new char* [num_refs];
  num_refs = fasta_loader.LoadFasta(dna_alphabet, param_set.f_reference, ref_header, ref_seq);
  //cout << "Finish loading data" << endl;
  
  string concat_seq; 
  Concatenator concat_obj(ref_seq, num_refs, concat_seq);
  //cout << concat_seq << endl;
  BWT bwt;
  bwt.Construct(dna_alphabet, concat_seq.c_str()); 
  //cerr << "Finish constructing index" << endl;
  bwt.WriteIndex(param_set.f_index);
  
  //cerr << "Finish constructing forward BWT: " << bwt.GetSize() << endl; 
  
  //concat_seq = string(concat_seq.rbegin(), concat_seq.rend());
  //BWT rev_bwt;
  //rev_bwt.ConstructNoPos(dna_alphabet, concat_seq.c_str());
  //rev_bwt.WriteBWTIndex("index.rev.bwt"); 
  //cout << "Finish constructing reverse BWT: " << bwt.GetSize() << endl; 

  for(int i = 0; i < num_refs; ++ i) {
    delete [] ref_header[i]; delete [] ref_seq[i];
  }
  delete [] ref_header; delete [] ref_seq;
  
  bwt.Purge();
  return 0;
}
