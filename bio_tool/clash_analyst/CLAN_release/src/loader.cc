#include "loader.h"

using namespace std;

Loader::Loader()  {
  return;
}

Loader::~Loader() {
  return;
}

bool Loader::IsFASTA(const char *file_name) {
  std::ifstream ifstrm(file_name, std::ios_base::in);
	if (!ifstrm) {
    std::cerr << "Loader::Cannot open file: " << file_name << std::endl;
    exit (1);
  }
  std::string line;
  std::getline(ifstrm, line);
  ifstrm.close();
  if(line[0] == '>') return true;
  return false;
}

bool Loader::IsFASTQ(const char *file_name) {
  std::ifstream ifstrm(file_name, std::ios_base::in);
	if (!ifstrm) {
    std::cerr << "Loader::Cannot open file: " << file_name << std::endl;
    exit (1);
  }
  std::string line;
  std::getline(ifstrm, line);
  ifstrm.close();
  if(line[0] == '@') return true;
  return false;
}

int Loader::CountFastaNumSeqs(const char *file_name)  {
  std::ifstream ifstrm(file_name, std::ios_base::in);
	if (!ifstrm) {
    std::cerr << "Loader::Cannot open file: " << file_name << std::endl;
    exit (1);
  }
  int count = 0;
  std::string line;
  while (std::getline(ifstrm, line)) {
    if(line[0] == '>') count++;
  }
	
  ifstrm.close();
  return count;
}

int Loader::CountFastqNumSeqs(const char *file_name)  {
  std::ifstream ifstrm(file_name, std::ios_base::in);
	if (!ifstrm) {
    std::cerr << "Loader::Cannot open file: " << file_name << std::endl;
    exit (1);
  }
  int count = 0;
  std::string line;
  while (std::getline(ifstrm, line)) {
    if(line[0] == '@') count++;
  }
  //cerr << "CountFASTQSeq: " << count << endl;
  ifstrm.close();
  return count;
}

void Loader::RecordSequence(
  char **header, char **seq, 
  std::string &single_header, std::string &single_seq, 
  int index
) {
  assert(single_header.length() > 0);
  assert(single_seq.length() > 0);
  
  header[index] = new char[single_header.length() + 1];
  strcpy(header[index], single_header.c_str());
  seq[index] = new char[single_seq.length() + 1];
  strcpy(seq[index], single_seq.c_str());
  return;
}

bool Loader::CheckSpecialChar(BioAlphabet &alphabet, std::string &sseq, float freq_cutoff) {
  if(sseq.length() <= 0)  return false;
  int num_special = 0;
  // special treatment to convert RNA sequence to DNA sequence
  for(int i = 0; i < sseq.length(); ++ i)  {
    if(sseq[i] == 'u' || sseq[i] == 'U')  {
        sseq[i] = 'T';
    }
    if(!alphabet.IsValid(sseq[i]))  {
      if(alphabet.IsValid(toupper(sseq[i])))  {
        sseq[i] = toupper(sseq[i]);
      } else  {
        ++ num_special;
        sseq[i] = alphabet.RandomChar();
      }
    }
  }
  if(num_special / sseq.length() > 1 - freq_cutoff) return false;
  return true;  
}

int Loader::LoadFasta(BioAlphabet &alphabet, const char *file_name, char **header, char **seq, const int min_len) {
  // opens the file and read line-by-line
  std::ifstream ifstrm(file_name, std::ios_base::in);
  std::string line, fasta_tag, fasta_seq;
  int count = 0;
  while (std::getline(ifstrm, line)) {
    if (line[0] == '>') {
      if (fasta_tag != "" && fasta_seq != "") {
        CheckSpecialChar(alphabet, fasta_seq);
        if(fasta_seq.length() > min_len)  {
          RecordSequence(header, seq, fasta_tag, fasta_seq, count);
          ++ count;
        }
      }
      fasta_tag = line.substr(1, line.length() - 1); fasta_seq = ""; 
    } else fasta_seq += line;
  }
  ifstrm.close();	
  // handle the last sequence		
  if (fasta_tag != "" && fasta_seq != "") {
    CheckSpecialChar(alphabet, fasta_seq);
    if(fasta_seq.length() > min_len)  {
      RecordSequence(header, seq, fasta_tag, fasta_seq, count);
      ++ count;
    }
  }		
  return count;
}

int Loader::LoadFastq(BioAlphabet &alphabet, const char *file_name, char **header, char **seq, const int min_len) {
  // opens the file and read line-by-line
  std::ifstream ifstrm(file_name, std::ios_base::in);
  std::string line, fasta_tag, fasta_seq;
  int count = 0;
  while (std::getline(ifstrm, line)) {
    if (line[0] == '@') {
      line = line.substr(1);
      std::string read_seq;
      std::getline(ifstrm, read_seq);
      //cerr << "read_seq:  " << read_seq << endl;
      CheckSpecialChar(alphabet, read_seq);
      if(read_seq.length() > min_len)  {
        RecordSequence(header, seq, line, read_seq, count);
        ++ count;
      }
      //cerr << "after recording" << endl;
      // skipping the quality information
      std::getline(ifstrm, line);
      std::getline(ifstrm, line);
      //cerr << "count: " << count << endl;
    } 
  }
  ifstrm.close();		
  return count;
}


// Obsolete code block
/*
void Loader::SortFasta(int n, char **header, char **seq)  {
  // make a copy of the sequences
  char **header_holder = new char* [n];
  char **seq_holder = new char* [n];
  unordered_map<string, list<int> > orders;
  int mer_len = 6;
  int i, j;
  for(i = 0; i < n; ++ i) {
    header_holder[i] = header[i];
    seq_holder[i] = seq[i];
    string s = seq[i];
    string min_str = string(mer_len, 'Z');
    for(j = 0; j < s.length() - mer_len + 1; ++ j) {
      if(s.substr(j, mer_len) < min_str)
        min_str = s.substr(j, mer_len);
    }
    orders[min_str].push_back(i);
  }
  // copy the sequences
  int nn = 0;
  for(auto it = orders.begin(); it != orders.end(); ++ it) {
    for(auto it_l = it->second.begin(); it_l != it->second.end(); ++ it_l) {
      header[nn] = header_holder[*it_l];
      seq[nn] = seq_holder[*it_l];
      ++ nn;
    }
  }
  // release the memory
  delete [] header_holder;
  delete [] seq_holder;
  return; 
}
*/
