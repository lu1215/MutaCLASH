#include "bwt.h"

using namespace std;


void BWT::BuildIndex(const BWTCHAR* text)  {
  str_ = new BWTCHAR [size_ + 1];
  strcpy((char*) str_, (char*) text);
  bwt_ = new BWTCHAR [size_ + 1];
  BWTIDX fm_size = (size_ / BWT_FM_GAP) + 1;
  BWTIDX i, j;
  // allocate space for FM-index
  for(int i = 0; i < alphabet_.GetSize(); ++ i)  {
    fm_index_[alphabet_.GetInvCharMap(i)] = new BWTIDX [fm_size];
    for(j = 0; j < fm_size; ++ j) fm_index_[alphabet_.GetInvCharMap(i)][j] = 0;
  }
  fm_index_[DELIM] = new BWTIDX [fm_size]; 
  for(j = 0; j < fm_size; ++ j) fm_index_[DELIM][j] = 0;
  // allocate space for position
  // count frequencies
  for(i = 0; i < size_; ++ i) {
    char c;
    if(position_[i] == 0) c = DELIM;
    else c = text[position_[i] - 1];
    // pad on the "bwt_" array
    bwt_[i] = c;  
    if(i % BWT_FM_GAP == 0)  {
      // copy the frequency to the FM index
      BWTIDX ix = i / BWT_FM_GAP;
      for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
        if(fm_index_.find(it->first) == fm_index_.end())  {
          cerr << "Error: BWT::BuildIndex: Unrecognized character in text/sequence: " << it->first << "; abort." << endl; exit(1);
        }
        // add the frequency to the FM-index
        fm_index_[it->first][ix] = it->second;
      }
      // record gapped position information
    }
    // add frequency
    ++ char_freq_[c];
  }
  bwt_[i] = '\0';
  // handle spacial case
  if(i % BWT_FM_GAP == 0) {
    // copy the frequency to the FM index
    BWTIDX ix = i / BWT_FM_GAP;
    for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
      if(fm_index_.find(it->first) == fm_index_.end())  {
        cerr << "Error: BWT::BuildIndex: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
        exit(1);
      }
      // add the frequency to the FM-index
      fm_index_[it->first][ix] = it->second;
    }
    // record gapped position information
  }  
  // compute accumulated index
  acc_freq_ = new BWTIDX [alphabet_.GetSize() + 2]; // note the trailing character '$'
  // the begin of '$' in the array is 0
  acc_freq_[0] = 0; acc_freq_[1] = char_freq_[DELIM];
  for(i = 0; i < alphabet_.GetSize(); ++ i) {
    acc_freq_[i + 2] = acc_freq_[i + 1] + char_freq_[alphabet_.GetInvCharMap(i)];
  }
  // check if character that do not found in the alphabet exists
  for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
    if(it->first != DELIM && !alphabet_.IsValid((char) it->first))  {
      cerr << "Error: BWT::BuildIndex: Unrecognized character found in the text/sequence: " << it->first << "; abort." << endl;
      exit(1);
    }
  }
  return;
}

void BWT::BuildIndexNoPos(const BWTCHAR* text)  {
  BWTIDX fm_size = (size_ / BWT_FM_GAP) + 1;
  BWTIDX i, j;
  // allocate space for FM-index
  for(i = 0; i < alphabet_.GetSize(); ++ i)  {
    fm_index_[alphabet_.GetInvCharMap(i)] = new BWTIDX [fm_size];
    for(j = 0; j < fm_size; ++ j) fm_index_[alphabet_.GetInvCharMap(i)][j] = 0;
  }
  fm_index_[DELIM] = new BWTIDX [fm_size];
  for(j = 0; j < fm_size; ++ j) fm_index_[DELIM][j] = 0;
  // count frequencies
  for(i = 0; i < size_; ++ i) { 
    if(i % BWT_FM_GAP == 0)  {
      // copy the frequency to the FM index
      BWTIDX ix = i / BWT_FM_GAP;
      for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
        if(fm_index_.find(it->first) == fm_index_.end())  {
          cerr << "Error: BWT::BuildIndexNoPos: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
          exit(1);
        }
        // add the frequency to the FM-index
        fm_index_[it->first][ix] = it->second;
      }
    }
    // add frequency
    ++ char_freq_[bwt_[i]];
  }
  // handle spacial case
  if(i % BWT_FM_GAP == 0) {
    // copy the frequency to the FM index
    BWTIDX ix = i / BWT_FM_GAP;
    for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
      if(fm_index_.find(it->first) == fm_index_.end())  {
        cerr << "Error: BWT::BuildIndexNoPos: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
        exit(1);
      }
      // add the frequency to the FM-index
      fm_index_[it->first][ix] = it->second;
    }
  }  
  // compute accumulated index
  acc_freq_ = new BWTIDX [alphabet_.GetSize() + 2]; // note the trailing character '$'
  // the begin of '$' in the array is 0
  acc_freq_[0] = 0; acc_freq_[1] = char_freq_[DELIM];
  for(i = 0; i < alphabet_.GetSize(); ++ i) {
    acc_freq_[i + 2] = acc_freq_[i + 1] + char_freq_[alphabet_.GetInvCharMap(i)];
  }
  // check if character that do not found in the alphabet exists
  for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
    if(it->first != DELIM && !alphabet_.IsValid((char) it->first))  {
      cerr << "Error: BWT::BuildIndexNoPos: Unrecognized character found in the text/sequence: " << it->first << "; abort." << endl;
      exit(1);
    }
  }
  return;
}

void BWT::ConstructLocationInfo(const int n, char **header, char **seq) {
  if(!build_success_) {
    cout << "BWT::BuildLocationInfo: BWT needs to be first built, abort." << endl;
    exit(1);
  }
  num_seqs_ = n;
  header_ = new BWTCHAR* [n];
  acc_len_ = new BWTIDX [n];
  for(int i = 0; i < n; ++ i) {
    header_[i] = new BWTCHAR [strlen(header[i]) + 1];
    strcpy((char*) header_[i], header[i]);
    // note that we need to add 1 to the sequence length because of the delimitor
    if(i > 0)  {  acc_len_[i] = acc_len_[i - 1] + strlen(seq[i]) + 1; } 
    else  { acc_len_[i] = strlen(seq[i]) + 2; }
  }
  if(acc_len_[n - 1] != size_)  {
    cout << "BWT::BuildLocationInfo: the sizes of the sequence set and the BWT are different, abort." << endl;
    cout << acc_len_[n - 1] << "  " << size_ << endl; 
    exit(1);
  }
  ref_success_ = true;
  return;
}

/*
BWTIDX BWT::GetRefLocation(
    const BWTIDX bwt_pos, std::string &ref_header, BWTIDX &ref_pos
)  {
  // check if the reference information has been built
  assert(ref_success_ && bwt_pos < acc_len_[num_seqs_ - 1]);
  BWTIDX pos = bwt_pos;
  BWTIDX traversed = 0;
  while(pos % BWT_FM_GAP != 0) {
    //cout << "pos: " << pos << endl;
    //cout << "char: " << bwt_[pos] << endl;
    BWTIDX r = CountOccurrence(bwt_[pos], pos);
    //cout << "occurrence: " << r << endl;
    int c_id = alphabet_.GetCharMap(bwt_[pos]);
    //cout << "cid: " << c_id << endl;
    pos = acc_freq_[c_id + 1] + r;
    //cout << "char begin:  " << acc_freq_[c_id + 1] << endl;
    //cout << "new pos: " << pos << endl;
    ++ traversed;
  }
  BWTIDX str_pos = (position_[pos / BWT_FM_GAP] + traversed) % size_;
  //cout << "pos info:  " << str_pos << " " << position_[pos % BWT_FM_GAP] << " " << traversed << endl;
  // perform binary search
  int left = 0, right = num_seqs_ - 1, mid = 0;
  while(left <= right) {
    mid = (left + right) / 2;
    //cout << "left:  " << left << "  right:  " << right << endl; 
    //cout << "mid: " << mid << " " << str_pos << " " << acc_len_[mid] << endl;
    if(str_pos >= acc_len_[mid]) {
      left = mid + 1;
    } else {
      if(right > 0 && mid >= 1 && str_pos < acc_len_[mid - 1]) right = mid - 1;
      else  break;
    }
  }
  //cout << "Identified ID: " << mid << " " << acc_len_[mid - 1] << endl;
  ref_header = string((char*) header_[mid]);
  if(mid == 0)  ref_pos = str_pos - 1; // because DILIM was added as the first char of the generalized string
  else ref_pos = str_pos - acc_len_[mid - 1];  
  return str_pos;
}
*/

/*
BWTIDX BWT::BWTToStrPos(const BWTIDX bwt_pos) {
  // check if the reference information has been built
  assert(ref_success_ && bwt_pos < acc_len_[num_seqs_ - 1]);
  BWTIDX pos = bwt_pos;
  BWTIDX traversed = 0;
  while(pos % BWT_FM_GAP != 0) {
    //cout << "pos: " << pos << endl;
    //cout << "char: " << bwt_[pos] << endl;
    BWTIDX r = CountOccurrence(bwt_[pos], pos);
    //cout << "occurrence: " << r << endl;
    int c_id = alphabet_.GetCharMap(bwt_[pos]);
    //cout << "cid: " << c_id << endl;
    pos = acc_freq_[c_id + 1] + r;
    //cout << "char begin:  " << acc_freq_[c_id + 1] << endl;
    //cout << "new pos: " << pos << endl;
    ++ traversed;
  }
  BWTIDX str_pos = (position_[pos / BWT_FM_GAP] + traversed) % size_;
  return str_pos;
}
*/

BWTIDX BWT::BWTToStrPos(const BWTIDX bwt_pos) {
  // check if the reference information has been built
  assert(ref_success_ && bwt_pos < acc_len_[num_seqs_ - 1]);
  return position_[bwt_pos];
}


void BWT::StrPosToRefPos(const BWTIDX str_pos, std::string &ref_header, BWTIDX &ref_pos)  {
  int left = 0, right = num_seqs_ - 1, mid = 0;
  while(left <= right) {
    mid = (left + right) / 2;
    //cout << "left:  " << left << "  right:  " << right << endl; 
    //cout << "mid: " << mid << " " << str_pos << " " << acc_len_[mid] << endl;
    if(str_pos >= acc_len_[mid]) {
      left = mid + 1;
    } else {
      if(right > 0 && mid >= 1 && str_pos < acc_len_[mid - 1]) right = mid - 1;
      else  break;
    }
  }
  //cout << "Identified ID: " << mid << " " << acc_len_[mid - 1] << endl;
  ref_header = string((char*) header_[mid]);
  if(mid == 0)  ref_pos = str_pos - 1; // because DILIM was added as the first char of the generalized string
  else ref_pos = str_pos - acc_len_[mid - 1]; 
  
  if(ref_pos < 0)  {
    cerr << str_pos << "  " << ref_header << "  " << ref_pos << endl;
    exit(1);
  }
   
  return;
}

void BWT::StrPosToRefPos(const BWTIDX str_pos, int &ref_header_idx, BWTIDX &ref_pos)  {
  int left = 0, right = num_seqs_ - 1, mid = 0;
  while(left <= right) {
    mid = (left + right) / 2;
    //cout << "left:  " << left << "  right:  " << right << endl; 
    //cout << "mid: " << mid << " " << str_pos << " " << acc_len_[mid] << endl;
    if(str_pos >= acc_len_[mid]) {
      left = mid + 1;
    } else {
      if(right > 0 && mid >= 1 && str_pos < acc_len_[mid - 1]) right = mid - 1;
      else  break;
    }
  }
  //cout << "Identified ID: " << mid << " " << acc_len_[mid - 1] << endl;
  ref_header_idx = mid;
  if(mid == 0)  ref_pos = str_pos - 1; // because DILIM was added as the first char of the generalized string
  else ref_pos = str_pos - acc_len_[mid - 1]; 
  
  if(ref_pos < 0)  {
    cerr << "CLAN:BWT::StrPosToRefPos out of range: " << str_pos << "  " << ref_header_idx << "  " << ref_pos << endl;
    exit(1);
  }
   
  return;
}

void BWT::Construct(const BioAlphabet &alphabet, const char *text) {
  alphabet_ = alphabet;
  size_ = string(text).length();
  //cout << text << " " << size_ << endl;
  // build suffix array 
  position_ = new BWTIDX [size_ + 1]; 
  int divsufsort_tag = divsufsort64((BWTCHAR*) text, position_, size_);
  //for(int i = 0; i < size_; ++ i) {
  //  cout << tmp_position_[i] << endl;
  //}
  assert(divsufsort_tag == 0);
  // build indexes for BWT
  //cout << "Finish sorting suffixes..." << endl;
  BuildIndex((BWTCHAR*) text);   
  build_success_ = true;
  return;
}

void BWT::ConstructNoPos(const BioAlphabet &alphabet, const char *text) {
  alphabet_ = alphabet;
  size_ = string(text).length();
  // build suffix array 
  bwt_ = new BWTCHAR [size_ + 1];
  int divsufsort_tag = divbwt64((BWTCHAR*) text, bwt_, NULL, size_);
  bwt_[size_] = '\0';
  // TODO: Handle sorting errors
  assert(divsufsort_tag >= 0);
  // build indexes for BWT
  BuildIndexNoPos((BWTCHAR*) text);
  build_success_ = true;
  return;
}

void BWT::WriteIndex(const std::string &file_prefix)  {
  WriteBWTIndex(file_prefix + ".bwt");
  WriteFMIndex(file_prefix + ".fmi");
  WritePosition(file_prefix + ".pos");
  return;  
}

void BWT::WriteIndexNoPos(const std::string &file_prefix) {
  WriteBWTIndex(file_prefix + ".bwt");
  WriteFMIndex(file_prefix + ".fmi");
  return;
}

void BWT::ConstructNoPosFromIndex(
    const BioAlphabet &alphabet, const char *text, 
    const std::string &file_prefix
)  {
  alphabet_ = alphabet;
  size_ = string(text).length();
  LoadBWTIndex(file_prefix + ".bwt");
  LoadFMIndex(file_prefix + ".fmi");
  // compute accumulated index
  for(int i = 0; i < size_; ++ i) {
    if(alphabet_.IsValid(bwt_[i]) || bwt_[i] == DELIM)  char_freq_[bwt_[i]] ++;
    else {
      cout << "BWT::ConstructNoPosFromIndex: unknow character encoutered: " << bwt_[i] << "; abort" << endl; exit(1);
    }
  }
  acc_freq_ = new BWTIDX [alphabet_.GetSize() + 2]; // note the trailing character '$'
  // the begin of '$' in the array is 0
  acc_freq_[0] = 0; acc_freq_[1] = char_freq_[DELIM];
  for(int i = 0; i < alphabet_.GetSize(); ++ i) {
    if(char_freq_.find(alphabet_.GetInvCharMap(i)) != char_freq_.end())  {
      acc_freq_[i + 2] = acc_freq_[i + 1] + char_freq_[alphabet_.GetInvCharMap(i)];
    } else  {
      acc_freq_[i + 2] = acc_freq_[i + 1];
    }
  }
  build_success_ = true;
  return;
}

void BWT::ConstructFromIndex(
    const BioAlphabet &alphabet, const char *text, 
    const std::string &file_prefix
)  {
  alphabet_ = alphabet;
  size_ = string(text).length();
  str_ = new BWTCHAR [size_ + 1];
  strcpy((char*) str_, text);
  LoadBWTIndex(file_prefix + ".bwt");
  //cout << "Finished loading BWT index" << endl;
  LoadFMIndex(file_prefix + ".fmi");
  //cout << "Finished loading FM index" << endl;
  LoadPosition(file_prefix + ".pos");
  //cout << "Finished loading position" << endl;
  // compute accumulated index
  for(int i = 0; i < size_; ++ i) {
    if(alphabet_.IsValid(bwt_[i]) || bwt_[i] == DELIM)  char_freq_[bwt_[i]] ++;
    else {
      cout << "BWT::ConstructNoPosFromIndex: unknow character encoutered: " << bwt_[i] << "; abort" << endl; exit(1);
    }
  }
  acc_freq_ = new BWTIDX [alphabet_.GetSize() + 2]; // note the trailing character '$'
  // the begin of '$' in the array is 0
  acc_freq_[0] = 0; acc_freq_[1] = char_freq_[DELIM];
  for(int i = 0; i < alphabet_.GetSize(); ++ i) {
    if(char_freq_.find(alphabet_.GetInvCharMap(i)) != char_freq_.end())  {
      acc_freq_[i + 2] = acc_freq_[i + 1] + char_freq_[alphabet_.GetInvCharMap(i)];
    } else  {
      acc_freq_[i + 2] = acc_freq_[i + 1];
    }
  }
  
  build_success_ = true;
  return;
}

void BWT::WriteBWTIndex(const std::string &file_name) {
  ofstream out_fh;
  out_fh.open(file_name, ios::out);
  if(!out_fh.is_open())  {
    cout << "BWT::WriteBWTIndex: Error in writing BWT indexing file " << file_name << "; Abort." << endl;
  }
  // writing reduced_map
  if(size_ != string((char*) bwt_).length()) {
    cout << "BWT::WriteBWTIndex: sizes of text and BWT do not match, abort." << endl;
  }
  out_fh << (char*) bwt_ << endl;
  out_fh.close();
  return;
}

void BWT::LoadBWTIndex(const std::string &file_name)  {
  // loading BWT index from file
  string line;
  ifstream in_fh;
  in_fh.open(file_name, ios::in);
  if(!in_fh.is_open())  {
    cout << "BWT::ConstructNoPosFromIndex: Error in loading BWT indexing file " << file_name << "; Abort." << endl;
    exit(1);
  }
  getline(in_fh, line);
  // if bwt_ has been previously initialized, reallocate space
  if(build_success_)  { 
    cout << "BWT::ConstructNoPosFromIndex: conflicting BWT construction; abort" << endl;
    exit(1);
  }  
  if(line.length() != size_)  {
    cout << "BWT::ConstructNoPosFromIndex: sizes of text and BWT index do not match; abort." << endl;
    exit(1);
  }
  bwt_ = new BWTCHAR [line.length() + 1];
  strcpy((char*) bwt_, line.c_str());
  in_fh.close();
  return;
}

void BWT::WritePosition(const std::string &file_name) {
  ofstream out_fh;
  out_fh.open(file_name, ios::out);
  if(!out_fh.is_open())  {
    cout << "BWT::WritePosition: Error in writing BWT position file " << file_name << "; Abort." << endl;
  }
  out_fh.write((char*) position_, size_ * sizeof(BWTIDX));
  out_fh.close();
  return;
}

void BWT::LoadPosition(const std::string &file_name)  {
  int i, j;
  ifstream in_fh;
  in_fh.open(file_name, ios::in);
  if(!in_fh.is_open())  {
    cout << "BWT::LoadPosition: failed to open BWT position indexing file: " << file_name << "; Abort." << endl; exit(1);
  }
  // loading the files line by line
  position_ = new BWTIDX [size_];
  in_fh.read((char*) position_, size_ * sizeof(BWTIDX));
  in_fh.close();
  return;
}

void BWT::WriteFMIndex(const std::string &file_name) {
  ofstream out_fh;
  out_fh.open(file_name, ios::out);
  if(!out_fh.is_open())  {
    cout << "BWT::WriteFMIndex: Error in writing FM index file " << file_name << "; Abort." << endl;
  }
  // writing reduced_map
  BWTIDX gap_size = (size_ / BWT_FM_GAP) + 1;
  out_fh << gap_size << endl;
  for(auto it = fm_index_.begin(); it != fm_index_.end(); ++ it) {
    out_fh << it->first << ":";
    for(BWTIDX i = 0; i < gap_size; ++ i) {
      out_fh << it->second[i] << ":";
    }
    out_fh << endl;
  }
  out_fh.close();
  return;
}

void BWT::LoadFMIndex(const std::string &file_name) {
  string line;
  ifstream in_fh;
  in_fh.open(file_name, ios::in);
  if(!in_fh.is_open())  {
    cout << "BWT::LoadFMIndex: failed to open FM-index file: " << file_name << "; Abort." << endl; exit(1);
  }
  // this is to load the size of the position array
  getline(in_fh, line);
  BWTIDX gap_size = (BWTIDX) stol(line);
  if(gap_size != (BWTIDX) (size_ / BWT_FM_GAP) + 1)  {
    cout << "BWT::LoadFMIndex: size of FM index does not match the text, abort." << endl; exit(1);
  }
  // loading the files line by line
  unordered_map<char, bool> loaded_char;
  while(getline(in_fh, line))  {
    vector<int> range;
    range.push_back(-1);
    for(int i = 0; i < line.length(); ++ i) {
      if(line[i] == ':') range.push_back(i);
    }
    char c = line[0];
    if(!alphabet_.IsValid(c) && c != DELIM) {
      cout << "BWT::LoadFMIndex: invalid character found in FM index: " << c << "; abort." << endl; exit(1);
    }
    if(gap_size != range.size() - 2) {
      cout << "BWT::LoadFMIndex: FM index size does not match text length; abort." << endl; exit(1);
    }
    fm_index_[c] = new BWTIDX [gap_size];
    for(int i = 1; i < range.size() - 1; ++ i) {
      (fm_index_[c])[i - 1] = 
          (BWTIDX) stol(line.substr(range[i] + 1, range[i + 1] - range[i] - 1));
    }
    loaded_char[c] = true;
  }
  // check if all characters have been loaded
  for(int i = 0; i < alphabet_.GetSize(); ++ i) {
    char c = alphabet_.GetInvCharMap(i);
    if(loaded_char.find(c) == loaded_char.end()) {
      cout << "BWT::LoadFMIndex: character not found in FM index: " << c << "; abort." << endl; exit(1);
    }
  }
  if(loaded_char.find(DELIM) == loaded_char.end()) {
    cout << "BWT::LoadFMIndex: DELIM not found in FM index; abort." << endl; exit(1);
  }
  in_fh.close();
  return;
}

std::pair<BWTIDX, BWTIDX> BWT::UpdateRange(
  const char c, const std::pair<BWTIDX, BWTIDX> &range
)  {
  int c_id = alphabet_.GetCharMap(c);  // the ID of the character in the alphabet
  if(c != DELIM && c_id < 0)  {
    // make sure that the character is in the alphabet
    cerr << "BWT::UpdateRange: Unknown character being searched against the database: " << c << ", abort." << endl;
    exit(1);
  }
  // calculate the first occurrence
  BWTIDX occbegin = CountOccurrence(c, range.first);
  BWTIDX occend = CountOccurrence(c, range.second);
  // compute the new range
  pair<BWTIDX, BWTIDX> update_range;
  update_range.first = acc_freq_[c_id + 1] + occbegin;
  update_range.second = acc_freq_[c_id + 1] + occend;
  
  cerr << "UpdateRange: " << c << "  " << acc_freq_[c_id + 1] << " " << CountOccurrence(c, range.first) << "  " << CountOccurrence(c, range.second) << endl;
  if(range.second - range.first < 100)  {
    cerr << string((char*) bwt_).substr(range.first, range.second - range.first + 1) << endl;
  }
  return update_range;  
}

void BWT::UpdateRange(
    const char c, const BWTIDX old_begin, const BWTIDX old_end, 
    BWTIDX &new_begin, BWTIDX &new_end
) {
  int c_id = alphabet_.GetCharMap(c);  // the ID of the character in the alphabet
  if(c != DELIM && c_id < 0)  {
    // make sure that the character is in the alphabet
    cerr << "BWT::UpdateRange: Unknown character being searched against the database: " << c << ", abort." << endl;
    exit(1);
  }
  // calculate the first occurrence
  new_begin = acc_freq_[c_id + 1] + CountOccurrence(c, old_begin);
  new_end = acc_freq_[c_id + 1] + CountOccurrence(c, old_end);
  // cerr << "UpdateRange: " << acc_freq_[c_id + 1] << " " << CountOccurrence(c, old_begin) << "  " << CountOccurrence(c, old_end) << endl;
  return;  
}

BWTIDX BWT::CountLexicoLess(const char c, const BWTIDX pos) {
  BWTIDX occ = CountOccurrence(DELIM, pos);
  //cout << "CountLexicoLess: " << occ << endl;
  for(int i = 0; i < alphabet_.GetCharMap(c); ++ i) { 
    occ += CountOccurrence(alphabet_.GetInvCharMap(i), pos);
    //cout << "CountLexicoLess: " << alphabet_.GetInvCharMap(i) << "  " << occ << endl;
  }
  //cout << "End of function CountLexicoLess" << endl;
  return occ;
}
