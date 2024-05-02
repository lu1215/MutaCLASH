#include "output_printer.h"

using namespace std;

void OutputPrinter::OutputPrinterCLAN(
  char const* const* ref_header, char const* const* ref_seq,
  char const* const* ref_header_aux, char const* const* ref_seq_aux,
  char const* const* header, char const* const* seq,          // the header and the sequence of the read set
  const int begin, const int size, const int upper,   // the begin and the number of sequences to map in this batch
  std::vector< std::vector<FragType> > &result,      // the output 
  std::ofstream &out_fh               // the output file handle
) {
  out_fh << "#read_id\tsolution_id\tread_mapped_begin\tread_mapped_end\tread_length\tmapped_locations" << endl;
  for(int i = 0; i < result.size() && begin + i < upper; ++ i) {
    //cerr << i << ": " << result.size() << "  " << header[begin + i] << endl;
    for(int j = 0; j < result[i].size(); ++ j) {
      //cerr << " " << j << ": " << result[i].size() << endl; 
      ++ result[i][j].q_begin; ++ result[i][j].q_end;
      out_fh << header[begin + i] << "\t" << result[i][j].sol_id << "\t" << result[i][j].q_begin << "\t" << result[i][j].q_end << "\t" << strlen(seq[begin + i]) << "\t";
      //cerr << header[begin + i] << "\t" << result[i][j].sol_id << "\t" << result[i][j].q_begin << "\t" << result[i][j].q_end << "\t" << strlen(seq[begin + i]) << endl;
      //cerr << "condition check: " << result[i][j].db << "  " << (ref_header_aux != NULL) << endl;
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        // convert 0-based index system to 1-based index system
        ++ result[i][j].targets[k].begin; ++ result[i][j].targets[k].end;
        //cerr << "here1:  " << result[i][j].targets[k].strand << "  " << endl;
        if(result[i][j].targets[k].strand)  {
          //cerr << "check: " << ref_header[1619] << endl;
          //cerr << (long long int) result[i][j].targets[k].begin << endl;
          //cerr << (long long int) result[i][j].targets[k].end << endl;
          if(result[i][j].db)   {
            out_fh << ref_header[result[i][j].targets[k].sid] << ":" << (long long int) result[i][j].targets[k].begin << "-" << (long long int) result[i][j].targets[k].end;
          } else if(!result[i][j].db && ref_header_aux != NULL)    {
            out_fh << ref_header_aux[result[i][j].targets[k].sid] << ":" << (long long int) result[i][j].targets[k].begin << "-" << (long long int) result[i][j].targets[k].end;   
          }
        } else  {
          if(result[i][j].db)  {
            out_fh << ref_header[result[i][j].targets[k].sid] << ":" << (long long int) result[i][j].targets[k].end << "-" << (long long int) result[i][j].targets[k].begin;
          } else if(!result[i][j].db && ref_header_aux != NULL)  {
            out_fh << ref_header_aux[result[i][j].targets[k].sid] << ":" << (long long int) result[i][j].targets[k].end << "-" << (long long int) result[i][j].targets[k].begin; 
          }
        }
        //cerr << "here2:  " << result[i][j].targets[k].strand << "  " << endl;
        if(k < result[i][j].targets.size() - 1) out_fh << ";";
        //delete [] result[i][j].targets[k].cigar;
      }
      //cerr << "finished output" << endl;
      out_fh << "\t";
      
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        delete [] result[i][j].targets[k].cigar;
      }
      //cerr << "finished delete" << endl;
      out_fh << endl;
    }
  }
  return;
}

void OutputPrinter::ParseCIGAR(const char* cigar, int &align_len, int &num_gap)  {
  int i = 0, n = strlen(cigar);
  int num = 0;
  align_len = 0; num_gap = 0;
  while(i < n && cigar[i] != '\0') {
    if(cigar[i] >= (int) '0' && cigar[i] <= (int) '9')  {
      // number, update number values
      num = num * 10 + ((int) cigar[i] - 48);
    } else if(cigar[i] == 'M') {
      // match case, increase alignment length, do not increase gap
      align_len += num;
      num = 0;
    } else if(cigar[i] == 'I' || cigar[i] == 'D') {
      // insertion or deletion; increase both alignment length and gap size
      align_len += num; num_gap += num;
      num = 0;
    } else  {
      cerr << "CLAN Error: ParseCIGAR: unsupported CIGAR character \'" << cigar[i] << "\'. Parsing may be incorrect."<< endl;
      return;
    }
    ++ i;
  }
  return;
}

void OutputPrinter::OutputPrinterBLAST(
  char const* const* ref_header, char const* const* ref_seq,
  char const* const* header, char const* const* seq,          // the header and the sequence of the read set
  const int begin, const int size, const int upper,   // the begin and the number of sequences to map in this batch
  std::vector< std::vector<FragType> > &result,       // the output 
  const std::string &db_name, const long long int db_size,   // the name and size of the database
  std::ofstream &out_fh               // the output file handle
)  {
  for(int i = 0; i < result.size() && begin + i < upper; ++ i) {
    out_fh << "# CLAN v0.02" << endl;
    out_fh << "# Query: " << header[begin + i] << endl;
    out_fh << "# Database: " << db_name << endl;
    if(result[i].size() <= 0) {
      out_fh << "# 0 hits found" << endl;
      continue;
    } else  {
      out_fh << "# Fields: query acc.ver, subject acc.ver, \% identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score" << endl;
      out_fh << "# " << result[i].size() << " hits found" << endl;
    }
    for(int j = 0; j < result[i].size(); ++ j) {
      
      ++ result[i][j].q_begin; ++ result[i][j].q_end;
      // calculate bit score
      double bit_score = (NT_BLAST_LAMBDA * result[i][j].score - log(NT_BLAST_K)) / log(2);
      double evalue = log10(strlen(seq[begin + i])) + log10(db_size) - bit_score * log10(2);
      int e_exponent = floor(evalue);
      double e_const = pow(10, ((float) evalue) - e_exponent);
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        out_fh << header[begin + i] << "\t" << ref_header[result[i][j].targets[k].sid] << "\t";  // query and target accesion
        int align_len, num_gap;
        ParseCIGAR(result[i][j].targets[k].cigar, align_len, num_gap);
        double identity = (double) 100 * (align_len - num_gap - result[i][j].targets[k].num_mismatch) / align_len;
        out_fh << std::fixed << std::setprecision(3) << identity << "\t";
        out_fh << align_len << "\t" << result[i][j].targets[k].num_mismatch << "\t" << num_gap << "\t";
        out_fh << result[i][j].q_begin << "\t" << result[i][j].q_end << "\t";
        // convert 0-based index system to 1-based index system
        ++ result[i][j].targets[k].begin; ++ result[i][j].targets[k].end;
        if(result[i][j].targets[k].strand)  {
          out_fh << (long long int) result[i][j].targets[k].begin << "\t" << (long long int) result[i][j].targets[k].end << "\t";
        } else  {
          out_fh << (long long int) result[i][j].targets[k].end << "\t" << (long long int) result[i][j].targets[k].begin << "\t";
        }
        
        out_fh << std::fixed << std::setprecision(2) << e_const << "e" << e_exponent << "\t";
        out_fh << std::setprecision(2) << bit_score;
        out_fh << endl;
        //delete [] result[i][j].targets[k].cigar;
      }
      
      //out_fh << "\t";
      
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        delete [] result[i][j].targets[k].cigar;
      }
      
    }
  }
  return;
}

void OutputPrinter::OutputPrinterSAM(
  char const* const* ref_header, char const* const* ref_seq,
  char const* const* header, char const* const* seq,          // the header and the sequence of the read set
  const int begin, const int size, const int upper,   // the begin and the number of sequences to map in this batch
  std::vector< std::vector<FragType> > &result,       // the output 
  const long long int db_size,   // the name and size of the database
  std::ofstream &out_fh               // the output file handle
)  {
  for(int i = 0; i < result.size() && begin + i < upper; ++ i) {
    for(int j = 0; j < result[i].size(); ++ j) {
      
      ++ result[i][j].q_begin; ++ result[i][j].q_end;
      // calculate bit score
      double bit_score = (NT_BLAST_LAMBDA * result[i][j].score - log(NT_BLAST_K)) / log(2);
      
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        out_fh << header[begin + i] << "\t" << (result[i][j].targets[k].strand ? 0 : 16) << "\t";
        out_fh << ref_header[result[i][j].targets[k].sid] << "\t" << result[i][j].targets[k].begin << "\t" << (int) bit_score << "\t";
        out_fh << result[i][j].q_begin << "S" << result[i][j].targets[k].cigar << strlen(seq[begin + i]) - result[i][j].q_end << "S\t";
        out_fh << "*\t*\t*\t*\t" << endl;
      }
      
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        delete [] result[i][j].targets[k].cigar;
      }
      
    }
  }
  return;
}

void OutputPrinter::OutputEncoded(
  char const* const* ref_seq, const int num_refs,
  const int begin, const int size, const int upper,
  std::vector< std::vector<FragType> > &result, 
  std::ofstream &out_fh  
) {
  
  //for(int i = 0; i < result.size() && begin + i < upper; ++ i) {
  //  cerr << begin + i << endl;
  //  cerr <<
  //}
  
////////////////////////////////////////
  //vector<BWTIDX> ref_len; ref_len.resize(num_refs);
  //for(int k = 0; k < num_refs; ++ k) {
  //  ref_len[k] = string(ref_seq[k]).length();
  //}
////////////////////////////////////////

  //return;
  
  long long int total_hits = 0;
  for(int i = 0; i < result.size() && begin + i < upper; ++ i) {
    for(int j = 0; j < result[i].size(); ++ j) {
      total_hits += result[i][j].targets.size();
    }
  }
  
  out_fh.write((char*) &total_hits, sizeof(long long int));
  if(total_hits == 0){
    return;
  }
  
  for(int i = 0; i < result.size() && begin + i < upper; ++ i) {
    int rid = begin + i;
    int num_maps = result[i].size();
    
    // do not write anything if the read has no mapping
    if(num_maps <= 0) continue;
    
    out_fh.write((char*) &rid, sizeof(int));
    out_fh.write((char*) &num_maps, sizeof(int));
    
    for(int j = 0; j < result[i].size(); ++ j) {
      
      out_fh.write((char*) &result[i][j].q_begin, sizeof(BWTINT));
      out_fh.write((char*) &result[i][j].q_end, sizeof(BWTINT));
      out_fh.write((char*) &result[i][j].sol_id, sizeof(int));
      out_fh.write((char*) &result[i][j].score, sizeof(int));
      out_fh.write((char*) &result[i][j].db, sizeof(bool));
      //cerr << "write db condition:  " <<  result[i][j].db << endl;
      int target_count = result[i][j].targets.size();
      out_fh.write((char*) &target_count, sizeof(int));
      
      if(target_count > 1000) {
        cerr << "Output error: super large number of targets: " << target_count << "  " << rid << endl;
      }
      
      for(int k = 0; k < result[i][j].targets.size(); ++ k) {
        out_fh.write((char*) &result[i][j].targets[k].sid, sizeof(int));
        out_fh.write((char*) &result[i][j].targets[k].begin, sizeof(BWTIDX));
        out_fh.write((char*) &result[i][j].targets[k].end, sizeof(BWTIDX));
        //if(result[i][j].targets[k].sid >= num_refs)  {
        //  cerr << "ref ID error!!!" << result[i][j].targets[k].sid << "  " << num_refs << endl;
        //}
        //if(result[i][j].targets[k].end >= ref_len[result[i][j].targets[k].sid])  {
        //  cerr << "Write error!!! " << rid << " " << result[i][j].targets[k].end << " " << result[i][j].targets[k].sid << endl;
        //}
             
        out_fh.write((char*) &result[i][j].targets[k].strand, sizeof(bool));
        out_fh.write((char*) &result[i][j].targets[k].num_mismatch, sizeof(int));
        //cerr <<rid<<"    "<<num_maps<<"    "<<result[i][j].q_begin << "    "<<result[i][j].targets.size()<<"   "<<result[i][j].targets[k].sid <<"    " <<result[i][j].targets[k].begin << "   " <<result[i][j].targets[k].end  <<"    " << result[i][j].targets[k].num_mismatch << endl;
        
        // DEBUG
        //cerr << i << "  " << j << " " << k << endl;
        //cerr << result[i][j].targets[k].sid << "    " << result[i][j].targets[k].begin << " " << result[i][j].targets[k].end << endl;
        //cerr << result[i][j].targets[k].strand << "    " << result[i][j].targets[k].num_mismatch << endl;
        //cerr << result[i][j].targets[k].cigar << endl;
        
        int cigar_len = strlen(result[i][j].targets[k].cigar);
        if(cigar_len < 2)  {
          cerr << "Output error: super short cigar: " << cigar_len << "  " << result[i][j].targets[k].cigar << "    " << result[i][j].targets[k].begin << " " << result[i][j].targets[k].end << endl;
        }
        
        out_fh.write((char*) &cigar_len, sizeof(int));
        out_fh.write(result[i][j].targets[k].cigar, cigar_len);
        //cerr << result[i][j].targets[k].cigar << endl;
        //delete [] result[i][j].targets[k].cigar;
      }
    }
  }
  return;
}

void ExitWithMessage(const char* message, const int exit_code)  {
  cerr << message << endl;
  exit(exit_code);
}

// TODO: need to verify the consistency of the reference and sequence database before processing
void OutputPrinter::ReadEncoded(
  const int max_entry, const std::string &file_name,
  std::vector< std::vector<FragType> > &result
) {

  ifstream in_fh;
  in_fh.open(file_name, ios::in);
  if(!in_fh.is_open()) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: Cannot open file. Abort.", 1);
  result.resize(max_entry);
  int num_chunks; in_fh.read((char *) &num_chunks, sizeof(int));
  int current_chunk = 0;
  //cerr << num_chunks << endl;
  while(!in_fh.eof()) {
    long long int total_hits; in_fh.read((char *) &total_hits, sizeof(long long int));
    //if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
    long long int current_hits = 0;
    if(total_hits==0){
      continue;
    }
    while(!in_fh.eof()) {
      int rid;  in_fh.read((char *) &rid, sizeof(int));
      if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
      int num_maps; in_fh.read((char *) &num_maps, sizeof(int));
      if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
      result[rid].resize(num_maps);
      //cerr << "rid: " << rid << endl; 
      //cerr << "num maps:  " << num_maps << endl;
      for(int i = 0; i < num_maps; ++ i) {
        in_fh.read((char*) &result[rid][i].q_begin, sizeof(BWTINT));
        if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
        in_fh.read((char*) &result[rid][i].q_end, sizeof(BWTINT));
        if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
        in_fh.read((char*) &result[rid][i].sol_id, sizeof(int));
        if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
        in_fh.read((char*) &result[rid][i].score, sizeof(int));
        if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
        in_fh.read((char*) &result[rid][i].db, sizeof(bool));
        //cerr << "read db condition: " << result[rid][i].db << endl;
        if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
        int num_targets; in_fh.read((char *) &num_targets, sizeof(int));
        if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
        result[rid][i].targets.resize(num_targets);
       // cerr << "num targets: " << num_targets << endl;
        for(int j = 0; j < num_targets; ++ j) {
          in_fh.read((char*) &result[rid][i].targets[j].sid, sizeof(int));
          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          in_fh.read((char*) &result[rid][i].targets[j].begin, sizeof(BWTIDX));
          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          in_fh.read((char*) &result[rid][i].targets[j].end, sizeof(BWTIDX));
          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          //if(result[rid][i].targets[j].end >= ref_len[result[rid][i].targets[j].sid])  {
          //  cerr << "Read error!!!  " << rid << " " << result[rid][i].targets[j].end << " " << result[rid][i].targets[j].sid << endl;
          //}
          
          
          in_fh.read((char*) &result[rid][i].targets[j].strand, sizeof(bool));
          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          in_fh.read((char*) &result[rid][i].targets[j].num_mismatch, sizeof(int));


          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          int cigar_len;
          in_fh.read((char*) &cigar_len, sizeof(int));
          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          result[rid][i].targets[j].cigar = new char [cigar_len + 1];
          in_fh.read(result[rid][i].targets[j].cigar, cigar_len);
          if(!in_fh) ExitWithMessage("CLAN:OutputPrinter::ReadEncoded: file could be corrupted. Abort.", 1);
          result[rid][i].targets[j].cigar[cigar_len] = '\0';
          //cerr << "CIGAR string:" << result[rid][i].targets[j].cigar << endl;
          ++ current_hits;
        }
        //cerr << "current vs total hits: " << current_hits << "  " << total_hits << endl;
      }
      if(current_hits >= total_hits) break;
    }
    //cerr << num_chunks << endl;
    ++ current_chunk; 
    //cerr << "current chunk: " << current_chunk << " " << num_chunks << endl;
    if(current_chunk >= num_chunks-1) break;
  }
  return;
}


void OutputPrinter::OutputAnnotationTable(
  const bool print_simplex, char const* const* ref_header, 
  std::vector<GPosType> &islands, std::vector<DuplexType> &duplex_annot, 
  std::ofstream &out_fh
) {
  // TODO: add option to output RPKM or raw reads

  char sr = '+'; int stl = 0;
  if(!print_simplex)  {
    // print header information
    out_fh << "#P1_ID\tP1_begin\tP1_end\tP1_strand\tP1_reads\tP1_annotation\t";
    out_fh << "P2_ID\tP2_begin\tP2_end\tP2_strand\tP2_reads\tP2_annotation\t";
    out_fh << "Duplex_reads\tDimer_MFE\tStack_MFE\tStack_len\n";  
    for(int i = 0; i < duplex_annot.size(); ++ i) {
    
      int p = duplex_annot[i].island_ID1, q = duplex_annot[i].island_ID2;
    
      //if(p == 250 && q == 2996)  {
      //  cerr << "I am printing out" << endl;
      //  cerr << islands[p].begin << " " << islands[p].end << endl;
      //}
      
      //cerr << islands[p].coverage << "  " << islands[q].coverage << endl;
      // info regarding the first partner
      out_fh << ref_header[islands[p].sid] << "\t" << islands[p].begin << "\t" << islands[p].end << "\t";
      sr = duplex_annot[i].strand_1 ? '+' : '-';
      out_fh << sr << "\t" << islands[p].coverage << "\t" << "NA" << "\t";
      // info regarding the second partner
      out_fh << ref_header[islands[q].sid] << "\t" << islands[q].begin << "\t" << islands[q].end << "\t";
      sr = duplex_annot[i].strand_2 ? '+' : '-';
      out_fh << sr << "\t" << islands[q].coverage << "\t" << "NA" << "\t";
      // info regarding the duplex
      out_fh << duplex_annot[i].info.coverage << "\t" << duplex_annot[i].info.dimer_energy << "\t";
      out_fh << duplex_annot[i].info.stack_info.energy << "\t";
      stl = duplex_annot[i].info.stack_info.end - duplex_annot[i].info.stack_info.begin + 1;
      out_fh << stl << endl;
    }
  } else  {
    // only print the summary of the detected islands
    out_fh << "#P1_ID\tP1_begin\tP1_end\tP1_strand\tP1_reads\tP1_annotation\n";
    for(int i = 0; i < islands.size(); ++ i) {
      out_fh << ref_header[islands[i].sid] << "\t" << islands[i].begin << "\t" << islands[i].end << "\t";
      sr = islands[i].strand ? '+' : '-';
      out_fh << sr << "\t" << islands[i].coverage << "\t" << "NA" << endl;
    }
  }
  return;
}
