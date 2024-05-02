#ifndef _MISC_H_
#define _MISC_H_

// common small utility functions

#include <string>

namespace misc  {

  std::string RevComplement(const std::string &s) {
    // assuming all capitalized chars
    char *t = new char [s.length() + 1];
    int k = 0;
    for(int i = s.length() - 1; i >= 0; -- i) {
      switch(s[i])  {
        case 'A':
          t[k ++] = 'T'; break;
        case 'C':
          t[k ++] = 'G'; break;
        case 'G':
          t[k ++] = 'C'; break;
        case 'T':
          t[k ++] = 'A'; break;
        case 'N':
          t[k ++] = 'N'; break;
        default:
          std::cerr << "MISC::RevComplement: unrecognized char: " << s[i] << "  " << s << std::endl;
          exit(1);
          t[k ++] = s[i]; 
          break;
      }
    }
    t[k] = '\0';
    std::string st(t);
    delete [] t;
    return std::string(st);
  }

  double MyTime (void)  {
      int flag;
      clockid_t cid = CLOCK_REALTIME; // CLOCK_MONOTONE might be better
      timespec tp;
      double timing;
	
      flag = clock_gettime(cid, &tp);
      if (flag == 0) timing = tp.tv_sec + 1.0e-9*tp.tv_nsec;
      else           timing = -17.0;         // If timer failed, return non-valid time
	  
      return(timing);
  }

  void PrintElapsed( double s, double e, const char *task ) {
	  double elapsed = e - s ;
	  printf ("[%02.0f:%02.0f:%05.2f]\t%s\n", 
	  		floor(elapsed/3600.0), 
	  		floor(fmod(elapsed,3600.0)/60.0), 
	  		fmod(elapsed,60.0),
	  		task);
	  return;
  }

}

#endif
