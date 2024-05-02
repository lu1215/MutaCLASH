#ifndef _FREE_ENERGY_H_
#define _FREE_ENERGY_H_

#ifndef SCALE_TO_INT
#define SCALE_TO_INT 10
#endif

#include <unordered_map>
#include <string>

namespace free_energy
{

  static std::unordered_map<std::string, int> basepair_map = 
  {
    {"AU", 0}, 
    {"CG", 1}, 
    {"GC", 2}, 
    {"GU", 3}, 
    {"UA", 4}, 
    {"UG", 5}
  };

  static int stack[6][6] =
  {
    //AU   CG   GC   GU   UA   UG 
    { -9, -22, -21,  -6, -11, -14}, // AU
    {-21, -33, -24, -14, -21, -21}, // CG
    {-24, -34, -33, -15, -22, -25}, // GC
    {-13, -25, -21,  -5, -14,  13}, // GU
    {-13, -24, -21, -10,  -9, -13}, // UA
    {-10, -15, -14,   3,  -6,  -5}  // UG
  };

  static int bulge[30] = 
  {
    38, 28, 32, 36, 40, 44,
    46, 47, 48, 49, 50, 51,
    52, 53, 54, 54, 55, 55,
    56, 57, 57, 58, 58, 58,
    59, 59, 60, 60, 60, 61
  };

  static int internal[30] = 
  {
    999, 999, 999,  11,  20,  20,
     21,  23,  24,  25,  26,  27,
     28,  29,  29,  30,  31,  31,
     32,  33,  33,  34,  34,  35,
     35,  35,  36,  36,  37,  37
  };


  static int internal_closure_AUGU = 7;
  static int internal_asymmetric = 6;



/*
static int terminal_mismatch[4][4] = 
{
  // A    C    G    U
  {999, 999, -10, 999}, // A
  {999, 999, 999, 999}, // C
  { -8, 999, -12, 999}, // G
  {999, 999, 999,  -7}  // U
}
*/
}

#endif
