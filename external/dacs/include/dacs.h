/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <nlehmann@dcc.uchile.cl> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Nicolás Lehmann
 * ----------------------------------------------------------------------------
 *
 * This file contain declaration of function implementing DAC to use inside C++.
 * It also contains implementations of the LoadFT and SaveFT functions using
 * file streams.
 */

#ifndef INCLUDE_DACS_H_
#define INCLUDE_DACS_H_

extern "C" {
typedef unsigned int uint;

struct sFTRep;
typedef struct sFTRep FTRep;

FTRep* createFT(uint *list, uint listLength);
uint accessFT(FTRep * listRep, uint param);
uint * decompressFT(FTRep * listRep, uint n);
void destroyFT(FTRep * listRep);
}
#include <fstream>


/**
 * Saves DAC to file
 *
 * @param out Output stream.
 * @param rep DAC representation.
 */
void SaveFT(std::ostream& out, FTRep *rep);
/**
 * Loads DAC from file
 *
 * @param in Input stream.
 * @return Pointer to representation. The caller must take the responsibility
 * to free the memory with destroyFT.
 */
FTRep *LoadFT(std::istream& in);
bool equalsFT(FTRep *lhs, FTRep *rhs);
uint getListLength(FTRep *ftrep);
#endif  // INCLUDE_DACS_H_
