#include "stdafx.h"
//
//  Shift the bit stream to correct for dropped bits in transmission.
//  D. P. Mitchell  05/19/2003.
//
//  The "bit stream" consists of 9 pixel bits (from low to high order) and
//  a 10th bit creating odd parity.  That stream is shifted, and the missing
//  bit is reconstructed from parity.  Sections of image that may look like
//  noise are sometimes just shifted so the parity bit is in a high order position.
//
//  Venera 9 and 10 used 6 bit pixels in a 7-bit block.
//
#include "Venus.h"

static int NBLOCK = 10;

static inline int
Parity(int n)
{
    n = n ^ (n >> 16);
    n = n ^ (n >> 8);
    n = n ^ (n >> 4);
    n = n ^ (n >> 2);
    n = n ^ (n >> 1);
    n ^= 1;             // Venera uses odd parity
    return n & 1;
}

void
BitShift(unsigned short rgnData[], int nShift, int nLength)
{
    char *pMemory, *pData, *pShifted;
    int i, iBit, nPixel, nBad;

    Assert(abs(nShift) < 999);
    pMemory = new char[NBLOCK*(nLength + abs(nShift) + 1)];
    if (nShift >= 0) {
        pData = pMemory;
    } else {
        pData = pMemory - nShift;
    }
    //
    //  Expand pixels into blocks of 10 bits (7 bits for Venera 9 and 10)
    //
    for (i = 0; i < nLength; i++) {
        nPixel = rgnData[i] >> (17 - NBLOCK);
        for (iBit = 0; iBit < NBLOCK; iBit++)
            pData[NBLOCK*i + iBit] = (nPixel >> iBit) & 1;  // parity bit == 0
    }
    pShifted = pData + nShift;
    //
    //  Pack bits into pixels, restoring missing bit via parity
    //
    for (i = 0; i < nLength; i++) {
        nPixel = 0;
        for (iBit = 0; iBit < NBLOCK; iBit++)
            nPixel |= pShifted[NBLOCK*i + iBit] << iBit;
        nPixel |= Parity(nPixel) << ((1000*NBLOCK - nShift - 1) % NBLOCK);
        rgnData[i] = nPixel << (17 - NBLOCK);
    }
    //
    //  Garbage is left where bits were shifted away
    //
    nBad = (abs(nShift) + NBLOCK - 1)/NBLOCK;
    if (nShift < 0)
        for (i = 0; i < nBad; i++)
            rgnData[i] = BADPIXEL;
    else
        for (i = 0; i < nBad; i++)
            rgnData[nLength - i - 1] = BADPIXEL;
    delete [] pMemory;
}

void
BitShiftPixels(unsigned short rgnData[MAXLINES][252], int nShift,
               int jFirst, int jLimit, int iFirst, int iLimit)
{
    NBLOCK = 10;
    BitShift(&rgnData[jFirst][iFirst], nShift, 
                    252*(jLimit - jFirst - 1) + iLimit - iFirst);
}
//
//  Venera 9 and 10 seem to use a 7-bit block
//
static void
BitShiftVenera9(unsigned short rgnData[], int nShift, int nLength)
{
    char *pMemory, *pData, *pShifted;
    int i, iBit, nPixel, nBad;

    NBLOCK = 7;             // NBLOCK = 6, no parity did not work
    Assert(abs(nShift) < 999);
    pMemory = new char[NBLOCK*(nLength + abs(nShift) + 1)];
    if (nShift >= 0) {
        pData = pMemory;
    } else {
        pData = pMemory - nShift;
    }
    //
    //  Expand pixels into blocks of 10 bits (7 bits for Venera 9 and 10)
    //
    for (i = 0; i < nLength; i++) {
        nPixel = rgnData[i] >> 10;
        //nPixel = (rgnData[i] >> 9) & 0x007E;    // try parity bit on other end
        for (iBit = 0; iBit < NBLOCK; iBit++)
            pData[NBLOCK*i + iBit] = (nPixel >> iBit) & 1;  // parity bit == 0
    }
    pShifted = pData + nShift;
    //
    //  Pack bits into pixels, restoring missing bit via parity
    //
    for (i = 0; i < nLength; i++) {
        nPixel = 0;
        for (iBit = 0; iBit < NBLOCK; iBit++)
            nPixel |= pShifted[NBLOCK*i + iBit] << iBit;
        nPixel |= (Parity(nPixel)^0) << ((1000*NBLOCK - nShift - 1) % NBLOCK);
        rgnData[i] = nPixel << 10;
        //nPixel |= (Parity(nPixel)^0) << ((1000*NBLOCK + nShift) % NBLOCK);
        //rgnData[i] = (nPixel & 0x007E) << 9;
    }
    //
    //  Garbage is left where bits were shifted away
    //
    nBad = (abs(nShift) + NBLOCK - 1)/NBLOCK;
    if (nShift < 0)
        for (i = 0; i < nBad; i++)
            rgnData[i] = BADPIXEL;
    else
        for (i = 0; i < nBad; i++)
            rgnData[nLength - i - 1] = BADPIXEL;
    delete [] pMemory;
}

void
BitShiftVenera9Pixels(unsigned short rgnData[512][128], int nShift,
               int jFirst, int jLimit, int iFirst, int iLimit)
{
    BitShiftVenera9(&rgnData[jFirst][iFirst], nShift, 
                    252*(jLimit - jFirst - 1) + iLimit - iFirst);
}