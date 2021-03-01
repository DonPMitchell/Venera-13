#pragma once
//
//  Simple routines for writing BMP format image files with 8-bit channels
//  D. P. Mitchell  04/10/2003.
//

extern int      ReadTIFF(char *szName, unsigned short rgnData[][252]);
extern int      WriteTIFF(char *szName, int nLength, unsigned short rgnData[][252]);
extern int      WriteHeaderBMP(HANDLE hFile, int nWidth, int nHeight, int nChannels);
extern int      WriteFloatImage(float *pfImage, char *szName, int nHigh, int nWide, double fGamma = 2.2);
extern int      WriteShortImage(unsigned short *pnImage, char *szName, int nHigh, int nWide);
extern int      WriteCharImage(char *pnImage, char *szName, int nHigh, int nWide);
extern int      WriteLongImage(unsigned *pnImage, char *szName, int nHigh, int nWide);
extern int      WriteHeaderWAVE(HANDLE hFile, int nSamples, int nSampleRate);
extern unsigned ColorPixel(int nRed, int nGreen, int nBlue);
