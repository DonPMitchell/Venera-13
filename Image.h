#pragma once
//
//  Simple image-file routines (Windows BMP format)
//  D. P. Mitchell  03/31/2003.
//
extern int  WriteHeaderBMP(HANDLE hFile, int nWidth, int nHeight, int nChannels);
