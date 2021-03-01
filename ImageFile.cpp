#include "stdafx.h"
#include "Venus.h"
#pragma intrinsic(pow)
//
//  Read the Russian Venera TIFF files.  Can't really parse any other TIFF files.
//
struct IFD {
    IFD(int n1, int n2, int n3, int n4) : nTag(n1), nType(n2), nCount(n3), nValue(n4) {}
    unsigned short  nTag;
    unsigned short  nType;
    unsigned int    nCount;
    unsigned int    nValue;
};
struct Tiff {
    unsigned short  nByteOrder;
    unsigned short  nMagic;
    unsigned int    nOffsetIFD;
};

static char     s_rgchTiffHeader[202];
static Tiff     *s_pHeader = (Tiff *)(s_rgchTiffHeader);
static short    *s_pnTags = (short *)(s_rgchTiffHeader + 8);
static IFD      *s_rgTags = (IFD *)(s_rgchTiffHeader + 10);
static char     *s_szSoftware = (s_rgchTiffHeader + 170);
static unsigned *s_pnNextIFD = (unsigned *)(s_rgchTiffHeader + 166);

int
ReadTIFF(char *szName, unsigned short rgnData[][252])
{
    HANDLE hFile;
    DWORD nBytesRead;
    int nLength;

    hFile = CreateFile(szName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                                   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                                   NULL);
    Assert(hFile != INVALID_HANDLE_VALUE);
    ReadFile(hFile, s_rgchTiffHeader, 202, &nBytesRead, NULL);
    //for (i = 0; i < *s_pnTags; i++)
    //    printf("%2d: %3d %1d %2d %d\n", i, s_rgTags[i].nTag, s_rgTags[i].nType,
    //            s_rgTags[i].nCount, s_rgTags[i].nValue);
    Assert(s_pHeader->nByteOrder == 0x4949);
    Assert(s_pHeader->nMagic == 42);                                // The Douglas Adams number
    Assert(*s_pnTags == 13);
    Assert(s_rgTags[1].nTag == 256 && s_rgTags[1].nValue == 252);   // image width
    Assert(s_rgTags[3].nTag == 258 && s_rgTags[3].nValue == 16);    // bits per pixel
    Assert(s_rgTags[2].nTag == 257);                                // image length
    nLength = s_rgTags[2].nValue;
    Assert(s_rgTags[6].nTag == 273 && s_rgTags[6].nValue == 202);   // image begins
    ReadFile(hFile, rgnData, 252*2*nLength, &nBytesRead, NULL);
    Assert(nBytesRead == 252*2*nLength);
    CloseHandle(hFile);
    return nLength;
}

int
WriteTIFF(char *szName, int nLength, unsigned short rgnData[][252])
{
    HANDLE hFile;
    DWORD nBytesWritten;

    //
    //  Build TIFF header, based on template used by Russian software
    //
    s_pHeader->nByteOrder = 0x4949;
    s_pHeader->nMagic = 42;
    s_pHeader->nOffsetIFD = 8;
    *s_pnTags = 13;
    s_rgTags[ 0] = IFD(254, 4,   1,       0);           // NewSubfileType
    s_rgTags[ 1] = IFD(256, 4,   1,     252);           // Image Width
    s_rgTags[ 2] = IFD(257, 4,   1,     nLength);       // Image Length
    s_rgTags[ 3] = IFD(258, 3,   1,      16);           // Bits Per Pixel
    s_rgTags[ 4] = IFD(259, 3,   1,       1);           // No Compression
    s_rgTags[ 5] = IFD(262, 3,   1,       1);           // 0 means black
    s_rgTags[ 6] = IFD(273, 4,   1,     202);           // Start of Pixels
    s_rgTags[ 7] = IFD(277, 3,   1,       1);           // Channels per Pixel
    s_rgTags[ 8] = IFD(278, 4,   1,     nLength);       // Rows per Strip
    s_rgTags[ 9] = IFD(279, 4,   1,     504*nLength);   // Bytes per Strip
    s_rgTags[10] = IFD(284, 3,   1,       1);           // Planes Chunky
    s_rgTags[11] = IFD(305, 2,  16,     170);           // Software String
    s_rgTags[12] = IFD(317, 3,   1,       1);           // No Prediction Coding
    *s_pnNextIFD = 0;
    strcpy(s_szSoftware, "Don P. Mitchell Venus Project");
    s_rgTags[11].nCount = unsigned(strlen(s_szSoftware) + 1);
    hFile = CreateFile(szName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
                               FILE_ATTRIBUTE_NORMAL, NULL);
    Assert(hFile != INVALID_HANDLE_VALUE);
    WriteFile(hFile, s_rgchTiffHeader, 202, &nBytesWritten, NULL);
    WriteFile(hFile, rgnData, nLength*504, &nBytesWritten, NULL);
    Assert(nBytesWritten == nLength*504);
    CloseHandle(hFile);
    return nLength;
}

int
WriteHeaderBMP(HANDLE hFile, int nWidth, int nHeight, int nChannels)
{
	unsigned char	bfType[2];
	struct {
		unsigned int	bfSize;
		unsigned short	bfReserved1;
		unsigned short	bfReserved2;
		unsigned int	bfOffBits;

		unsigned int	biSize;
		unsigned int	biWidth;
		unsigned int	biHeight;
		unsigned short	biPlanes;
		unsigned short	biBitCount;
		unsigned int	biCompression;
		unsigned int	biSizeImage;
		unsigned int	biXPelsPerMeter;
		unsigned int	biYPelsPerMeter;
		unsigned int	biClrUsed;
		unsigned int	biClrImportant;
	} H;
    struct {
        unsigned char   b, g, r, a;
    } rgColors[256];
    unsigned long nBytesWritten, nLineSize;
    int i;

    //
    //  Build header structure
    //
    H.biSize = 40;
    H.biWidth = nWidth;
    H.biHeight = nHeight;
    H.biPlanes = 1;
    if (nChannels == 1)
        H.biBitCount = 8;       // Grayscale image, needs a colormap
    else if (nChannels == 3)
        H.biBitCount = 24;      // RGB color image
    else
        return 0;
    H.biCompression = 0;
	H.biSizeImage = 0;
	H.biXPelsPerMeter = 0;
	H.biYPelsPerMeter = 0;
	H.biClrUsed = 0;
	H.biClrImportant = 0;
	bfType[0] = 'B';
	bfType[1] = 'M';
	nLineSize = 4*((H.biWidth*H.biBitCount + 31) >> 5);   // pad to 32-bit boundary
	H.bfSize = 40 + 14 + nHeight*nLineSize;
	H.bfReserved1 = 0;
	H.bfReserved2 = 0;
	H.bfOffBits = 40 + 14;
    //
    //  2 channels would be an FFT magnitude/phase image.
    //  The magnitude image will be written
    //
    if (nChannels <= 2) {
        H.bfSize += sizeof(rgColors);
        H.bfOffBits += sizeof(rgColors);
    }
    //
    //  Write header and colormap
    //
	if (WriteFile(hFile, bfType, 2, &nBytesWritten, NULL) == 0)
        return 0;
    if (WriteFile(hFile, &H, sizeof(H), &nBytesWritten, NULL) == 0)
        return 0;
    if (nChannels <= 2) {
        for (i = 0; i < 256; i++)
            rgColors[i].r = rgColors[i].g = rgColors[i].b = i;
        if (WriteFile(hFile, rgColors, sizeof(rgColors), &nBytesWritten, NULL) == 0)
            return 0;
    }
    return 1;
}
//
//  Write out an array of floating-point numbers as gamma-corrected bmp image file.
//
int
WriteFloatImage(float *pfImage, char *szName, int nHigh, int nWide, double fGamma)
{
    HANDLE hFile;
    char *pch;
    DWORD nBytesWritten;
    int i, j, nLine;
    double f;

    hFile = CreateFile(szName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
                               FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("Could not create %s\n", szName);
        return 0;
    }
    nLine = (nWide + 3) & 0xFFFFFFFC;
    pch = new char[nLine];
    for (i = 0; i < nLine; i++)
        pch[i] = 0;
    WriteHeaderBMP(hFile, nWide, nHigh, 1);
    fGamma = 1.0/fGamma;
    for (j = nHigh - 1; j >= 0; --j) {  // bmp format is upside down
        for (i = 0; i < nWide; i++) {
            f = pfImage[i + j*nWide];
            if (f > 1.0)
                f = 1.0;
            if (f < 0.0)
                f = 0.0;
            f = pow(f, fGamma);
            pch[i] = int(f*255.0 + 0.5);
        }
        WriteFile(hFile, pch, nLine, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
    delete [] pch;
    return 1;
}

int
WriteShortImage(unsigned short *pnImage, char *szName, int nHigh, int nWide)
{
    HANDLE hFile;
    char *pch;
    DWORD nBytesWritten;
    int i, j, nLine, nPixel, nData;

    hFile = CreateFile(szName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                               FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("Could not create %s\n", szName);
        return 0;
    }
    nLine = (nWide + 3) & 0xFFFFFFFC;
    pch = new char[nLine];
    for (i = 0; i < nLine; i++)
        pch[i] = 0;
    if (nHigh < 0) {
        for (j = -nHigh - 1; j >= 0; --j) {
            nData = 0;
            for (i = 0; i < nWide; i++)
                if (pnImage[i + j*nWide] & 0xFFF0)
                    nData++;
            if (nData > 0)
                break;
        }
        nHigh = j + 1;
    }
    WriteHeaderBMP(hFile, nWide, nHigh, 1);
    for (j = nHigh - 1; j >= 0; --j) {  // bmp format is upside down
        for (i = 0; i < nWide; i++) {
            nPixel = pnImage[i+j*nWide];
            if (nPixel & 1)
                pch[i] = 0;     // bad-pixel flag
            else
                pch[i] = nPixel >> 8;
        }
        WriteFile(hFile, pch, nLine, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
    delete [] pch;
    return 1;
}

int
WriteCharImage(char *pnImage, char *szName, int nHigh, int nWide)
{
    HANDLE hFile;
    char *pch;
    DWORD nBytesWritten;
    int i, j, nLine;

    hFile = CreateFile(szName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                               FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("Could not create %s\n", szName);
        return 0;
    }
    nLine = (nWide + 3) & 0xFFFFFFFC;
    pch = new char[nLine];
    for (i = 0; i < nLine; i++)
        pch[i] = 0;
    WriteHeaderBMP(hFile, nWide, nHigh, 1);
    for (j = nHigh - 1; j >= 0; --j) {  // bmp format is upside down
        for (i = 0; i < nWide; i++) {
            pch[i] = pnImage[i+j*nWide];
        }
        WriteFile(hFile, pch, nLine, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
    delete [] pch;
    return 1;
}
//
//  Color rgb images
//
int
WriteLongImage(unsigned *pnImage, char *szName, int nHigh, int nWide)
{
    HANDLE hFile;
    char *pch;
    DWORD nBytesWritten;
    int i, j, nLine;

    hFile = CreateFile(szName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                               FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("Could not create %s\n", szName);
        return 0;
    }
    nLine = (3*nWide + 3) & 0xFFFFFFFC;
    pch = new char[nLine];
    for (i = 0; i < nLine; i++)
        pch[i] = 0;
    WriteHeaderBMP(hFile, nWide, nHigh, 3);
    for (j = nHigh - 1; j >= 0; --j) {  // bmp format is upside down
        for (i = 0; i < nWide; i++) {
            pch[3*i + 2] = (pnImage[i+j*nWide] >>  0) & 0xFF;   // blue
            pch[3*i + 1] = (pnImage[i+j*nWide] >>  8) & 0xFF;   // green
            pch[3*i + 0] = (pnImage[i+j*nWide] >> 16) & 0xFF;   // red
        }
        WriteFile(hFile, pch, nLine, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
    delete [] pch;
    return 1;
}

unsigned
ColorPixel(int nRed, int nGreen, int nBlue)
{
    if (nRed < 0) nRed = 0;
    if (nRed > 255) nRed = 255;
    if (nGreen < 0) nGreen = 0;
    if (nGreen > 255) nGreen = 255;
    if (nBlue < 0) nBlue = 0;
    if (nBlue > 255) nBlue = 255;
    return unsigned(nRed + (nGreen << 8) + (nBlue << 16));
}

//
//  Write out test images for debugging and inspection of data.
//      Column 1:       any-difference image
//      Columns 2 to 5: raw image transmissions
//      Column 6:       master version
//
static inline char
Pixel(unsigned short n)
{
    if (n & BADPIXEL) {
        if ((n >> 7) == 0)
            return 0;
        else
            return 1;
    } else
        return n >> 8;
}

void
DebugImage(unsigned short rgnDiff[MAXLINES][252],
           unsigned short rgnMaster[MAXLINES][252], char *szTestImage)
{
    HANDLE hFile;
    char rgch[252*7];
    DWORD nBytesWritten;
    int i, j, k;
    //
    //  Write out test images for debugging and inspection of data.
    //      Column 1:       any-difference image
    //      Columns 2 to 5: raw image transmissions
    //      Column 6:       master version
    //
    hFile = CreateFile(szTestImage, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                                    FILE_ATTRIBUTE_NORMAL, NULL);
    WriteHeaderBMP(hFile, 252*7, MAXLINES, 1);
    for (i = 0; i < 252*7; i++)
        rgch[i] = 0;
    for (j = 0; j < MAXLINES; j++) {
        for (i = 0; i < 252; i++)
            rgch[i + 252*4] = Pixel(rgnDiff[j][i]);
        for (k = 0; k < 4; k++) {
            for (i = 0; i < 252; i++)
                rgch[i + 252*k] = Pixel(g_rgnRaw[k][j][i]);
        }
        for (i = 0; i < 252; i++)
            rgch[i + 252*5] = Pixel(rgnMaster[j][i]);
        for (i = 0; i < 252; i++)
            rgch[i + 252*6] = Pixel(g_rgnBrown[j][i]);
        WriteFile(hFile, rgch, 252*7, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
}

void
DebugFourImages(unsigned short rgn1[MAXLINES][252], 
                 unsigned short rgn2[MAXLINES][252],
                 unsigned short rgn3[MAXLINES][252],
                 unsigned short rgn4[MAXLINES][252],
                 char *szTestImage)
{
    HANDLE hFile;
    char rgch[252*4];
    DWORD nBytesWritten;
    int i, j;

    hFile = CreateFile(szTestImage, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                                    FILE_ATTRIBUTE_NORMAL, NULL);
    WriteHeaderBMP(hFile, 252*4, MAXLINES, 1);
    for (i = 0; i < 252*4; i++)
        rgch[i] = 0;
    for (j = 0; j < MAXLINES; j++) {
        for (i = 0; i < 252; i++)
            rgch[i + 252*0] = Pixel(rgn1[j][i]);
        for (i = 0; i < 252; i++)
            rgch[i + 252*1] = Pixel(rgn2[j][i]);
        for (i = 0; i < 252; i++)
            rgch[i + 252*2] = Pixel(rgn3[j][i]);
        for (i = 0; i < 252; i++)
            rgch[i + 252*3] = Pixel(rgn4[j][i]);
        WriteFile(hFile, rgch, 252*4, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
}

int
WriteHeaderWAVE(HANDLE hFile, int nSamples, int nSampleRate)
{
    DWORD nBytesWritten, n4;
    unsigned short n2;

    //
    //  RIFF header
    //
    WriteFile(hFile, "RIFF", 4, &nBytesWritten, NULL);
    n4 = 36 + 2*nSamples;
    WriteFile(hFile, &n4, 4, &nBytesWritten, NULL);
    WriteFile(hFile, "WAVE", 4, &nBytesWritten, NULL);
    //
    //  WAVE header
    //
    WriteFile(hFile, "fmt ", 4, &nBytesWritten, NULL);
    n4 = 16;
    WriteFile(hFile, &n4, 4, &nBytesWritten, NULL);   // PCM
    n2 = 1;
    WriteFile(hFile, &n2, 2, &nBytesWritten, NULL);   // uncompressed linear quantize
    n2 = 1;
    WriteFile(hFile, &n2, 2, &nBytesWritten, NULL);   // one channel
    n4 = nSampleRate;
    WriteFile(hFile, &n4, 4, &nBytesWritten, NULL);   // sample rate (8000, 44100,...)
    n4 = nSampleRate*1*2;
    WriteFile(hFile, &n4, 4, &nBytesWritten, NULL);   // byte rate
    n2 = 2;
    WriteFile(hFile, &n2, 2, &nBytesWritten, NULL);   // block alignment
    n2 = 16;
    WriteFile(hFile, &n2, 2, &nBytesWritten, NULL);   // bits per sample
    //
    //  Start of data subchunk, just write 16-bit samples after all this
    //
    WriteFile(hFile, "data", 4, &nBytesWritten, NULL);
    n4 = nSamples * 2;
    WriteFile(hFile, &n4, 4, &nBytesWritten, NULL);   // data chunk size
    return 1;
}










