#include "stdafx.h"
//
//  Process 8-bit image data from Brown University
//  D. P. Mitchell  05/15/2003.
//
#include "Venus.h"

//
//  The Brown images were separate panoramas.  I concatinated them into
//  four single strips, with calibration data restored to its correct location,
//  so they would resemble the original telemetry format.
//
static unsigned char s_rgnBuffer[MAXLINES][252];

int
LoadBrownImage(unsigned short rgnRaw[MAXLINES][252], char *szFile, int jFirst)
{
    HANDLE hFile;
    BOOL bResult;
    DWORD nBytesRead, nBytesToRead;
    int i, j, nLines, jShift;

    for (j = 0; j < MAXLINES; j++) {
        for (i = 0; i < 252; i++) {
            rgnRaw[j][i] = BADPIXEL | BURSTPIXEL;
            g_rgnTest1[j][i] = BADPIXEL;
        }
    }
    hFile = CreateFile(szFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                    FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    Assert(hFile != INVALID_HANDLE_VALUE);
    nBytesToRead = 1*252*(MAXLINES - jFirst);
    bResult = ReadFile(hFile, s_rgnBuffer[jFirst], nBytesToRead, &nBytesRead, NULL);
    Assert(bResult && nBytesRead != 0 && nBytesRead%252 == 0);
    nLines = nBytesRead/252;
    CloseHandle(hFile);
    for (j = 0; j < nLines; j++)
        for (i = 0; i < 252; i++)
            rgnRaw[jFirst + j][i] = s_rgnBuffer[jFirst + j][i] << 8;
    if (jFirst > 4000)
        jShift = 5000;
    else
        jShift = 0;
    printf("     %s %4d:", szFile, nLines);
    for (i = 0; i < 14; i++)
        printf(" %03o", rgnRaw[jShift + 800+11*i][50+10*i] >> 7);
    printf("\n");
    return nLines;
};
//
//  Some panoramas are rotated by 1, my photoshop mistake
//
void
FixBrown()
{
    int j, jFirst, jLimit, nZero1, nZero2;

    printf("    Fix up Brown University images\n");
    jFirst = 0;
    while (IsBad(g_rgnBrown[jFirst][50]))
        jFirst++;
    for (;;) {
        //
        //  The individual Brown University images were processed in photoshop to
        //  put them and their calibration sections into telemetry position.  There
        //  is a two-column gap of zero between each image.
        //
        while (g_rgnBrown[jFirst][100] == 0 && g_rgnBrown[jFirst][50] == 0)
            jFirst++;
        if (IsBad(g_rgnBrown[jFirst][100]))
            return;
        jLimit = jFirst + 1;
        while (g_rgnBrown[jLimit][100] != 0 || g_rgnBrown[jLimit][50] != 0)
            jLimit++;
        nZero1 = nZero2 = 0;
        for (j = jFirst; j < jFirst + 44; j++) {
            nZero1 += ((g_rgnBrown[j][0]>>7) == 0);
            nZero2 += ((g_rgnBrown[j][1]>>7) == 0);
        }
        Assert(nZero1 > 22 || nZero2 > 22);
        if (nZero1 < 22 && nZero2 > 22)
            RotateColumns(g_rgnBrown, jFirst, jLimit, 1, ROT_ROTATE);
        jFirst = jLimit;
    }
}
//
//  Find Brown pixel corresponding to real telemetry.  This is complicated by
//  the accidental mangling of the data when it was originaly extracted from
//  the Russian tapes.  About 2 percent of the pixels are missing.
//
#define BROWN_MATCH(B, P) (IsGood(P) && abs((B) - (P)) <= 0x0100)
#define NB      3
#define NBLIM   15
#define BROWN_VERBOSE   1
#define BROWN_REPAIR    2

unsigned short
BrownPixel(int jRaw, int iRaw, int nMode)
{
    int i, j, kk, k, kBrown, jBrown, jMax, nQuad1, nQuad2, nQuad3, nQuad4, nTmp;
    unsigned nMaxScore;
    int nPixel, rgnNHood[2*NB+1][2*NB+1], rgnMatch[2*NB+1][2*NB+1];
    union {
        unsigned nScore;
        char     rgnQuad[4];
    } u;

    //
    //  Find the transmission that matches the Brown University version, and save a
    //  neighborhood region of it around (jRaw, iRaw).
    //
    kBrown = -1;
    for (kk = 0; kk < 2; kk++) {
        if ((k = g_rgkBrown[kk]) == -1)
            continue;
        if (IsGood(g_rgnRaw[k][jRaw-1][iRaw]) || IsGood(g_rgnRaw[k][jRaw+1][iRaw]))
            kBrown = k;
    }
    if (kBrown == -1)
        return BADPIXEL;
    if (jRaw < NB || jRaw >= MAXLINES-NB)
        return BADPIXEL;
    for (j = -NB; j <= NB; j++) {
        for (i = -NB; i <= NB; i++) {
            rgnNHood[j+NB][i+NB] = 0;
            if (i+iRaw < 0 || i+iRaw >= 252)
                continue;
            rgnNHood[j + NB][i + NB] = g_rgnRaw[kBrown][jRaw + j][iRaw + i];
        }
    }
    //
    //  Brown data does not sync up horizontally with real data because of
    //  lost lines, but it is vertically aligned.
    //
    nMaxScore = 0;
    jMax = -1;
    for (jBrown = jRaw - 190; jBrown < jRaw + 5; jBrown++) {
        if (!BROWN_MATCH(g_rgnBrown[jBrown-1][iRaw], rgnNHood[NB-1][NB]) &&
            !BROWN_MATCH(g_rgnBrown[jBrown+1][iRaw], rgnNHood[NB+1][NB]))
            continue;
        for (j = -NB; j <= NB; j++) {
            for (i = -NB; i <= NB; i++) {
                rgnMatch[j+NB][i+NB] = 0;
                if (i+iRaw < 0 || i+iRaw >= 252)
                    continue;
                rgnMatch[j+NB][i+NB] = BROWN_MATCH(g_rgnBrown[jBrown + j][iRaw + i], rgnNHood[j + NB][i + NB]);
            }
        }
        //
        //  Subcolumns may be missing to the left or right, or the
        //  rows above or below can be horizontally shifted.
        //
        nQuad1 = nQuad2 = nQuad3 = nQuad4 = 0;
        for (j = 0; j <= NB; j++) {
            for (i = 0; i <= NB; i++) {
                nQuad1 += rgnMatch[j][i];
                nQuad2 += rgnMatch[j + NB][i];
                nQuad3 += rgnMatch[j][i + NB];
                nQuad4 += rgnMatch[j + NB][i + NB];
            }
        }
        //
        //  Find best match, based on lexical order of sorted quad matches
        //
        u.rgnQuad[0] = nQuad1;
        u.rgnQuad[1] = nQuad2;
        u.rgnQuad[2] = nQuad3;
        u.rgnQuad[3] = nQuad4;
        for (j = 1; j < 4; j++) {
            nTmp = u.rgnQuad[j];
            for (i = j; i > 0 && nTmp <u.rgnQuad[i - 1]; i = i - 1)
                u.rgnQuad[i] = u.rgnQuad[i - 1];
            u.rgnQuad[i] = nTmp;
        }
        if (u.nScore > nMaxScore) {
            nMaxScore = u.nScore;
            jMax = jBrown;
        }
    }
    if (nMode == BROWN_VERBOSE)
        printf("%4d -> %4d %08x %3o - %3o = %3o\n", jRaw, jMax, nMaxScore,
            rgnNHood[NB][NB] >> 7,
            g_rgnBrown[jMax][iRaw] >> 7,
            abs(g_rgnBrown[jMax][iRaw] - rgnNHood[NB][NB]));
    if (nMaxScore < 0x0E0E0000)
        return BADPIXEL;
    else {
        nPixel = g_rgnBrown[jMax][iRaw];
        if (nMode == BROWN_REPAIR && !BROWN_MATCH(nPixel, g_rgnRaw[kBrown][jRaw][iRaw])
            && ((nQuad1>=15 && nQuad2>=15)||(nQuad3>=15 && nQuad4>=15))) {
            g_rgnRaw[kBrown][jRaw][iRaw] = nPixel;
            g_rgnTest2[jRaw][iRaw] = nPixel;
        }
        return nPixel;
    }
}

int
BrownReplace(int k, int j, int i)
{
    unsigned short nPixel, nAbove;

    if (g_rgnFlags[k] & FLAGS_BROWN) {
        if (IsBad(nPixel = BrownPixel(j, i)))
            return 0;
        if (i > 0) {
            nAbove = BrownPixel(j, i-1);
            if (!BROWN_MATCH(nAbove, g_rgnRaw[k][j][i-1]))
                return 0;
        }
        g_rgnRaw[k][j][i] = nPixel;
        return 1;
    }
    return 0;
}
//
//  Repair a strip of raw pixels with 8-bit data.  A last resort.  bBadOnly should be
//  used when good pixels from non-brown-transmission are merged into the image.
//
int
BrownRepair(int k, int j, int iFirst, int iLimit, int bBadOnly)
{
    int i, nPixel, nFixed;

    if ((g_rgnFlags[k] & FLAGS_BROWN) == 0)
        return 0;
    nFixed = 0;
    for (i = iLimit-1; i >= iFirst; --i) {
        if (bBadOnly && IsGood(g_rgnRaw[k][j][i]))
            continue;
        nPixel = BrownPixel(j, i);
        if (IsBad(nPixel))
            continue;
        if (BROWN_MATCH(nPixel, g_rgnRaw[k][j][i]))
            continue;                               // don't replace a good 9-bit pixel
        g_rgnRaw[k][j][i] = nPixel;
        nFixed++;
    }
    for (i = iFirst; i < iLimit; i++) {
        if (bBadOnly && IsGood(g_rgnRaw[k][j][i]))
            continue;
        nPixel = BrownPixel(j, i);
        if (IsBad(nPixel))
            continue;
        if (BROWN_MATCH(nPixel, g_rgnRaw[k][j][i]))
            continue;                               // don't replace a good 9-bit pixel
        g_rgnRaw[k][j][i] = nPixel;
        nFixed++;
    }

    printf("Brown Repair %d of %d pixels\n", nFixed, iLimit - iFirst);
    return nFixed;
}

void
VisualizeBrown(int jFirst, int jLimit, int nMode)
{
    int i, j, nPixel;

    BadPixels(g_rgnTest2, 0, MAXLINES);
    for (j = jFirst; j < jLimit; j++)
        for (i = 251; i >= 0; --i) {
            nPixel = BrownPixel(j, i, nMode);
            if (nMode != BROWN_REPAIR)
                g_rgnTest2[j][i] = nPixel;
        }
}

void
CopyBrown(unsigned short rgnData[MAXLINES][252], int jData, int jBrown,
          int iFirst, int iLimit)
{
    int i;

    for (i = iFirst; i < iLimit; i++)
        rgnData[jData][i] = g_rgnBrown[jBrown][i];
}