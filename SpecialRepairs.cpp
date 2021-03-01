#include "stdafx.h"
//
//  Repair various special regions of the Venera telemetry
//  D. P. Mitchell  05/15/2003.
//
#pragma intrinsic(log)
#include "Venus.h"
#include "ImageProcessing.h"
#include "ImageFile.h"

//
//  Visualization aids
//
static int
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

static int
HammingDistance(int n1, int n2)
{
    int i, n, nHam;

    n = n1 ^ n2;
    nHam = 0;
    for (i = 0; i < 16; i++)
        if (n & (1 << i))
            nHam++;
    return nHam;
}

static char *
Binary(int n, char *sz = 0)
{
    static char rgch[16];
    int i;

    if (sz == 0)
        sz = rgch;
    sz[9] = 0;
    for (i = 8; i >= 0; --i)
        sz[i] = '0' + ((n >> i) & 1);   // in order of transmission (LSB first)
    return sz;
}
//
//  Here is typical off-by-one bit sync error found in Venera 14-I:
//
//  D  29  45  39  33 101110000_ 101101000_ 111001000_ 100001000_
//  C 270 278 275 272 0111000011 0110100011 1100100011 0000100011
//  D  39  27  33  25 111001000_ 110110000_ 100001000_ 100110000_
//  C 275 269 272 268 1100100011 1011000011 0000100011 0011000010
//
//  A few lines are off-by-two, which proves the 10-bit-block theory:
//
//  D 316 422 476   0 001111001_ 011001011_ 001110111_ 000000000_
//  C 207 233 255 256 1111001101 1001011100 1111111101 0000000010
//  D  33  55  35  35 100001000_ 111011000_ 110001000_ 110001000_
//  C 264 269 264 264 0001000011 1011000011 0001000011 0001000011
//
//  Some regions are complemented, as here in 13-I, line 2433:
//
//  D 188 199 199 200 001111010_ 111000110_ 111000110_ 000100110_
//  C 323 312 312 311 1100001011 0001110011 0001110011 1110110011
//  D 198 194 312 303 011000110_ 010000110_ 000111001_ 111101001_
//  C 313 317 312 303 1001110010 1011110011 0001110011 1111010011
//
static void
Display(unsigned short rgnData[MAXLINES][252], unsigned short rgnCheck[MAXLINES][252],
        int j, int iFirst = 0, int iLimit = 252)
{
    int i, i4, nData;

    for (i = iFirst; i < iLimit; i += 4) {
        printf("D %3d ", i%252);
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnData[j][i4] >> 7;
            printf("%3d ", nData);
        }
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnData[j][i4] >> 7;
            printf("%s_ ", Binary(nData));
        }
        printf("\nC     ");
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnCheck[j][i4] >> 7;
            printf("%3d ", nData);
        }
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnCheck[j][i4] >> 7;
            printf("%s%d ", Binary(nData), Parity(nData));
        }
        printf("\n");
    }
}
//
//  Trying to get a handle on the mysterious "secondary damage".  The older Brown Univ.
//  data shows bitsync errors (but since it is not 9-bit data, we cannot repair it by
//  rebuilding the 10th parity bit).  The Russian 9-bit data shows damage in the same
//  position, but changed and apparently not simply a different or additional bit shift.
//
//  Here is an example from line 1439 of 14-II.
//  
//  D 107 257 311 305 307 100000001_ 111011001_ 100011001_ 110011001_
//  B      14 444 396 410 011100000_ 001111011_ 001100011_ 010110011_
//  C     325 320 315 312 1010001011 0000001011 1101110011 0001110011
//  D 111 308 312 316 257 001011001_ 000111001_ 001111001_ 100000001_
//  B     418 450 482  12 010001011_ 010000111_ 010001111_ 001100000_
//  C     313 314 316 318 1001110010 0101110010 0011110010 0111110011
//
static void
Display3(unsigned short rgnData[MAXLINES][252], unsigned short rgnCheck[MAXLINES][252],
         unsigned short rgnBrown[MAXLINES][252],
        int j, int jBrown, int iFirst = 0, int iLimit = 252)
{
    int i, i4, nData;

    for (i = iFirst; i < iLimit; i += 4) {
        printf("D %3d ", i%252);
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnData[j][i4] >> 7;
            printf("%3d ", nData);
        }
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnData[j][i4] >> 7;
            printf("%s_ ", Binary(nData));
        }
        printf("\nB     ", i%252);
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnBrown[jBrown][i4] >> 7;
            printf("%3d ", nData);
        }
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnBrown[jBrown][i4] >> 7;
            printf("%s_ ", Binary(nData));
        }
        printf("\nC     ");
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnCheck[j][i4] >> 7;
            printf("%3d ", nData);
        }
        for (i4 = i; i4 < i+4; i4++) {
            nData = rgnCheck[j][i4] >> 7;
            printf("%s%d ", Binary(nData), Parity(nData));
        }
        printf("\n");
    }
}

static void
BitSyncImage(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit)
{
    unsigned short *pn;
    int i, j, nShift, jWide;

    jWide = jLimit - jFirst;
    pn = new unsigned short[252*10*jWide];
    for (nShift = -4; nShift < 6; nShift++) {
        for (j = jFirst; j < jLimit; j++) {
            for (i = 0; i < 252; i++)
                g_rgnTest2[j][i] = rgnData[j][i];
        }
        BitShiftPixels(g_rgnTest2, nShift, jFirst, jLimit, 0, 252);
        for (j = jFirst; j < jLimit; j++) {
            for (i = 0; i < 252; i++)
                pn[j - jFirst + jWide*(i + 252*(nShift+4))] = g_rgnTest2[j][i];
        }
    }
    WriteShortImage(pn, "BitSync.bmp", 252*10, jWide);
    delete [] pn;
}

static void
AutoBitShift(unsigned short rgnData[MAXLINES][252], int jData,
             unsigned short rgnCheck[MAXLINES][252], int jCheck,
             int iFirst=0, int iLimit=252, int nMode = 0)
{
    double fCorr;
    int i, nShift, nShift2, jDraw, nFirst2, nLimit2;
    unsigned short rgnSave[252];

    printf("Bit-Sync Correlations at %4d:\n", jData);
    for (i = 0; i < 252; i++)
        rgnSave[i] = rgnData[jData][i];
    //
    //  nShift2 is to check for regions which have undergone two shifts
    //
    jDraw = jData - 50;
    if (0 && nMode == 1) {
        nFirst2 = 0;
        nLimit2 = 1;
    } else {
        nFirst2 = -4;
        nLimit2 = 6;
    }
    for (nShift2 = nFirst2; nShift2 < nLimit2; nShift2++) {
        for (nShift = -4; nShift < 6; nShift++) {
            BitShiftPixels(rgnData, nShift2, jData, jData+1, iFirst, iLimit);
            BitShiftPixels(rgnData, nShift, jData, jData+1, iFirst, iLimit);
            fCorr = Correlation(&rgnData[jData][iFirst], &rgnCheck[jCheck][iFirst], iLimit-iFirst);
            if (nMode || fCorr > 0.33 || fCorr < -0.33) {
                printf("%3d %3d:", nShift2, nShift);
                printf(" %6.3f\n", fCorr);
                for (i = 0; i < 252; i++)
                    rgnData[jDraw][i] = rgnData[jData][i];
                --jDraw;
            }
            for (i = 0; i < 252; i++)
                rgnData[jData][i] = rgnSave[i];
        }
    }
}
//
//  Some light-gray noise in 14-II seems to have low-order bits set to zero
//  or one:
//
//  Bad:    256 000000001 287 111110001 256 000000001 258 010000001
//  Good:   283 110110001 278 011010001 285 101110001 286 011110001
//  Bad:    284 001110001 264 000100001 256 000000001 287 111110001
//  Good:   276 001010001 280 000110001 275 110010001 284 001110001
//  Bad:    287 111110001 256 000000001 258 010000001 256 000000001
//  Good:   259 110000001 266 010100001 271 111100001 280 000110001
//  Bad:    288 000001001 288 000001001 108 001101100 224 000001110
//  Good:   292 001001001 314 010111001 124 001111100 255 111111110
//
void
ExamineErrors(unsigned short rgnData[MAXLINES][252], unsigned short rgnCheck[MAXLINES][252],
              int jFirst, int jLimit)
{
    int i, j, n, rgnBad[4], rgnGood[4];
    int rgnHammingHist[16];

    n = 0;
    for (i = 0; i < 16; i++)
        rgnHammingHist[i] = 0;
    for (j = jFirst; j < jLimit; j++) {
        for (i = 0; i < 252; i++) {
            if (IsGood(rgnCheck[j][i])) {
                if (rgnData[j][i] != rgnCheck[j][i]) {
                    if (n >= 4) {
                        while (--n >= 0)
                            printf("%3d %s ", rgnBad[n], Binary(rgnBad[n]));
                        printf("\n");
                        n = 4;
                        while (--n >= 0)
                            printf("%3d %s ", rgnGood[n], Binary(rgnGood[n]));
                        printf("\n\n");
                        n = 0;
                    }
                    rgnBad[n] = rgnData[j][i] >> 7;
                    rgnGood[n++] = rgnCheck[j][i] >> 7;
                    rgnHammingHist[HammingDistance(rgnData[j][i] >> 7, rgnCheck[j][i] >> 7)]++;
                } else
                    rgnHammingHist[0]++;
            }
        }
    }
    printf("Hamming Distances of Errors:\n");
    for (i = 0; i < 16; i++)
        printf("%2d %5d\n", i, rgnHammingHist[i]);
}

//
//  Fix up a few bitstream errors in Venera 13 Camera I
//
static void
Complement(unsigned short rgnData[MAXLINES][252], int bBadFirst, int bBadLimit,
           int j, int iFirst, int iLimit)
{
    int i;

    Assert(iLimit <= 252);
    Assert(iFirst < 252);
    for (i = iFirst; i < iLimit; i++)
        rgnData[j][i] = (~rgnData[j][i]) & NINEBITS;
    if (bBadFirst)
        SetBad(rgnData[j][iFirst]);
    if (bBadLimit)
        SetBad(rgnData[j][iLimit-1]);
}

#define NCOMPSEG 42

static void
AutoComplement(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit)
{
    int i, iSeg, j, nDark;

    for (j = jFirst; j < jLimit; j++) {
        for (iSeg = 0; iSeg < 252; iSeg += NCOMPSEG) {
            nDark = 0;
            for (i = iSeg; i < iSeg+NCOMPSEG; i++)
                if ((rgnData[j][i] >> 7) < 256)
                    nDark++;
            if (nDark < NCOMPSEG/2)
                for (i = iSeg; i < iSeg+NCOMPSEG; i++)
                    rgnData[j][i] = (~rgnData[j][i]) & NINEBITS;
        }
    }
}

#define NQUANT (1 << 3)

static double
NegEntropy(unsigned short rgnData[], int n)
{
    int i, nQuantum, rgnHist[NQUANT];
    double fProb, fEntropy;

    for (i = 0; i < NQUANT; i++)
        rgnHist[i] = 0;
    for (i = 0; i < n; i++) {
        nQuantum = rgnData[i] >> (7 + 9 - 3);
        rgnHist[nQuantum]++;
    }
    fEntropy = 0.0;
    for (i = 0; i < NQUANT; i++) {
        fProb = double(rgnHist[i])/double(n);
        if (fProb)
            fEntropy += double(rgnHist[i])*log(fProb);
    }
    return fEntropy;
}

static void
MaxEntropy(unsigned short rgnData[], int n)
{
    double fEntropy, fMaxEntropy;
    int i, nShift, nMaxShift, bComplement;
    unsigned short rgn[252];

    Assert(n < 252);
    fMaxEntropy = NegEntropy(rgnData, n);
    bComplement = 0;
    nMaxShift = 0;
    for (nShift = -4; nShift < 6; nShift++) {
        for (i = 0; i < n; i++)
            rgn[i] = rgnData[i];
        BitShift(rgn, nShift, n);
        fEntropy = NegEntropy(rgn, n);
        if (fEntropy > fMaxEntropy) {
            fMaxEntropy = fEntropy;
            nMaxShift = nShift;
        }
    }
    BitShift(rgnData, nMaxShift, n);
}

#define NSEGENT NCOMPSEG

static void
MaxEntropyRepair(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit)
{
    int i, j;

    for (j = jFirst; j < jLimit; j++)
        for (i = 0; i < 252; i += NSEGENT)
            MaxEntropy(&rgnData[j][i], NSEGENT);
}
//
//  The strange transform:
//  Some regions have an out-of-sync bitstream, but are further damaged.  Often
//  the Brown University version does not have the "further damage", but it only
//  contains the high order 7 bits of correct information.  May be able to combine
//  Brown and Russian data to assemble the needed 9 bits.
//
//  D 268 267 257 272 001100001_ 110100001_ 100000001_ 000010001_
//  C 280 278 258 289 0001100010 0110100011 0100000011 1000010010
//  D 279 276 268 269 111010001_ 001010001_ 001100001_ 101100001_
//  C 302 296 280 282 0111010010 0001010010 0001100010 0101100011
//
static void
StrangeTransform(unsigned short rgnData[MAXLINES][252], int j, int iFirst, int iLimit)
{
    int i;

    BadPixels(rgnData, j, j+1, iFirst, iLimit); //REVIEW: not adequately solved yet
    return;
    //
    //  This is just an approximate fix, the low order bit can be wrong
    //
    BitShiftPixels(rgnData, -1, j, j+1, iFirst, iLimit);
    for (i = iFirst; i < iLimit; i++)
        if (rgnData[j][i] < 200 || rgnData[j][i] > 340)
            rgnData[j][i] = rgnData[j][i] | 0x8000;
}
//
//  13-I is the cleanest and most extensive telemetry.  Some data is masked out by
//  hand where it seemed to be confusing the MasterVersion routine.  A few parts can
//  be repaired by complementing data ("half-bit out of sync").
//
void
SpecialV13C1()
{
    printf("  Fix bitstream errors in 13-I\n");
    //AutoBitShift(g_rgnRaw[2], 5374, g_rgnRaw[2], 5373, 98, 120);
    //Display(g_rgnRaw[0], g_rgnRaw[1], 2433);
    BadPixels(g_rgnRaw[1], 926, 927, 68, 101);
    BadPixels(g_rgnRaw[0], 1146, 1159, 24, 163);
    BadPixels(g_rgnRaw[1], 1195, 1196, 80, 129);
    BadPixels(g_rgnRaw[1], 1207, 1253);             // subtle errors in red panorama
    BadPixels(g_rgnRaw[0], 2272, 2273, 47, 65);
    BadPixels(g_rgnRaw[1], 2429, 2430, 241, 252);
    BadPixels(g_rgnRaw[0], 2690, 2692, 107, 70);
    BadPixels(g_rgnRaw[1], 2630, 2631, 234, 249);
    BadPixels(g_rgnRaw[0], 3243, 3244, 52, 72);     // 3rd clear, noisy area in k=1
    BadPixels(g_rgnRaw[0], 3376, 3377, 169, 181);
    BadPixels(g_rgnRaw[0], 3425, 3426, 136, 154);
    BadPixels(g_rgnRaw[1], 3428, 3486, 412-252, 252);   // many faint errors
    BadPixels(g_rgnRaw[1], 3426, 3427, 141, 142);
    Complement(g_rgnRaw[0], 1, 0, 3426, 148, 169);

    Complement(g_rgnRaw[2], 0, 0, 5374, 36, 40);
    Complement(g_rgnRaw[2], 1, 0, 5374, 43, 47);
    Complement(g_rgnRaw[2], 0, 0, 5374, 138, 140);
    Complement(g_rgnRaw[2], 0, 1, 5374, 142, 148);
    Complement(g_rgnRaw[2], 0, 1, 5374, 151, 154);
    Complement(g_rgnRaw[2], 0, 1, 5374, 155, 170);
    Complement(g_rgnRaw[2], 1, 1, 5375, 56, 84);
    Complement(g_rgnRaw[2], 0, 1, 5376, 16, 26);
    BadPixels(g_rgnRaw[2], 6071, 6072, 1688%252, 1734%252); // subtle but destroyed

    Complement(g_rgnRaw[2], 0, 0, 6481, 128, 180);
    Complement(g_rgnRaw[2], 0, 0, 6481, 189, 252);
    CopyBrown(g_rgnRaw[2], 6481, 6316, 1309%252, 1331%252);
    Complement(g_rgnRaw[2], 1, 1, 6481, 1309%252, 1331%252);

    Complement(g_rgnRaw[0], 1, 0, 2430, 197, 207);
    Complement(g_rgnRaw[0], 1, 1, 2430, 225, 249);
    Complement(g_rgnRaw[0], 1, 1, 2431, 13, 24);
    Complement(g_rgnRaw[0], 1, 1, 2431, 32, 42);
    Complement(g_rgnRaw[0], 0, 1, 2431, 53, 67);
    Complement(g_rgnRaw[0], 1, 1, 2431, 141, 152);
    Complement(g_rgnRaw[0], 0, 1, 2431, 157, 166);
    Complement(g_rgnRaw[0], 0, 1, 2431, 192, 209);
    Complement(g_rgnRaw[0], 0, 0, 2431, 213, 222);
    Complement(g_rgnRaw[0], 0, 0, 2431, 481-252, 491-252);
    Complement(g_rgnRaw[0], 0, 0, 2432, 9, 15);
    Complement(g_rgnRaw[0], 0, 0, 2432, 20, 25);
    Complement(g_rgnRaw[0], 0, 0, 2432, 28, 33);
    Complement(g_rgnRaw[0], 0, 1, 2432, 37, 44);
    BadPixels(g_rgnRaw[0], 2432, 2442, 43, 43);     // Lots of complemented segments
    BadPixels(g_rgnRaw[1], 2431, 2432, 191);
}
//
//  13-II was full of zero errors.  The last red and green panoramas end in very
//  bad noise.  A little of it can be cleaned up by bit-stream resync and time-base
//  correction.  There must exist other forms of bit-stream resyncs, perhaps involving
//  complementing as well as shifting (e.g., "2.5 bits out of sync").
//
void
SpecialV13C2()
{
    printf("  Fix bitstream errors in 13-II\n");
    //BitSyncImage(g_rgnRaw[1], 5000, 6620);
    //
    //  Clean up the noisy stuff a little in the last two panoramas
    //
    BitShiftPixels(g_rgnRaw[1], -4, 5818, 5823, 0, 252);
    BitShiftPixels(g_rgnRaw[1], -4, 6314, 6316, 1033%252, 1065%252);
    BitShiftPixels(g_rgnRaw[1],  3, 6199, 6202, 0, 1151%252);
    BitShiftPixels(g_rgnRaw[1], -4, 6140, 6141, 1097%252, 252);
    BitShiftPixels(g_rgnRaw[1],  2, 6002, 6006, 0, 252);
    BitShiftPixels(g_rgnRaw[1],  4, 6272, 6273, 0, 252);
    BitShiftPixels(g_rgnRaw[1],  4, 6253, 6254, 0, 252);
    BitShiftPixels(g_rgnRaw[1],  4, 6260, 6261, 0, 252);
    BitShiftPixels(g_rgnRaw[1],  4, 6262, 6263, 2152%252, 252);
    BitShiftPixels(g_rgnRaw[1],  4, 6259, 6260, 2194%252, 252);
    BitShiftPixels(g_rgnRaw[1],  4, 6383, 6389, 0, 252);
    BitShiftPixels(g_rgnRaw[1],  2, 6378, 6381, 0, 252);
    BitShiftPixels(g_rgnRaw[1], -2, 6406, 6413, 0, 252);
    TimeBaseCorrection(g_rgnRaw[1], 6180, 6550, 6179, TBC_NOSHIFTING);
    CopyBrown(g_rgnRaw[1], 6315, 6188, 1454%252, 252);
    ShiftSubColumn(g_rgnRaw[1], 6316, 1454%252, 252);
    BadPixels(g_rgnRaw[1], 6515, MAXLINES);           // Noise near end
    BadPixels(g_rgnRaw[1], 6438, 6452, 0, 1433%252);
    BadPixels(g_rgnRaw[1], 6391, 6398);
    BadPixels(g_rgnRaw[1], 5988, 5997);
    BadPixels(g_rgnRaw[1], 5916, 5979);
    BadPixels(g_rgnRaw[1], 5893, 5904);
    BadPixels(g_rgnRaw[1], 5860, 5874, 0, 1310%252);
    BadPixels(g_rgnRaw[1], 5829, 5853);
}
//
//  14-I has a patch of bit-stream sync errors on the left horizion of the 3rd clear
//  panorama.  It ends with a noisy red (and green?) panorama with severe bit sync
//  and time sync errors.
//
void
SpecialV14C1()
{
    int k;

    printf("  Fix bitstream errors in 14-I\n");
    //BitSyncImage(g_rgnRaw[0], 3000, 4300);
    //AutoBitShift(g_rgnRaw[0], 3781, g_rgnRaw[2], 3781);
    //
    //  Repair the streaks on 3rd clear horizon
    //
    for (k = 0; k < 2; k++) {
        BadPixels(g_rgnRaw[k], 3397, 3398, 240, 252);

        BitShiftPixels(g_rgnRaw[k],   1, 3349, 3350, 0, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3350, 3351, 0, 88);
        BitShiftPixels(g_rgnRaw[k],   1, 3350, 3351, 102, 201);
        BitShiftPixels(g_rgnRaw[k],   1, 3350, 3351, 218, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3353, 3354, 1424%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3354, 3355, 1440%252, 225);
        BitShiftPixels(g_rgnRaw[k],   1, 3354, 3355, 239, 252);
        BitShiftPixels(g_rgnRaw[k],   3, 3354, 3355, 225, 239);

        BitShiftPixels(g_rgnRaw[k],   1, 3355, 3356, 1413%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3356, 3357, 1413%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3359, 3360, 1455%252, 1493%252);
        BitShiftPixels(g_rgnRaw[k],   1, 3359, 3360, 1506%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3360, 3361, 0, 1411%252);
        BitShiftPixels(g_rgnRaw[k],  -9, 3362, 3363, 1430%252, 252);
        BitShiftPixels(g_rgnRaw[k],  -9, 3363, 3364, 1397%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3364, 3365, 1396%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3364, 3365, 1056%252, 1064%252);

        BitShiftPixels(g_rgnRaw[k],   1, 3369, 3370, 1362%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3369, 3370, 1292%252, 1307%252);
        BitShiftPixels(g_rgnRaw[k],   1, 3370, 3371, 1362%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3374, 3375, 1396%252, 1486%252);
        BitShiftPixels(g_rgnRaw[k],   1, 3374, 3375, 1496%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3381, 3382, 1366%252, 1738%252);   //1491
        BitShiftPixels(g_rgnRaw[k],   1, 3385, 3386, 1435%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3386, 3387, 1337%252, 1405%252);

        BadPixels(g_rgnRaw[k], 3348, 3349, 236, 252);
        BitShiftPixels(g_rgnRaw[k],   2, 3375, 3376, 1648%252, 252);
        BitShiftPixels(g_rgnRaw[k],   1, 3375, 3376, 26, 43);
        BitShiftPixels(g_rgnRaw[k],   2, 3382, 3383, 1618%252, 252);
        BitShiftPixels(g_rgnRaw[k],   2, 3381, 3382, 1738%252, 252);
        BitShiftPixels(g_rgnRaw[k],   3, 3354, 3355, 1908%252, 1941%252);
        BitShiftPixels(g_rgnRaw[k],   3, 3360, 3361, 1928%252, 1993%252);   //2246
        BitShiftPixels(g_rgnRaw[k],   4, 3360, 3362, 2246%252, 252);
        BitShiftPixels(g_rgnRaw[k],   5, 3386, 3387, 2413%252, 252);
        Complement(g_rgnRaw[k], 1, 0, 3386, 2413%252, 252);
        BitShiftPixels(g_rgnRaw[k],   5, 3387, 3388, 2344%252, 252);
        BitShiftPixels(g_rgnRaw[k],   5, 3387, 3388, 2305%252, 2335%252);
        Complement(g_rgnRaw[k], 1, 0, 3387, 2344%252, 252);
        Complement(g_rgnRaw[k], 1, 0, 3387, 2305%252, 2335%252);
        BitShiftPixels(g_rgnRaw[k], -10, 3364, 3365, 1056%252, 1257%252);
        BitShiftPixels(g_rgnRaw[k],  10, 3382, 3383, 0, 80);

        BadPixels(g_rgnRaw[k], 3927, 3930, 27, 179);
        Complement(g_rgnRaw[k], 0, 0, 3929, 179, 235);
        BadPixels(g_rgnRaw[k], 3929, 3931, 235, 130);
        Complement(g_rgnRaw[k], 0, 0, 3930, 130, 252);
        //                      |||
        //  ^^^  BEFORE TBC,    vvv  AFTER TBC
        //  |||
        TimeBaseCorrection(g_rgnRaw[k], 3310, 3370, 3306, TBC_NOSHIFTING);  // danaged horizon
        TimeBaseCorrection(g_rgnRaw[k], 3559, 3927, 3559, TBC_NOSHIFTING);  // Last clear panorama
        TimeBaseCorrection(g_rgnRaw[k], 3931, 4155, 3559, TBC_NOSHIFTING);  // Last clear panorama
        TimeBaseCorrection(g_rgnRaw[k], 4155, 4255, 3559);  // L-W-EW event at 4200, shifts 1
        TimeBaseCorrection(g_rgnRaw[k], 4278, 4308, 1920);  // Last red panorama

        StrangeTransform(g_rgnRaw[k], 3382, 83, 105);
        StrangeTransform(g_rgnRaw[k], 3375, 62, 114);
        StrangeTransform(g_rgnRaw[k], 3354, 103, 127);
        StrangeTransform(g_rgnRaw[k], 3375, 62, 114);

        BadPixels(g_rgnRaw[k], 3349, 3350, 151, 175);
        CopyBrown(g_rgnRaw[k], 3349, 3293, 0, 97);      // Fill in after rotations
        CopyBrown(g_rgnRaw[k], 3320, 3264, 0, 167);
        CopyBrown(g_rgnRaw[k], 3320, 3265, 167, 252);
        CopyBrown(g_rgnRaw[k], 3781, 3717, 0, 252);
        CopyBrown(g_rgnRaw[k], 3817, 3753, 141, 252);
        CopyBrown(g_rgnRaw[k], 3880, 3815, 154, 252);
        CopyBrown(g_rgnRaw[k], 3909, 3843, 187, 207);   // Brown lost below 207
        CopyBrown(g_rgnRaw[k], 3925, 3859, 138, 211);
        CopyBrown(g_rgnRaw[k], 4090, 4021, 0, 56);
        CopyBrown(g_rgnRaw[k], 4090, 4021, 64, 75);
        CopyBrown(g_rgnRaw[k], 4090, 4021, 93, 110);
        CopyBrown(g_rgnRaw[k], 4147, 4077, 0, 158);     // Shear in Brown data
        CopyBrown(g_rgnRaw[k], 4147, 4078, 158, 252);
        CopyBrown(g_rgnRaw[k], 4152, 4083, 184, 252);
        CopyBrown(g_rgnRaw[k], 4166, 4097, 239, 252);
        CopyBrown(g_rgnRaw[k], 4175, 4105, 170, 252);
        CopyBrown(g_rgnRaw[k], 4178, 4108, 148, 252);
        CopyBrown(g_rgnRaw[k], 4194, 4124, 230, 252);
        CopyBrown(g_rgnRaw[k], 4197, 4127, 186, 252);
        CopyBrown(g_rgnRaw[k], 4200, 4129, 0, 105);
        CopyBrown(g_rgnRaw[k], 4202, 4131, 0, 126);
        CopyBrown(g_rgnRaw[k], 4203, 4133, 126, 252);
        CopyBrown(g_rgnRaw[k], 4214, 4144, 237, 252);
        CopyBrown(g_rgnRaw[k], 4231, 4160, 0, 86);

        //
        //  The noisy data at the end can be roughly re-synced, just to get
        //  a general idea of where the lander was in its program.
        //
        //  5578-5623 (-2)
        //  5623-5784 (-4)
        //  5784-5812 (+5)
        //  5812-5840 (+4)
        //  5855-5894 (+2)
        //
        //MaxEntropyRepair(g_rgnRaw[k], 4308, 5150);
        //AutoComplement(g_rgnRaw[k], 4308, 5150);
        //TimeBaseCorrection(g_rgnRaw[k], 4308, 5150, 1920, TBC_ROTATEONLY, 0.1, 2260);
    }
    return;
    BitShiftPixels(g_rgnRaw[0], -2, 4578, 4623, 0, 252);
    BitShiftPixels(g_rgnRaw[0], -4, 4623, 4784, 0, 252);
    BitShiftPixels(g_rgnRaw[0], +5, 4784, 4812, 0, 252);
    BitShiftPixels(g_rgnRaw[0], +4, 4812, 4840, 0, 252);
    BitShiftPixels(g_rgnRaw[0], +2, 4855, 4895, 0, 252);
    AutoComplement(g_rgnRaw[0], 4308, 4895);
    TimeBaseCorrection(g_rgnRaw[0], 4308, 4600, 1920, TBC_ROTATEONLY, 0.1);
    TimeBaseCorrection(g_rgnRaw[0], 4600, 4895, 2260, TBC_ROTATEONLY, 0.1);
    //Display(g_rgnRaw[0], g_rgnRaw[2], 3354, 100, 250);

    return;
    TimeBaseCorrection(g_rgnRaw[1], 3559, 4012);            // Last clear panorama
    TimeBaseCorrection(g_rgnRaw[0], 3559, 4012);
    TimeBaseCorrection(g_rgnRaw[1], 4024, 4255, 3559);      // Repairs L-W-EW event at 4200, shifts 1
    TimeBaseCorrection(g_rgnRaw[0], 4024, 4255, 3559);
    TimeBaseCorrection(g_rgnRaw[1], 4278, 4308, 1920);  // Last red panorama
    TimeBaseCorrection(g_rgnRaw[0], 4278, 4308, 1920);
    TimeBaseRecode(g_rgnRaw[0], 4308, 4615, 1920);      // Very noisy data (red panorama)
    TimeBaseRecode(g_rgnRaw[0], 4615, 5065, 2000);      // Is a green panorama buried in the noise?
    TimeBaseRecode(g_rgnRaw[1], 4308, 4745, 1920);
    TimeBaseCorrection(g_rgnRaw[0], 4308, 4615, 1920, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[0], 4615, 5065, 2000, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[1], 4308, 4745, 1920, TBC_NOSHIFTING);
    ShiftColumns(g_rgnRaw[0], 4357, 1);                 // Experiment with last red panorama
    ShiftColumns(g_rgnRaw[1], 4357, 1);
    ShiftColumns(g_rgnRaw[0], 4415, 1);
    ShiftColumns(g_rgnRaw[1], 4415, 1);

}
//
//  The left horizon in the last clear panorama is the only view of that region, so
//  cleaning it up is a top priority.  Each transmission seems to have different
//  bit-syncronization errors.
//
void
SpecialV14C2()
{
    int k;

    printf("  Fix bitstream errors in 14-II\n");
    Complement(g_rgnRaw[0], 1, 1, 1028, 232, 252);
    //
    //  Fix last panorama in transmission k=1
    //
    k = 1;
    BitShiftPixels(g_rgnRaw[k], -1, 4594, 4595, 938-756, 252);      // secondary damage
    BitShiftPixels(g_rgnRaw[k], -1, 4603, 4604, 953-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4612, 4613, 956-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4619, 4620, 839-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4620, 4622, 950-756, 946-756);
    BitShiftPixels(g_rgnRaw[k], -1, 4622, 4623, 811-756, 930-756);
    BitShiftPixels(g_rgnRaw[k], +3, 4634, 4636, 0, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4639, 4641, 772-756, 971-756);
    BitShiftPixels(g_rgnRaw[k], -1, 4648, 4649, 988-756, 252);      // secondary damage above
    BitShiftPixels(g_rgnRaw[k], -1, 4649, 4650, 988-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4650, 4651, 988-756, 252);
    BadPixels(g_rgnRaw[k], 4648, 4649, 0, 988-756);                 // secondary damage
    BadPixels(g_rgnRaw[k], 4650, 4651, 0, 988-756);
    BadPixels(g_rgnRaw[k], 4662, 4664, 485-252, 252);               // secondary damage
    BitShiftPixels(g_rgnRaw[k], -6, 4679, 4682, 1144-1008, 1038-1008);  //was +4
    BitShiftPixels(g_rgnRaw[k], -1, 4692, 4694, 915-756, 664-504);
    BitShiftPixels(g_rgnRaw[k], -2, 4693, 4697, 664-504, 535-504);
    BitShiftPixels(g_rgnRaw[k], -1, 4733, 4739, 970-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4739, 4740, 1156-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4768, 4770, 924-756, 873-756);
    BitShiftPixels(g_rgnRaw[k], -2, 4769, 4771, 873-756, 1054-1008);
    BadPixels(g_rgnRaw[k], 4786, 4787);                             // secondary damage
    BadPixels(g_rgnRaw[k], 4789, 4790);                             // secondary damage
    BitShiftPixels(g_rgnRaw[k], +4, 4791, 4795, 1133-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4809, 4810, 1102-1008, 638-504);
    BitShiftPixels(g_rgnRaw[k], -2, 4809, 4818, 638-504, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4858, 4869, 1026-1008, 994-756);
    BitShiftPixels(g_rgnRaw[k], +5, 4868, 4871, 994-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4897, 4902, 1177-1008, 252);
    BadPixels(g_rgnRaw[k], 4909, 4912, 1063-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4932, 4939, 1233-1008, 1053-1008);
    BitShiftPixels(g_rgnRaw[k], +3, 4975, 4977, 1107-1008, 1991%252);
    BitShiftPixels(g_rgnRaw[k], -1, 4979, 4980, 852-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4980, 4981, 839-756, 252);
    BitShiftPixels(g_rgnRaw[k], +4, 4788, 4789, 1222-1008, 252);    // secondary damage above
    BitShiftPixels(g_rgnRaw[k], +4, 4789, 4790, 1222-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4991, 4995, 932-756, 880-756);
    BitShiftPixels(g_rgnRaw[k], -1, 5022, 5025, 1092-1008, 1150-1008);
    BitShiftPixels(g_rgnRaw[k], -1, 5033, 5044, 0, 252);
    BitShiftPixels(g_rgnRaw[k], -3, 5044, 5048, 0, 457-252);
    BitShiftPixels(g_rgnRaw[k], +1, 5048, 5054, 457-252, 252);
    BitShiftPixels(g_rgnRaw[k], +4, 5056, 5059, 1120-1008, 2125%252);
    BitShiftPixels(g_rgnRaw[k], -2, 5059, 5061, 652-504, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 5064, 5079, 901-756, 873-756);
  //ShiftColumns(g_rgnRaw[1], 4813, 1);
    //BitSyncImage(g_rgnRaw[k], 4000, 5320);
    //TimeBaseCorrection(g_rgnRaw[k], 4624, 4909, 1700, TBC_NOSHIFTING);
    //TimeBaseCorrection(g_rgnRaw[k], 4912, 5079, 1700, TBC_NOSHIFTING);
  NewTimeBaseCorrection(g_rgnRaw[k], 4624, 4847, 1700);
    /*
    NegShiftSubColumn(g_rgnRaw[k], 4630, 261-252+1, 252);
    NegShiftSubColumn(g_rgnRaw[k], 4649, 1, 252);
    NegShiftSubColumn(g_rgnRaw[k], 4792, 1, 252);
    NegShiftSubColumn(g_rgnRaw[k], 4901+1, 1, 252);
    NegShiftSubColumn(g_rgnRaw[k], 4994+1, 1, 252);
    */

    //
    //  Repair bit sync damage in k=2 transmission
    //
    k=2;
    BitShiftPixels(g_rgnRaw[k], +3, 4548, 4549, 1090-1008, 252);    // blue panorama
    BitShiftPixels(g_rgnRaw[k], +3, 4549, 4550, 1090-1008, 252);
    BitShiftPixels(g_rgnRaw[k], +3, 4549, 4550, 1065-1008, 1081-1008);
    BitShiftPixels(g_rgnRaw[k], -1, 4594, 4595, 938-756, 252);      // second clear panorama
    BitShiftPixels(g_rgnRaw[k], -1, 4612, 4613, 956-756, 252);
    BitShiftPixels(g_rgnRaw[k], +3, 4617, 4618, 1887%252, 1970%252);
    BitShiftPixels(g_rgnRaw[k], -4, 4621, 4622, 53, 252);           // scan reversal region
    BitShiftPixels(g_rgnRaw[k], -1, 4622, 4623, 811-756, 930-756);  // recheck this
    BitShiftPixels(g_rgnRaw[k], -1, 4633, 4636, 955-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4639, 4641, 772-756, 971-756);
    BitShiftPixels(g_rgnRaw[k], +3, 4644, 4647, 1956%252, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4647, 4648, 1219-1008, 1233-1008);
    BitShiftPixels(g_rgnRaw[k], -1, 4648, 4649, 1219-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4649, 4651, 988-756, 873-756);
    BitShiftPixels(g_rgnRaw[k], -4, 4650, 4652, 873-756, 81);       // ending is uncertain
    BitShiftPixels(g_rgnRaw[k], +3, 4651, 4653, 1180-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -3, 4657, 4662, 398-252, 252);
    BitShiftPixels(g_rgnRaw[k], +3, 4672, 4676, 1143-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -4, 4692, 4694, 1167-1008, 222);
    BitShiftPixels(g_rgnRaw[k], +5, 4693, 4696, 222, 252);
    BadPixels(g_rgnRaw[k], 4795, 4796);                             // secondary damage
    BadPixels(g_rgnRaw[k], 4699, 4700);                             // secondary damage
    BitShiftPixels(g_rgnRaw[k], -3, 4702, 4709, 1078-1008, 252);
    BadPixels(g_rgnRaw[k], 4720, 4721);                             // secondary damage
    BadPixels(g_rgnRaw[k], 4724, 4725, 544-504, 252);
    BitShiftPixels(g_rgnRaw[k], +2, 4733, 4739, 970-756, 252);
    BadPixels(g_rgnRaw[k], 4737, 4738);                             // secondary damage
    BitShiftPixels(g_rgnRaw[k], -1, 4739, 4740, 1156-1008, 252);
    BadPixels(g_rgnRaw[k], 4747, 4748, 1063-1008, 252);             // secondary damage
    BitShiftPixels(g_rgnRaw[k], -1, 4768, 4770, 924-756, 873-756);
    BitShiftPixels(g_rgnRaw[k], -2, 4769, 4771, 664-504, 1056-1008);
    BitShiftPixels(g_rgnRaw[k], -3, 4774, 4775, 1088-1008, 1135-1008);
    BitShiftPixels(g_rgnRaw[k], -3, 4774, 4776, 399-252, 252);
    BitShiftPixels(g_rgnRaw[k], -3, 4784, 4786, 1109-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4786, 4787, 0, 252);
    BitShiftPixels(g_rgnRaw[k], -3, 4798, 4804, 1055-1008, 252);
    BitShiftPixels(g_rgnRaw[k], +1, 4809, 4812, 1102-1008, 602-504);
    BitShiftPixels(g_rgnRaw[k], -2, 4811, 4819, 602-504, 252);
    BitShiftPixels(g_rgnRaw[k], +3, 4822, 4827, 1087-1008, 252);
    BitShiftPixels(g_rgnRaw[k], -3, 4842, 4844, 1080-1008, 1080-1008);
    BitShiftPixels(g_rgnRaw[k], -3, 4856, 4859, 490-252, 92);
    BitShiftPixels(g_rgnRaw[k], -4, 4858, 4871, 92, 799-756);
    BitShiftPixels(g_rgnRaw[k], -1, 4870, 4873, 799-756, 839-756);
    BitShiftPixels(g_rgnRaw[k], -3, 4874, 4878, 262-252, 252);
    BadPixels(g_rgnRaw[k], 4883, 4884, 1051-1008, 252);         // zipper noise
    BadPixels(g_rgnRaw[k], 4884, 4885, 0, 1092-1008);           // secondary damage?
    BitShiftPixels(g_rgnRaw[k], -3, 4889, 4892, 1025-1008, 1101-1008);
    BadPixels(g_rgnRaw[k], 4895, 4900, 0, 805-756);             // zipper noise
    BitShiftPixels(g_rgnRaw[k], -1, 4899, 4901, 805-756, 252);
    BadPixels(g_rgnRaw[k], 4901, 4905, 0, 1030-1008);           // zipper noise
    BitShiftPixels(g_rgnRaw[k], -1, 4911, 4914, 1212-1008, 1117-1008);
    BitShiftPixels(g_rgnRaw[k], -3, 4924, 4929, 361-252, 373-252);
    BitShiftPixels(g_rgnRaw[k], -1, 4934, 4936, 818-756, 252);
    //4957-4960 is very noisy, maybe secondary damage
    BitShiftPixels(g_rgnRaw[k], -1, 4976, 4978, 807-756, 882-756);
    BitShiftPixels(g_rgnRaw[k], -2, 4977, 4978, 882-756, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 4979, 4980, 870-756, 252);
    //4976, 4977 secondary damage
    BitShiftPixels(g_rgnRaw[k], -1, 4980, 4981, 0, 252);    // a bad section
    BitShiftPixels(g_rgnRaw[k], -1, 4991, 4995, 957-756, 875-756);
    BitShiftPixels(g_rgnRaw[k], -1, 5022, 5025, 840-756, 1150-1008);
    BitShiftPixels(g_rgnRaw[k], -1, 5033, 5044, 0, 252);
    BitShiftPixels(g_rgnRaw[k], -3, 5044, 5054, 0, 252);
    BitShiftPixels(g_rgnRaw[k], -4, 5054, 5056, 0, 252);
    //5055-5057 ?
    BitShiftPixels(g_rgnRaw[k], +4, 5057, 5061, 2169%252, 252);
    BitShiftPixels(g_rgnRaw[k], -1, 5065, 5076, 1153-1008, 920-756);
    //some 0 and some -1 shift inbetween here
    BitShiftPixels(g_rgnRaw[k], -1, 5078, 5080, 1031-1008, 857-756);
    BitShiftPixels(g_rgnRaw[k], -3, 5080, 5088, 408-252, 252);
    //a telemetry burst followed by a couple unshifted lines
    //BitShiftPixels(g_rgnRaw[k], -1, 5098, 5104, 884-756, 252);  // before shiftcolumns !!!
    /* do this over again after bitsync errors are fixed
    //ShiftColumns(g_rgnRaw[2], 4810, -1, 94);  // no, k=1 needs to shift 1 at 4813
    ShiftColumns(g_rgnRaw[3], 4775, 6, 0);
    ShiftColumns(g_rgnRaw[2], 4977, 1, 249);
    ShiftColumns(g_rgnRaw[2], 5045, -1);
    */
    //BitSyncImage(g_rgnRaw[k], 4000, 5320);
    //TimeBaseCorrection(g_rgnRaw[k], 4624, 5105, 1700, TBC_NOSHIFTING);
  NewTimeBaseCorrection(g_rgnRaw[k], 4624, 4847, 1700);
    /*
    NegShiftSubColumn(g_rgnRaw[k], 4649, 1, 252);
    ShiftSubColumn(g_rgnRaw[k], 4650, 1, 252);
    NegShiftSubColumn(g_rgnRaw[k], 4651, 1, 252);
    ShiftSubColumn(g_rgnRaw[k], 4786, 1, 252);
    ShiftSubColumn(g_rgnRaw[k], 4659, 1, 252);
    ShiftSubColumn(g_rgnRaw[k], 4660, 1, 252);
    ShiftSubColumn(g_rgnRaw[k], 4661, 1, 252);
    NegShiftSubColumn(g_rgnRaw[k], 4694, 1, 252);
    */

    //MaskBurst(4846, 4856+1, 462%252, 404%252);
    //MaskBurst(5085, 5093+1);

    //TimeBaseCorrection(g_rgnRaw[0], 4362, 4369, 3680);
    //TimeBaseCorrection(g_rgnRaw[1], 4408, 4540, 3680, TBC_NOSHIFTING);  // blue panorama
    //TimeBaseCorrection(g_rgnRaw[2], 4379, 4518, 3680, TBC_NOSHIFTING);
    //TimeBaseCorrection(g_rgnRaw[3], 4379, 4518, 3680, TBC_NOSHIFTING);

    //Display3(g_rgnRaw[1], g_rgnRaw[0], g_rgnBrown, 1439, 1428, 307-252, 252);
}
//
//  Old versions (delete after new versions are stable)
//
void
AlignHorizon14II()
{
    int k;

    //TimeBaseCorrection(g_rgnRaw[1], 4627, 4644, 1700, TBC_NOSHIFTING);  // lefthand
    //TimeBaseCorrection(g_rgnRaw[2], 4637, 4644, 1700, TBC_NOSHIFTING);  // horizion of
    //TimeBaseCorrection(g_rgnRaw[3], 4637, 4644, 1700, TBC_NOSHIFTING);  // last clear
    //RotateColumns(g_rgnRaw[1], 4630, 4631, -1, ROT_SHIFT);
    //RotateColumns(g_rgnRaw[1], 4631, 4634, 91, ROT_HELICAL);
    for (k = 1; k < 4; k++) {
        //RotateColumns(g_rgnRaw[k], 4640, 4642, -38, ROT_HELICAL);
        //RotateColumns(g_rgnRaw[k], 4642, 4643, 1, ROT_SHIFT);
    }
    ShiftColumns(g_rgnRaw[2], 4810, -1, 94);
    ShiftColumns(g_rgnRaw[3], 4775, 6, 0);
    /*
    BadPixels(g_rgnRaw[1], 4692, 4697, 663-504);    // gray noise
    BadPixels(g_rgnRaw[1], 4733, 4739, 708-504);
    BadPixels(g_rgnRaw[1], 4809, 4818, 598-504);
    BadPixels(g_rgnRaw[1], 4858, 4871, 522-504);
    BadPixels(g_rgnRaw[3], 4858, 4871, 522-504);
    BadPixels(g_rgnRaw[1], 4897, 4902, 674-504);
    BadPixels(g_rgnRaw[3], 4897, 4902, 674-504);
    BadPixels(g_rgnRaw[2], 4873, 4878, 765-756, 762-756);
    BadPixels(g_rgnRaw[2], 4855, 4872, 993-756, 842-756);
    */
    ShiftColumns(g_rgnRaw[2], 4977, 1, 249);
    /*
    TimeBaseCorrection(g_rgnRaw[1], 4855, 4872, 1700, TBC_NOSHIFTING);  // middle of
    TimeBaseCorrection(g_rgnRaw[2], 4855, 4872, 1700, TBC_NOSHIFTING);  // last clear
    TimeBaseCorrection(g_rgnRaw[3], 4855, 4872, 1700, TBC_NOSHIFTING);  // panorama
    TimeBaseCorrection(g_rgnRaw[1], 4831, 4846, 1700, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[2], 4831, 4846, 1700, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[3], 4831, 4846, 1700, TBC_NOSHIFTING);
    */
    //BadPixels(g_rgnRaw[2], 5033, 5061);
    ShiftColumns(g_rgnRaw[2], 5045, -1);
    //BadPixels(g_rgnRaw[1], 5078, MAXLINES, 620-504);
    //BadPixels(g_rgnRaw[2], 5104, MAXLINES, 814-756);
}

//
//  The second clear panorama in 14-II is the only extended view of the lefthand
//  horizon.  Recover as much as possible from noisy transmission near end of mission.
//
//  Some of this data is bit shifted by three:
//
//  D 223 203 231 239 111110110_ 110100110_ 111001110_ 111101110_
//  C 475 473 476 477 1101101110 1001101111 0011101111 1011101110
//  D 243 243 219 235 110011110_ 110011110_ 110110110_ 110101110_
//  C 478 478 475 413 0111101110 0111101110 1101101110 1011100111
//
/*
void
RepairHorizon14II()
{
    int k;

    BitSyncImage(g_rgnRaw[3], 4560, 5320);
    //
    //  After about 4670, we can do correlation testing against the first
    //  panorama at j - 3897.  Before that, use neighboring lines or clean
    //  lines in differing transmissions of the last panorama.
    //
    for (k = 1; k < 2; k++) {
        ResyncBits(g_rgnRaw[k], 4633, 67+91, 252, -6); // recode this
        ResyncBits(g_rgnRaw[k], 4634, 0, 252, 3);
        ResyncBits(g_rgnRaw[k], 4635, 1, 252, 3);
        BadPixels(g_rgnRaw[k], 4629, 4630, 13, 252);
    }
    for (k = 2; k < 4; k++) {
        ResyncBits(g_rgnRaw[k], 4633, 198, 252, -1);
        ResyncBits(g_rgnRaw[k], 4634, 4, 252, -1);
        ResyncBits(g_rgnRaw[k], 4635, 4, 252, 9);
        ResyncBits(g_rgnRaw[k], 4644, 192, 251, 3);
        ResyncBits(g_rgnRaw[k], 4645, 2, 252, 3);
        ResyncBits(g_rgnRaw[k], 4646, 2, 252, 3);
        BadPixels(g_rgnRaw[k], 4624, 4625, 23, 252);
    }
    for (k = 1; k < 4; k++) {
        ResyncBits(g_rgnRaw[k], 4639, 4, 252, -1);
        ResyncBits(g_rgnRaw[k], 4640, 4, 252, -1);
    }
    //
    //  Apparently video sync errors happen before bit-stream sync errors,
    //  so they should be corrected after.
    //
    RotateColumns(g_rgnRaw[1], 4630, 4631, -1, ROT_SHIFT);
    RotateColumns(g_rgnRaw[1], 4631, 4634, 91, ROT_HELICAL);
    RotateColumns(g_rgnRaw[1], 4634, 4636, 42, ROT_HELICAL);
    for (k = 2; k < 4; k++) {
        RotateColumns(g_rgnRaw[k], 4632, 4633, 91, ROT_SHIFT);
    }
    for (k = 1; k < 4; k++) {
        RotateColumns(g_rgnRaw[k], 4640, 4642, -38, ROT_HELICAL);
        RotateColumns(g_rgnRaw[k], 4642, 4643, 1, ROT_SHIFT);
    }
    //AutoBitShift(g_rgnRaw[1], 4633, g_rgnRaw[1], 4634, 0, 160);
    return;

    for (k = 1; k < 4; k++) {
        ResyncBits(g_rgnRaw[k], 4639, 4, 252, -1);
        ResyncBits(g_rgnRaw[k], 4640, 4, 252, -1);
    }
    for (k = 1; k < 2; k++) {
        BadPixels(g_rgnRaw[k], 4629, 4630, 13, 252);    // Resync 4 or -2, -4
        ResyncBits(g_rgnRaw[k], 4634, 0, 252, 3);
        ResyncBits(g_rgnRaw[k], 4635, 1, 252, 3);
    }
    for (k = 2; k < 4; k++) {
        BadPixels(g_rgnRaw[k], 4624, 4625, 23, 252);
        RotateColumns(g_rgnRaw[k], 4632, 4633, 91, ROT_SHIFT);
        BadPixels(g_rgnRaw[k], 4632, 4633, 793-756, 862-756);
        BadPixels(g_rgnRaw[k], 4633, 4634, 199, 252);
        ResyncBits(g_rgnRaw[k], 4634, 4, 252, -1);
        ResyncBits(g_rgnRaw[k], 4635, 4, 252, 9);
        ResyncBits(g_rgnRaw[k], 4644, 192, 251, 3);
        ResyncBits(g_rgnRaw[k], 4645, 2, 252, 3);
        ResyncBits(g_rgnRaw[k], 4646, 2, 252, 3);
        ResyncBits(g_rgnRaw[k], 4650, 2, 116, 9);
        ResyncBits(g_rgnRaw[k], 4650, 116, 927-756, 3);
        ResyncBits(g_rgnRaw[k], 4650, 927-756, 252, 6);

        //ResyncBits(g_rgnRaw[k], 4651, 2, 252, -1);
        //ResyncBits(g_rgnRaw[k], 4652, 2, 252, -1);

        //ResyncBits(g_rgnRaw[k], 4734, 2, 252, 5);
        //ResyncBits(g_rgnRaw[k], 4735, 2, 252, 5);
        //ResyncBits(g_rgnRaw[k], 4736, 2, 252, 5);
        //ResyncBits(g_rgnRaw[k], 4737, 2, 252, 5);
        //ResyncBits(g_rgnRaw[k], 4738, 2, 252, 5);
    }
    RotateColumns(g_rgnRaw[1], 4630, 4631, -1, ROT_SHIFT);
    RotateColumns(g_rgnRaw[1], 4631, 4634, 91, ROT_HELICAL);
    RotateColumns(g_rgnRaw[1], 4634, 4636, 42, ROT_HELICAL);
    for (k = 1; k < 4; k++) {
        RotateColumns(g_rgnRaw[k], 4640, 4642, -38, ROT_HELICAL);
        RotateColumns(g_rgnRaw[k], 4642, 4643, 1, ROT_SHIFT);
    }
}
*/