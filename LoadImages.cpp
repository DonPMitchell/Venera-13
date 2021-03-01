#include "stdafx.h"
//
//  Early image processing of undecoded 9-bit telemetry
//  D. P. Mitchell  03/25/2003.
//
//  1. Read raw telemetry transmissions
//  2. Align transmissions, fix up video sync (time-base correction)
//  3. Mask out telemetry bursts and obvious bad data
//  4. Mask out statistically rejected noise
//  5. Re-decode pixel values in some anomalous transmissions or regions
//  6. Assemble best master version from multiple transmissions
//
#pragma intrinsic(abs)
#include "Venus.h"
#include "ImageFile.h"
#include "ImageProcessing.h"
#include "Recoding.h"

unsigned short  g_rgnRaw[4][MAXLINES][252];         // versions of raw telemetry
unsigned short  g_rgnFlags[4];                      // version info
unsigned short  g_rgnBrown[MAXLINES][252];          // 8-bit Brown Univ. image
unsigned short  g_rgnMasterV13C1[MAXLINES][252];    // Venera 13, Camera 1
unsigned short  g_rgnMasterV13C2[MAXLINES][252];    // Venera 13, Camera 2
unsigned short  g_rgnMasterV14C1[MAXLINES][252];    // Venera 14, Camera 1
unsigned short  g_rgnMasterV14C2[MAXLINES][252];    // Venera 14, Camera 2
unsigned short  g_rgnTest1[MAXLINES][252];          // Visualizeations
unsigned short  g_rgnTest2[MAXLINES][252];
unsigned short  g_rgnBackup[4][MAXLINES][252];      // Backup of anomalous images
unsigned short  g_rgnMisPlaced[400][252];           // 14-II data form a 14-I tape
char            g_rgnSource[MAXLINES][252];         // Origin of master pixels
int             g_rgkBrown[2];                      //REVIEW: delete
unsigned short  g_rgnTelemetry[128][4096];          // Telemetry-burst data
int             g_jTelemetry = 0;
double          g_fShiftLimit = 0.85;               // Repair zero shifts without check
int             g_kVerbose = -1;
//
//  Video front porch [4] == 33152, 32896, 33664, 33408 for the four cameras
//
unsigned short  g_rgnPorch[5] = { 0, 51072, 11136, 26880, 33408};

//
//  Read raw data into temperary buffers.  This is the 16-bit data from the Russians,
//  containing 9-bit image data in high order, and 0 or 127 in the low order.
//
//  Note, Photoshop mangled the Russian and Brown data when it was used to convert
//  into 16-bit raw format.  It may dither, round, extend and scale bits, so it was
//  taken out of the loop except to convert 8-bit bmp into 8-bit raw BU data.  The
//  original Russian TIFF files are now loaded directly.
//
static int
LoadRawImage(unsigned short rgnRaw[MAXLINES][252], char *szFile, int jFirst)
{
    int i, j, nLines, jShift;

    for (j = 0; j < MAXLINES; j++) {
        for (i = 0; i < 252; i++) {
            rgnRaw[j][i] = BADPIXEL | BURSTPIXEL;
            g_rgnTest1[j][i] = BADPIXEL;
        }
    }
    nLines = ReadTIFF(szFile, g_rgnTest2);
    for (j = 0; j < nLines; j++) {
        if (j+jFirst < 0 || j+jFirst >= MAXLINES)
            continue;
        for (i = 0; i < 252; i++)
            rgnRaw[j+jFirst][i] = g_rgnTest2[j][i];
    }
    for (j = 0; j < nLines; j++)
        for (i = 0; i < 252; i++) {
            g_rgnTest1[jFirst + j][i] = rgnRaw[jFirst + j][i] << 9;
            rgnRaw[jFirst + j][i] &= NINEBITS;
        }
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
//  Set and test pixel flags.
//
void
SetBad(unsigned short &nPixel)
{
    nPixel = nPixel | BADPIXEL;
}

void
SetBurst(unsigned short &nPixel)
{
    nPixel = nPixel | BURSTPIXEL | BADPIXEL;    // e.g., telemetry burst pixels
}

void
SetRecode(unsigned short &nPixel)
{
    nPixel = nPixel | RECODEPIXEL | BADPIXEL;   // failed to recode
}

void
BadPixels(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
          int iFirst, int iLimit)
{
    int ij, ijFirst, ijLimit;
    unsigned short *pn;

    ijFirst = iFirst + 252*jFirst;
    ijLimit = iLimit + 252*(jLimit - 1);
    pn = rgnData[0];
    for (ij = ijFirst; ij < ijLimit; ij++)
        //SetBad(pn[ij]);
        pn[ij] = 1;
}

void
BurstPixels(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
            int iFirst = 0, int iLimit = 252)
{
    int ij, ijFirst, ijLimit;
    unsigned short *pn;

    ijFirst = iFirst + 252*jFirst;
    ijLimit = iLimit + 252*(jLimit - 1);
    pn = rgnData[0];
    for (ij = ijFirst; ij < ijLimit; ij++)
        SetBurst(pn[ij]);
}

int
IsBad(unsigned short nPixel)
{
    return ((nPixel & BADPIXEL) != 0);
}

int
IsBurst(unsigned short nPixel)
{
    return ((nPixel & BURSTPIXEL) != 0);
}

int
IsRecode(unsigned short nPixel)
{
    return ((nPixel & RECODEPIXEL) != 0);
}

int
IsGood(unsigned short nPixel)
{
    return ((nPixel & BADPIXEL) == 0);
}

int
IsPorch(unsigned short rgnData[])
{
    int i;

    for (i = 0; i < 5; i++)
        if (rgnData[i] != g_rgnPorch[i])    // video front porch
            return 0;
    return 1;
}

//
//  Routines to manually or automatically align the telemetry versions.
//
void
VisualizeDifferences(int nCopies)
{
    int i, j, k, nDiff;

    for (j = 0; j < MAXLINES; j++) {
        for (i = 0; i < 252; i++) {
            nDiff = 0x00FF;
            for (k = 0; k < nCopies; k++) {
                if (IsBad(g_rgnRaw[k][j][i]))
                    continue;
                if (nDiff == 0x00FF)
                    nDiff = g_rgnRaw[k][j][i];
                else if (g_rgnRaw[k][j][i] != nDiff)
                    nDiff = 0xFF00;
            }
            g_rgnTest1[j][i] = nDiff;
        }
    }
}

void
ShiftColumns(unsigned short rgnData[MAXLINES][252], int jFirst, int nShift,
             int iFirst)
{
    int ij, ijFirst, ijShift;
    unsigned short *pn;

    ijFirst = iFirst + 252*jFirst;
    ijShift = 252*nShift;
    pn = rgnData[0];
    if (ijShift < 0) {
        for (ij = ijFirst; ij < MAXLINES*252; ij++) {
            if (ij + ijShift < 0)
                continue;
            pn[ij + ijShift] = pn[ij];
            //SetBad(pn[ij]);
            pn[ij] = BADPIXEL|BURSTPIXEL;
        }
    } else {
        for (ij = MAXLINES*252-1; ij >= ijFirst; --ij) {
            if (ij + ijShift >= MAXLINES*252)
                continue;
            pn[ij + ijShift] = pn[ij];
            //SetBad(pn[ij]);
            pn[ij] = BADPIXEL|BURSTPIXEL;
        }
    }
}

void
RotateColumns(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit, int nRot,
              int bShift, int kBrown)
{
    int i, j, jStart, jStop, jDelta;
    unsigned short rgn[252];

    if (nRot >= 0) {        // Helical rotations must be done in proper order
        jStart = jFirst;
        jStop = jLimit;
        jDelta = 1;
    } else {
        jStart = jLimit-1;
        jStop = jFirst-1;
        jDelta = -1;
    }
    for (j = jStart; j != jStop; j += jDelta) {
        for (i = 0; i < 252; i++)
            rgn[i] = rgnData[j][i];
        for (i = 0; i < 252; i++)
            rgnData[j][i] = rgn[(i + nRot + 252) % 252];
        if (bShift != ROT_ROTATE) {
            if (nRot > 0)
                for (i = 252 - nRot; i < 252; i++)
                    //SetBad(rgnData[j][i]);
                    rgnData[j][i] = BADPIXEL|BURSTPIXEL;
            else
                for (i = 0; i < -nRot; i++)
                    //SetBad(rgnData[j][i]);
                    rgnData[j][i] = BADPIXEL|BURSTPIXEL;
        }
        //
        //  Helical rotations repair video-sync errors.
        //
        if (bShift == ROT_HELICAL && nRot > 0) {
            for (i = 0; i < nRot; i++)
                if (rgnData[j-1][i + 252 - nRot] & BADPIXEL)
                    rgnData[j-1][i + 252 - nRot] = rgn[i];
        }
        if (bShift == ROT_HELICAL && nRot < 0) {
            for (i = 0; i < -nRot; i++)
                if (rgnData[j+1][i] & BADPIXEL)
                    rgnData[j+1][i] = rgn[252 + nRot + i];
        }
    }
}
//
//  Main routine to do all loading and early image processing.
//  Images are loaded into g_rgnRaw[] roughly in order of quality.
//
extern void DebugImage(unsigned short rgnDiff[MAXLINES][252],
                        unsigned short rgnMaster[MAXLINES][252],
                        char *szTestImage);
extern void DebugFourImages(unsigned short rgn1[MAXLINES][252], 
                        unsigned short rgn2[MAXLINES][252],
                        unsigned short rgn3[MAXLINES][252],
                        unsigned short rgn4[MAXLINES][252],
                        char *szTestImage);

void
LoadImages()
{
    int i, j;

    printf("II. Load, Align and Mask Venera 13 & 14 Images\n");
    //goto V14II;
    //
    //  V-13-I
    //
    printf("  A. Venera 13, Camera 1\n");
    LoadRawImage(g_rgnRaw[0], "B-13-5.tif",    0);
    LoadRawImage(g_rgnRaw[1], "B-13-3.tif",   47);
    LoadRawImage(g_rgnRaw[2], "B-13-4.tif", 5145);  // Positioned based on image alignment
    LoadBrownImage(g_rgnBrown, "B-13-I-BrownUniv.raw", 242);  // 13-3 and 13-4
    g_rgnFlags[0] = FLAGS_TELEMETRY;
    g_rgnFlags[1] = FLAGS_TELEMETRY | FLAGS_BROWN;
    g_rgnFlags[2] = FLAGS_TELEMETRY | FLAGS_BROWN;
    g_rgnFlags[3] = FLAGS_NODATA;
    g_rgkBrown[0] = 1;
    g_rgkBrown[1] = 2;
    g_rgnPorch[4] = 33152;
    FixBrown();
    BadPixels(g_rgnRaw[3], 0, MAXLINES);
    BurstPixels(g_rgnRaw[0], 0, 242);                       // line 241 stretched.
    BurstPixels(g_rgnRaw[1], 0, 242);
    BadPixels(g_rgnRaw[0], 2066, 2067, 128, 229);           // noise in blue panorama
    ShiftSubColumn(g_rgnRaw[1], 2112, 86+1);
    BadPixels(g_rgnRaw[1], 307, 308);                       // subtle miscoding, might be missed so nuke it
    BadPixels(g_rgnRaw[1], 621, 623);                       // TBC doesn't fix this
    BurstPixels(g_rgnRaw[2], 5145, 5151);                   // Noise just before the mis-synced bit
    ShiftColumns(g_rgnRaw[2], 6870, 58+3, 62);              // Missing lines in last red,green images
    RotateColumns(g_rgnRaw[2], 6928+3, 6929+3, 132-41, ROT_SHIFT);
    ShiftSubColumn(g_rgnRaw[2], 6929+3, 5);                 // These are rare in 13-I
    BrownReplace(2, 6929+5, 251);
    ShiftSubColumn(g_rgnRaw[2], 6283, 1345%252);
    BrownReplace(2, 6283, 251);
    ShiftColumns(g_rgnRaw[2], 5145, -35);                   // The scrambled section
    ShiftColumns(g_rgnRaw[2], 5119, 35);
    AutoMaskBurst(g_rgnRaw[1], 713, 0);
    AutoMaskBurst(g_rgnRaw[2], 5290, 1);
    MaskBurst(2874, 2885, 334%252, 278%252);
    MaskBurst(6485, 6497, 874%252, 818%252);
    MaskBurst(7205+3, 7216+3, 774%252, 252);
    TimeBaseCorrection(g_rgnRaw[2], 6872+62, 7145+62, -1);      // last green image (4520 ramp)
    TimeBaseCorrection(g_rgnRaw[2], 7164+62, 7388+62, 4850);      // last blue image
    BrownRepair(2, 7142, 1, 1407%252);
    //SpecialV13C1();
    g_kVerbose = 0;
    MasterVersion(g_rgnMasterV13C1, 3, SpecialV13C1, 6700);
    DebugImage(g_rgnTest1, g_rgnMasterV13C1, "Test1.bmp");
    WriteTIFF("..//ProjectResults//Master1.tif", 7449-242, &g_rgnMasterV13C1[242]);
    return;
    //
    //  V-13-II
    //
V13II:
    printf("  B. Venera 13, Camera 2\n");
    LoadRawImage(g_rgnRaw[0], "B-13-1.tif",    0);
    LoadRawImage(g_rgnRaw[1], "B-13-2.tif", 5036-1);          // Positioned by telemetry burst
    LoadRawImage(g_rgnRaw[2], "B-13-6.tif",  190-285);          // Anomalous image
    LoadBrownImage(g_rgnBrown, "B-13-II-BrownUniv.raw", 208); // 13-1 and 13-2
    g_rgnFlags[0] = FLAGS_TELEMETRY | FLAGS_BROWN;
    g_rgnFlags[1] = FLAGS_TELEMETRY | FLAGS_BROWN;
    g_rgnFlags[2] = FLAGS_TELEMETRY | FLAGS_ANOMALOUS;
    g_rgnFlags[3] = FLAGS_NODATA;
    g_rgkBrown[0] = 0;
    g_rgkBrown[1] = 1;
    g_rgnPorch[4] = 32896;
    FixBrown();
    BadPixels(g_rgnRaw[3], 0, MAXLINES);
    RotateColumns(g_rgnRaw[0], 208, 209, 4, ROT_SHIFT);     // The first line seems fixable
    RotateColumns(g_rgnRaw[2], 208, 209, 4, ROT_SHIFT);
    BurstPixels(g_rgnRaw[0], 0, 208);
    BurstPixels(g_rgnRaw[2], 0, 208);
    BurstPixels(g_rgnRaw[1], 5035, 5036);
    ShiftColumns(g_rgnRaw[2], 2843, -1);                // Extra line within a telemetry burst
    BurstPixels(g_rgnRaw[0], 3134, 3135);
    BurstPixels(g_rgnRaw[0], 3991, 3992);
    BurstPixels(g_rgnRaw[2], 1057, 1058, 119);
    BurstPixels(g_rgnRaw[2], 1058, 1059);
    BurstPixels(g_rgnRaw[2], 1071, 1072, 134);
    BurstPixels(g_rgnRaw[2], 1072, 1073);
    BurstPixels(g_rgnRaw[2], 1113, 1114, 59);
    BurstPixels(g_rgnRaw[2], 1114, 1115);
    BurstPixels(g_rgnRaw[2], 1350, 1352);    
    AutoMaskBurst(g_rgnRaw[0], 680, 0);
    AutoMaskBurst(g_rgnRaw[1], 5248, 1);
    MaskBurst(2116, 2126, 34, 252);
    MaskBurst(926, 927);               // repair an autoburst
    BurstPixels(g_rgnRaw[2], 4160, MAXLINES);
    ShiftSubColumn(g_rgnRaw[0], 1321, 302+1-252);
    BrownReplace(0, 1321, 251);
    ShiftSubColumn(g_rgnRaw[0], 2766, 337+1-252);
    BrownReplace(0, 2766, 251);
    ShiftSubColumn(g_rgnRaw[0], 2789, 323+1-252);
    BrownReplace(0, 2766, 251);
    ShiftSubColumn(g_rgnRaw[2], 1442, 1+1);
    ShiftSubColumn(g_rgnRaw[2], 2687, 63+1);
    ShiftSubColumn(g_rgnRaw[2], 3407, 802+1-756);
    ShiftSubColumn(g_rgnRaw[2], 3590, 927+1-756);   // shifted by 2
    ShiftSubColumn(g_rgnRaw[2], 3797, 929+1-756);
    ShiftSubColumn(g_rgnRaw[2], 3920, 28+1);
    ShiftSubColumn(g_rgnRaw[2], 3590, 927+1-756);
    ShiftSubColumn(g_rgnRaw[2], 3797, 929+1-756);
    ShiftSubColumn(g_rgnRaw[2], 3920, 28+1);
    ShiftSubColumn(g_rgnRaw[1], 5092, 66+1);
    BrownReplace(1, 5092, 251);                     //REVIEW: do in MasterVersion()
    ShiftSubColumn(g_rgnRaw[1], 5427, 50+1);        // shifted twice?, zero at 173
    BrownReplace(1, 5427, 251);
    ShiftSubColumn(g_rgnRaw[1], 5118, 149+1);
    BrownReplace(1, 5118, 251);
    ShiftSubColumn(g_rgnRaw[1], 5297, 242+1);
    BrownReplace(1, 5297, 251);
    ShiftSubColumn(g_rgnRaw[1], 6008, 518-504+1);
    BrownReplace(1, 6008, 251);
    BurstPixels(g_rgnRaw[0], 3327, 3328, 43);         // subtle error
    BurstPixels(g_rgnRaw[0], 3375, 3376, 43);
    //SpecialV13C2();
    MasterVersion(g_rgnMasterV13C2, 3, SpecialV13C2, 5720, 2);
    DebugImage(g_rgnTest1, g_rgnMasterV13C2, "Test2.bmp");
    WriteTIFF("..//ProjectResults//Master2.tif", 6515-208, &g_rgnMasterV13C2[208]);
    //
    //  V-14-I
    //
V14I:
    printf("  C. Venera 14, Camera 1\n");
    LoadRawImage(g_rgnRaw[0], "B-14-1.tif",    0);
    LoadRawImage(g_rgnRaw[1], "B-14-6.tif",  396);       // Data from Camera II at the end (?)
    LoadRawImage(g_rgnRaw[2], "B-14-9.tif",  558);          // Anomalous image
    LoadRawImage(g_rgnRaw[3], "B-14-13.tif", 566);          // Anomalous image
    LoadBrownImage(g_rgnBrown, "B-14-I-BrownUniv.raw", 660);  // 14-1
    g_rgnFlags[0] = FLAGS_TELEMETRY | FLAGS_BROWN;
    g_rgnFlags[1] = FLAGS_TELEMETRY;
    g_rgnFlags[2] = FLAGS_TELEMETRY | FLAGS_ANOMALOUS;
    g_rgnFlags[3] = FLAGS_TELEMETRY | FLAGS_ANOMALOUS;
    g_rgkBrown[0] = 0;
    g_rgkBrown[1] = -1;
    g_rgnPorch[4] = 33664;
    FixBrown();
    BurstPixels(g_rgnRaw[0], 0, 660);                     // 659 is stretched and noisy
    BurstPixels(g_rgnRaw[1], 0, 660);
    BurstPixels(g_rgnRaw[2], 0, 660);
    BurstPixels(g_rgnRaw[3], 0, 660);
    BurstPixels(g_rgnRaw[2], 4429, MAXLINES, 187);
    BurstPixels(g_rgnRaw[3], 4429, MAXLINES, 187);
    ShiftColumns(g_rgnRaw[0], 3295, -1);                // extra line within telemetry burst
    ShiftColumns(g_rgnRaw[1], 3295, -1);
    RotateColumns(g_rgnRaw[0], 3300, 3301, 218, ROT_ROTATE);    // fix line after telemetry burst
    RotateColumns(g_rgnRaw[1], 3300, 3301, 218, ROT_ROTATE);
    RotateColumns(g_rgnRaw[2], 3300, 3301, 74, ROT_SHIFT);      // distrupted differently in this transmission
    RotateColumns(g_rgnRaw[3], 3300, 3301, 74, ROT_SHIFT);
    ShiftColumns(g_rgnRaw[0], 4261, 1);                 // missing line within telemetry burst
    ShiftColumns(g_rgnRaw[1], 4261, 1);
    ShiftColumns(g_rgnRaw[0], 4274, -1);                // unknown sync problem during camera reverse
    ShiftColumns(g_rgnRaw[1], 4274, -1);
    NegShiftSubColumn(g_rgnRaw[2], 1652, 87);
    ShiftSubColumn(g_rgnRaw[0], 3590, 183+1);
    for (j = 4745; j < 4745+400; j++)                   
        for (i = 0; i < 252; i++)
            g_rgnMisPlaced[j-4745][i] = g_rgnRaw[1][j][i];  // This data is from camera II
    BadPixels(g_rgnRaw[1], 4745, MAXLINES);
    ShiftColumns(g_rgnRaw[3], 681, 12);                 // repair very-anomalous image
    BadPixels(g_rgnRaw[3], 1671, 1789);                 // buffer outputs same data
    ShiftColumns(g_rgnRaw[3], 1789, 21);
    ShiftColumns(g_rgnRaw[3], 2009, 1);
    ShiftColumns(g_rgnRaw[3], 3618, 5);
    AutoMaskBurst(g_rgnRaw[0], 1130, 0);                // No alignments below here
    MaskBurst(4256, 4267, 980%252);
    //SpecialV14C1();                                 // misc. special repairs
    g_fShiftLimit = 0.90;
    MasterVersion(g_rgnMasterV14C1, 3, SpecialV14C1, 3200, 2);    // 4th version is not helpful
    DebugImage(g_rgnTest1, g_rgnMasterV14C1, "Test3.bmp");
    WriteTIFF("..//ProjectResults//Master3.tif", 5068-660, &g_rgnMasterV14C1[660]);
    return;
    //
    //  V-14-II
    //
V14II:
    printf("  D. Venera 14, Camera 2\n");
    LoadRawImage(g_rgnRaw[0], "B-14-11.tif",  0);
    LoadRawImage(g_rgnRaw[1], "B-14-5.tif", 108);
    LoadRawImage(g_rgnRaw[2], "B-14-3.tif", 132);
    LoadRawImage(g_rgnRaw[3], "B-14-7.tif", 487);
    LoadBrownImage(g_rgnBrown, "B-14-II-BrownUniv.raw", 773); // 14-5
    g_rgnFlags[0] = FLAGS_TELEMETRY;
    g_rgnFlags[1] = FLAGS_TELEMETRY | FLAGS_BROWN;
    g_rgnFlags[2] = FLAGS_TELEMETRY;
    g_rgnFlags[3] = FLAGS_TELEMETRY;
    g_rgkBrown[0] = 1;
    g_rgkBrown[1] = -1;
    g_rgnPorch[4] = 33408;
    FixBrown();
    RotateColumns(g_rgnRaw[0], 773, 774, 469-504, ROT_SHIFT);
    RotateColumns(g_rgnRaw[1], 773, 774, 676-571-3, ROT_SHIFT);
    RotateColumns(g_rgnRaw[2], 773, 774, 898-1008, ROT_SHIFT);
    RotateColumns(g_rgnRaw[3], 773, 774, 1150-1260, ROT_SHIFT);
    BurstPixels(g_rgnRaw[0], 0, 774);   // Initial telemetry, 773 is too stetched
    BurstPixels(g_rgnRaw[1], 0, 774);
    BurstPixels(g_rgnRaw[2], 0, 774);
    BurstPixels(g_rgnRaw[3], 0, 774);
    RotateColumns(g_rgnRaw[1], 774, 775, 23-6, ROT_SHIFT);    
    RotateColumns(g_rgnRaw[0], 787, 788, 4, ROT_SHIFT);             // Stretched top bit
    ShiftColumns(g_rgnRaw[0], 864, 1);
    BadPixels(g_rgnRaw[0], 863, 864, 169);
    BadPixels(g_rgnRaw[2], 1285, 1293);
    BadPixels(g_rgnRaw[3], 1285, 1293);
    TimeBaseCorrection(g_rgnRaw[2], 1293, 1314, 1212, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[3], 1293, 1314, 1212, TBC_NOSHIFTING);
    BadPixels(g_rgnRaw[2], 1297, 1301);
    BadPixels(g_rgnRaw[3], 1297, 1301);
    BadPixels(g_rgnRaw[2], 1314, 1332);
    BadPixels(g_rgnRaw[3], 1314, 1332);
    BadPixels(g_rgnRaw[2], 1433, 1451);
    BadPixels(g_rgnRaw[3], 1433, 1451);
    ShiftColumns(g_rgnRaw[2], 1437, 1);
    ShiftColumns(g_rgnRaw[3], 1437, 1);
    BadPixels(g_rgnRaw[2], 1559, 1576, 234, 11);
    BadPixels(g_rgnRaw[3], 1559, 1576, 234, 11);
    ShiftColumns(g_rgnRaw[2], 1566, 1);
    ShiftColumns(g_rgnRaw[3], 1566, 1);
    BadPixels(g_rgnRaw[2], 2034, 2037);     // masking out very bad stuff gives
    BadPixels(g_rgnRaw[3], 2034, 2037);     // the automatic masking a better chance.
    BadPixels(g_rgnRaw[2], 2085, 2086);
    BadPixels(g_rgnRaw[3], 2085, 2086);
    BadPixels(g_rgnRaw[2], 2124, 2126);
    BadPixels(g_rgnRaw[3], 2124, 2126);
    BadPixels(g_rgnRaw[2], 2244, 2246);
    BadPixels(g_rgnRaw[3], 2244, 2246);
    BadPixels(g_rgnRaw[2], 2577, 2580);
    BadPixels(g_rgnRaw[3], 2577, 2580);
    TimeBaseCorrection(g_rgnRaw[2], 2030, 2040, 2030, TBC_NOSHIFTING);  //noisy
    TimeBaseCorrection(g_rgnRaw[3], 2030, 2040, 2030, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[0], 1786, 1788, 1783);
    TimeBaseCorrection(g_rgnRaw[2], 2202, 2205, 2030, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[3], 2202, 2205, 2030, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[0], 2494, 2498, 2493, TBC_NOSHIFTING);
    ShiftColumns(g_rgnRaw[1], 2212, -1);
    ShiftColumns(g_rgnRaw[2], 2690, -1);
    ShiftColumns(g_rgnRaw[3], 2690, -1);
    ShiftColumns(g_rgnRaw[0], 3135, 8);
    BadPixels(g_rgnRaw[0], 3134, 3145, 134);    //REVIEW: timebase correction stuff here
    BadPixels(g_rgnRaw[2], 3341, 3360, 77, 252);
    BadPixels(g_rgnRaw[3], 3341, 3360, 77, 252);
    BadPixels(g_rgnRaw[2], 3364, 3376, 154, 177);
    BadPixels(g_rgnRaw[3], 3364, 3376, 154, 177);
    BadPixels(g_rgnRaw[2], 3387, 3410, 213, 252);
    BadPixels(g_rgnRaw[3], 3387, 3410, 213, 252);
    BadPixels(g_rgnRaw[2], 3492, 3515, 122, 770-756);
    BadPixels(g_rgnRaw[3], 3492, 3515, 122, 770-756);
    ShiftColumns(g_rgnRaw[2], 3500, -1);
    ShiftColumns(g_rgnRaw[3], 3500, -1);
    BadPixels(g_rgnRaw[2], 3895, 3902);
    BadPixels(g_rgnRaw[3], 3895, 3902);
    ShiftColumns(g_rgnRaw[2], 4135, -1);
    ShiftColumns(g_rgnRaw[3], 4135, -1);
    //TimeBaseCorrection(g_rgnRaw[0], 4362, 4369, 3680);
    ShiftColumns(g_rgnRaw[0], 4369, 1);
    BadPixels(g_rgnRaw[2], 4423, 4426, 13, 82-56);  // this is destroyed, zipper noise
    BadPixels(g_rgnRaw[3], 4423, 4426, 13, 82-56);
    BadPixels(g_rgnRaw[1], 4398, 4408, 115);
    //TimeBaseCorrection(g_rgnRaw[1], 4408, 4540, 3680, TBC_NOSHIFTING);  // blue panorama
    //TimeBaseCorrection(g_rgnRaw[2], 4379, 4518, 3680, TBC_NOSHIFTING);
    //TimeBaseCorrection(g_rgnRaw[3], 4379, 4518, 3680, TBC_NOSHIFTING);
    BadPixels(g_rgnRaw[0], 4541, MAXLINES);
    BadPixels(g_rgnRaw[1], 4542, 4572, 153);    // not bitsync repairable
    /*
    //
    //  Alignment of last clear panorama is done in AlignHorizon14II().
    //
    TimeBaseCorrection(g_rgnRaw[1], 4627, 4644, 1700, TBC_NOSHIFTING);  // lefthand
    TimeBaseCorrection(g_rgnRaw[2], 4637, 4644, 1700, TBC_NOSHIFTING);  // horizion of
    TimeBaseCorrection(g_rgnRaw[3], 4637, 4644, 1700, TBC_NOSHIFTING);  // last clear
    ShiftColumns(g_rgnRaw[2], 4810, -1, 94);
    ShiftColumns(g_rgnRaw[3], 4775, 6, 0);
    BadPixels(g_rgnRaw[1], 4692, 4697, 663-504);    // gray noise
    BadPixels(g_rgnRaw[1], 4733, 4739, 708-504);
    BadPixels(g_rgnRaw[1], 4809, 4818, 598-504);
    BadPixels(g_rgnRaw[1], 4858, 4871, 522-504);
    BadPixels(g_rgnRaw[3], 4858, 4871, 522-504);
    BadPixels(g_rgnRaw[1], 4897, 4902, 674-504);
    BadPixels(g_rgnRaw[3], 4897, 4902, 674-504);
    BadPixels(g_rgnRaw[2], 4873, 4878, 765-756, 762-756);
    BadPixels(g_rgnRaw[2], 4855, 4872, 993-756, 842-756);
    ShiftColumns(g_rgnRaw[2], 4977, 1, 249);
    TimeBaseCorrection(g_rgnRaw[1], 4855, 4872, 1700, TBC_NOSHIFTING);  // middle of
    TimeBaseCorrection(g_rgnRaw[2], 4855, 4872, 1700, TBC_NOSHIFTING);  // last clear
    TimeBaseCorrection(g_rgnRaw[3], 4855, 4872, 1700, TBC_NOSHIFTING);  // panorama
    TimeBaseCorrection(g_rgnRaw[1], 4831, 4846, 1700, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[2], 4831, 4846, 1700, TBC_NOSHIFTING);
    TimeBaseCorrection(g_rgnRaw[3], 4831, 4846, 1700, TBC_NOSHIFTING);
    BadPixels(g_rgnRaw[2], 5033, 5061);
    ShiftColumns(g_rgnRaw[2], 5045, -1);
    BadPixels(g_rgnRaw[1], 5078, MAXLINES, 620-504);
    BadPixels(g_rgnRaw[2], 5104, MAXLINES, 814-756);
    */
    //AlignHorizon14II();
    AutoMaskBurst(g_rgnRaw[0], 1245, 0, 4500);
    MaskBurst(4368, 4377+1, 674%252);
    //MaskBurst(4845, 4856+1, 714%252, 656%252);
    //MaskBurst(5085, 5093+1);
    ShiftSubColumn(g_rgnRaw[2], 1138, 877-756+1);
    ShiftSubColumn(g_rgnRaw[2], 1238, 802-756+1);
    ShiftSubColumn(g_rgnRaw[2], 1461, 781-756+1);
    ShiftSubColumn(g_rgnRaw[2], 1463, 805-756);     // double zero
    ShiftSubColumn(g_rgnRaw[2], 1463, 805-756);
    //RepairHorizon14II();
    //SpecialV14C2();
    g_kVerbose = -1;
    MasterVersion(g_rgnMasterV14C2, 3, SpecialV14C2, 4300);
    DebugImage(g_rgnTest1, g_rgnMasterV14C2, "Test4.bmp");
    WriteTIFF("..//ProjectResults//Master4.tif", 5317 - 773, &g_rgnMasterV14C2[773]);
}
