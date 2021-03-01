#include "stdafx.h"
//
//  Process the Venera-9 and Venera-10 panoramas
//  D. P. Mitchell  04/11/2003.
//
//  This is sort of a warm-up for the more complex tasks of Venera 13 and 14.
//  The 13 lines of calibration data are currently missing, it would be nice to
//  have that to correct for the varying automatic gain.
//
#include "Venus.h"
#include "ImageFile.h"
#include "ImageProcessing.h"

extern void InPainting(float *pfImage, char *pnMask, int nHigh, int nWide, float fPhase);

static int s_nLinesV09, s_nLinesV10;
static unsigned short   s_rgnV09[512][128];     // raw images
static unsigned short   s_rgnV10[512][128];
static unsigned short   s_rgnTest[512][128];    // temp space
static float            s_rgfV09[512][128];     // linearized images
static float            s_rgfV10[512][128];
static char             s_rgnV09Mask[512][128]; // missing regions
static char             s_rgnV10Mask[512][128];

static unsigned short   s_rgnRaw2x[256][1024];   // workspaces
static float            s_rgfTemp[512][256];
static char             s_rgnTemp[512][256];
static float            s_rgfLaplace[512][128];
static float            s_rgfTranspose[256][1024];
static float            s_rgfWorkSpace[512];

//
//  Radiometric response functions (signal U to optical density D) [Selivanov76]
//
static double s_rgfV10Response[64] = {
    1.146000,  1.087939,  1.049477,  1.009567,  0.979406,  0.954016,
    0.933022,  0.914348,  0.900804,  0.888205,  0.875469,  0.863894,
    0.854154,  0.845258,  0.837549,  0.829635,  0.822588,  0.815614,
    0.808080,  0.798682,  0.787597,  0.775944,  0.764487,  0.753425,
    0.741561,  0.729764,  0.717522,  0.706133,  0.694274,  0.682613,
    0.668656,  0.654463,  0.636824,  0.621254,  0.606496,  0.593464,
    0.576129,  0.562603,  0.546256,  0.528534,  0.511020,  0.493814,
    0.476500,  0.456686,  0.436933,  0.419022,  0.397030,  0.376685,
    0.356925,  0.337931,  0.321265,  0.304176,  0.284768,  0.263936,
    0.243365,  0.221782,  0.201184,  0.181208,  0.152054,  0.126196,
    0.099266,  0.071775,  0.042124,  0.000000,
};
static double s_rgfV09Response[64] = {
    1.440000,  1.384433,  1.342176,  1.299087,  1.258459,  1.214433,
    1.178239,  1.145326,  1.112713,  1.078642,  1.054863,  1.031528,
    1.009859,  0.994823,  0.973186,  0.955512,  0.935826,  0.917643,
    0.896935,  0.877774,  0.858883,  0.840457,  0.824234,  0.807710,
    0.787541,  0.767834,  0.754582,  0.739670,  0.722197,  0.704667,
    0.688370,  0.665924,  0.647213,  0.628312,  0.610143,  0.590143,
    0.574350,  0.560762,  0.543769,  0.526771,  0.509195,  0.492149,
    0.474389,  0.456203,  0.436332,  0.418267,  0.399684,  0.376447,
    0.356859,  0.338345,  0.322471,  0.304605,  0.290557,  0.264504,
    0.244820,  0.223420,  0.202815,  0.180124,  0.153250,  0.127796,
    0.100881,  0.074146,  0.041678,  0.000000,
};
//
//  Read 6-bit data (high order of 16-bit integers, like the V-13 and V-14 raw data)
//  Set bad-pixel mask (noisy data was hand-painted to zero)
//
static int
LoadImages(unsigned short rgnRaw[512][128], char rgnMask[512][128], char *szFile)
{
    HANDLE hFile;
    DWORD nBytesRead, nBytesToRead;
    int i, j;

    hFile = CreateFile(szFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                       FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("Cannot open %s\n", szFile);
        return 0;
    }
    for (j = 0; j < 512; j++)
        for (i = 0; i < 128; i++) {
            rgnRaw[j][i] = 0;
            rgnMask[j][i] = 0;
        }
    nBytesToRead = 2*128*512;
    ReadFile(hFile, rgnRaw, nBytesToRead, &nBytesRead, NULL);
    CloseHandle(hFile);
    //
    //  Zero values mark noise in the manually-masked versions
    //
    for (j = 0; j < 512; j++) {
        for (i = 0; i < 128; i++) {
            rgnRaw[j][i] &= 0xFC00;
            if (rgnRaw[j][i] == 0) {
                rgnMask[j][i] = BADPIXEL;   // was char(255)
            }
        }
    }
    //
    //  Make a rotated blow up of raw image by pixel replication
    //
    for (j = 0; j < 256; j++)
        for (i = 0; i < 1024; i++)
            s_rgnRaw2x[j][i] = rgnRaw[511 - i/2][j/2];
    return nBytesRead / 256;
}
//
//  Apply the radiometric response (linearization) reported by the Russians.  Convert
//  signal U to optical density D = log10(T), then convert that to transparancy T.
//  Expand contrast by bringing down the minimum level, but leave room for aperture
//  correction, which will sharpen the image, to dip down further.
//
static void
RadiometricResponse()
{
    int i, j, nU;
    double fD, fT, fMinT;

    printf("  B. Correct for radiometric response function\n");

    fMinT = pow(10.0, -s_rgfV09Response[0]);
    printf("    1. Venera 9: MinT = %f\n", fMinT);
    for (j = 0; j < 512; j++) {
        for (i = 0; i < 128; i++) {
            nU = s_rgnV09[j][i] >> 10;
            fD = s_rgfV09Response[nU];
            fT = pow(10.0, -fD);
            s_rgfV09[j][i] = float(fT);
        }
    }

    fMinT = pow(10.0, -s_rgfV10Response[0]);
    printf("    2. Venera 10: MinT = %f\n", fMinT);
    for (j = 0; j < 512; j++) {
        for (i = 0; i < 128; i++) {
            nU = s_rgnV10[j][i] >> 10;
            fD = s_rgfV10Response[nU];
            fT = pow(10.0, -fD);
            s_rgfV10[j][i] = float(fT);
        }
    }
}

//
//  Isophote interpolation of single missing pixels or one-pixel-wide columns
//
extern int IsoPhote(float *pfImage, char *pnMask, int ij, int nWide);
extern int IsoPhote2(float *pfImage, char *pnMask, int ij, int nWide);

static void
IsophoteRepair(float rgfImage[512][128], char rgnMask[512][128], int bDiffuse)
{
    int i, j, ij, nRepairs, nRepairs2, nDiffuse;

    nRepairs = nRepairs2 = nDiffuse = 0;
    for (j = 1; j < 511; j++) {
        //
        //  repair bad boarder pixels with linear interpolation
        //
        if (rgnMask[j][0] && !rgnMask[j-1][0] && !rgnMask[j+1][0]) {
            rgfImage[j][0] = Lerp(rgfImage[j-1][0], rgfImage[j+1][0], 0.5);
            rgnMask[j][0] = 0;
        }
        if (rgnMask[j][127] && !rgnMask[j-1][127] && !rgnMask[j+1][127]) {
            rgfImage[j][127] = Lerp(rgfImage[j-1][127], rgfImage[j+1][127], 0.5);
            rgnMask[j][127] = 0;
        }
        //
        //  Remaining pixels fixed by interpolation along isophotes
        //
        for (i = 127-1; i >= 1; --i) {
            if (rgnMask[j][i]) {
                ij = i + j*128;
                nRepairs += IsoPhote(rgfImage[0], rgnMask[0], ij, 128);
            }
        }
    }
    //
    //  A 45-degree version of IsoPhote, fixes a few more pixels.
    //
    for (j = 1; j < 511; j++) {
        for (i = 127-1; i >= 1; --i) {
            if (rgnMask[j][i]) {
                ij = i + j*128;
                nRepairs2 += IsoPhote2(rgfImage[0], rgnMask[0], ij, 128);
            }
        }
    }
    //
    //  Diffusion extrapolation from boundaries
    //
    if (bDiffuse) {
        for (j = 1; j < 511; j++) {
            for (i = 127-1; i >= 1; --i) {
                if (rgnMask[j][i]) {
                    ij = i + j*128;
                    nDiffuse += DiffuseExtrapolate(rgfImage[0], rgnMask[0], ij, 128);
                }
            }
        }
    }
    //REVIEW: set filled pixels mask to 2
    printf("    Isophote repair of %d and %d and %d pixels\n", nRepairs, nRepairs2, nDiffuse);
}
//
//  In-Paint missing regions of panorama with Bertalmio's algorithm.
//  Doesn't really work well, but it prevents aperture correction from ringing
//  at border of missing data.
//

static void
CloseBoundary(float rgfImage[512][128], char rgnMask[512][128], int iTop)
{
    int i, i2, j, jDelta;

    //
    //  extend data by mirroring it
    //
    for (j = 0; j < 512; j++) {
        for (i = 0; i < 256; i++) {
            if (i < 64+iTop) {
                s_rgfTemp[j][i] = rgfImage[j][64 + 2*iTop - i];
                s_rgnTemp[j][i] =  rgnMask[j][64 + 2*iTop - i];
            } else if (i > 64+127) {
                s_rgfTemp[j][i] = rgfImage[j][64 + 2*127 - i];
                s_rgnTemp[j][i] =  rgnMask[j][64 + 2*127 - i];
            } else {
                s_rgfTemp[j][i] = rgfImage[j][i - 64];
                s_rgnTemp[j][i] =  rgnMask[j][i - 64];
            }
        }
    }
    //
    //  Close boundary by pinching the ends shut
    //
    for (i2 = -6; i2 < 7; i2++) {
        i = i2 & 255;
        for (j = 0; j < 511; j++)
            if (s_rgnTemp[j][i])
                s_rgfTemp[j][i] = s_rgfTemp[j+1][i];
        for (j = 0; j < 511; j++)
            if (s_rgnTemp[j][i])
                s_rgnTemp[j][i] = s_rgnTemp[j+1][i];
        for (j = 511; j >= 1; --j)
            if (s_rgnTemp[j][i])
                s_rgfTemp[j][i] = s_rgfTemp[j-1][i];
        for (j = 511; j >= 1; --j)
            if (s_rgnTemp[j][i])
                s_rgnTemp[j][i] = s_rgnTemp[j-1][i];
    }
    for (i2 = -5; i2 < 6; i2++) {
        i = i2 & 255;
        for (j = 0; j < 511; j++)
            if (s_rgnTemp[j][i])
                s_rgfTemp[j][i] = s_rgfTemp[j+1][i];
        for (j = 0; j < 511; j++)
            if (s_rgnTemp[j][i])
                s_rgnTemp[j][i] = s_rgnTemp[j+1][i];
        for (j = 511; j >= 1; --j)
            if (s_rgnTemp[j][i])
                s_rgfTemp[j][i] = s_rgfTemp[j-1][i];
        for (j = 511; j >= 1; --j)
            if (s_rgnTemp[j][i])
                s_rgnTemp[j][i] = s_rgnTemp[j-1][i];
    }
    for (i2 = -4; i2 < 5; i2++) {
        i = i2 & 255;
        for (j = 0; j < 511; j++)
            if (s_rgnTemp[j][i])
                s_rgfTemp[j][i] = s_rgfTemp[j+1][i];
        for (j = 0; j < 511; j++)
            if (s_rgnTemp[j][i])
                s_rgnTemp[j][i] = s_rgnTemp[j+1][i];
        for (j = 511; j >= 1; --j)
            if (s_rgnTemp[j][i])
                s_rgfTemp[j][i] = s_rgfTemp[j-1][i];
        for (j = 511; j >= 1; --j)
            if (s_rgnTemp[j][i])
                s_rgnTemp[j][i] = s_rgnTemp[j-1][i];
    }
    //
    //  Shear the upper portion to match the 45-degree slope of the perspective
    //
    for (i = 64 + iTop - 1; i >= 0; --i) {
        jDelta = 2*(64 + iTop - i);
        for (j = 340; j < 512 - jDelta; j++) {
            s_rgfTemp[j][i] = s_rgfTemp[j + jDelta][i];
            s_rgnTemp[j][i] = s_rgnTemp[j + jDelta][i];
        }
        for (j = 170; j >= jDelta; --j) {
            s_rgfTemp[j][i] = s_rgfTemp[j - jDelta][i];
            s_rgnTemp[j][i] = s_rgnTemp[j - jDelta][i];
        }
    }
}

static void
InPaintNoise()
{
    int i, j;

    printf("  D. In-Paint regions of missing image data\n");
        printf("    1. Venera 9\n");
        CloseBoundary(s_rgfV09, s_rgnV09Mask, 13);
        InPainting(s_rgfTemp[0], s_rgnTemp[0], 512, 256, 0.90f);
        for (j = 0; j < s_nLinesV09; j++)
            for (i = 13; i < 128; i++)
                if (s_rgnV09Mask[j][i])
                    s_rgfV09[j][i] = s_rgfTemp[j][i + 64];
        printf("    2. Venera 10\n");
        CloseBoundary(s_rgfV10, s_rgnV10Mask, 15);
        InPainting(s_rgfTemp[0], s_rgnTemp[0], 512, 256, 0.60f);
        for (j = 0; j < s_nLinesV10; j++)
            for (i = 15; i < 128; i++)
                if (s_rgnV10Mask[j][i])
                    s_rgfV10[j][i] = s_rgfTemp[j][i + 64];
}

//
//  Perform aperture correction
//
static void
Contrast(float rgf[], char rgnMask[], int n)
{
    float fMin;
    int i;

    fMin = 1.0;
    for (i = 0; i < n; i++)
        if (rgnMask[i] == 0 && rgf[i] < fMin)
            fMin = rgf[i];
    for (i = 0; i < n; i++)
        rgf[i] = (rgf[i] - fMin)/(1.0f - fMin);
}

static void
ApertureCorrection(float rgfImage[512][128], char rgnMask[512][128])
{
    int i, j;

    Contrast(rgfImage[0], rgnMask[0], 512*128);
    Laplacian(s_rgfLaplace[0], 0.00f, rgfImage[0], 512, 128);   // don't do aperture here
    for (j = 0; j < 512; j++) {
        for (i = 0; i < 128; i++)
            s_rgfTemp[j][i] = 0.0;
        for (; i < 256; i++)
            s_rgfTemp[j][i] = s_rgfLaplace[j][i - 128];
    }
    Interpolate(s_rgfTemp[0], s_rgfWorkSpace, 512, 128, 0.15f); // do aperture here
    for (j = 0; j < 256; j++) {
        for (i = 0; i < 512; i++)
            s_rgfTranspose[j][i] = 0.0;
        for (; i < 1024; i++)
            s_rgfTranspose[j][i] = s_rgfTemp[1023 - i][j];
    }
    Interpolate(s_rgfTranspose[0], s_rgfWorkSpace, 256, 512);
}

static void
MaskTelemetry(float rgfImage[512][128], char rgnMask[512][128])
{
    int i, j;

    for (j = 0; j < 256; j++)
        for (i = 0; i < 1024; i++)
            if (rgnMask[511 - i/2][j/2])
                s_rgfTranspose[j][i] = 0.0;
}
//
//  Experimental attempt to correct bit-stream syncronization errors
//
static void
BitSyncImage(unsigned short rgnData[512][128])
{
    unsigned short *pn;
    int i, j, nShift, jWide, jFirst, jLimit;

    jWide = 512;
    jFirst = 0;
    jLimit = 512;
    pn = new unsigned short[128*7*jWide];
    for (nShift = -3; nShift < 4; nShift++) {
        for (j = 0; j < 256; j++) {
            for (i = 0; i < 128; i++)
                s_rgnTest[j][i] = rgnData[jLimit - j - 1 + jFirst][i];
        }
        for (j = 256; j < 512; j++) {
            for (i = 0; i < 128; i++)
                s_rgnTest[j][i] = rgnData[jLimit - j - 1 + jFirst+256][127-i];
        }
        BitShiftVenera9Pixels(s_rgnTest, nShift, jFirst, jLimit, 0, 252);
        for (j = 0; j < 512; j++) {
            for (i = 0; i < 128; i++)
                pn[j - jFirst + jWide*(i + 128*(nShift+3))] = s_rgnTest[j][i];
        }
    }
    WriteShortImage(pn, "V9BitSync.bmp", 128*7, jWide);
    delete [] pn;
}

//
//  Process the single smaller panoramas from Venera 9 and Venera 10
//
//  1. Load 6-bit data, zero values are "no data" indicators.
//  2. Convert nonlinear signal values into linear radiance
//  3. Repair isolated bad pixels or bad columns with isophote interpolation
//  4. Repair larger missing regions with in-painting
//  5. Apply aperture correction
//  6. Double size with windowed sinc
//  7. Write gamma-corrected final version
//

void
ProcessVenera9and10()
{

    printf("I. Process Venera-9 and Venera-10 Panoramas\n");
    LoadImages(s_rgnV09, s_rgnV09Mask, "B-09-Raw.raw");
    BitSyncImage(s_rgnV09);
return;
    WriteShortImage(s_rgnRaw2x[0], "..\\ProjectResults\\Venera9RawData.bmp", 256, 1024);
    s_nLinesV09 = LoadImages(s_rgnV09, s_rgnV09Mask, "B-09-Masked.raw");
    LoadImages(s_rgnV10, s_rgnV10Mask, "B-10-Raw.raw");
    WriteShortImage(s_rgnRaw2x[0], "..\\ProjectResults\\Venera10RawData.bmp", 256, 1024);
    s_nLinesV10 = LoadImages(s_rgnV10, s_rgnV10Mask, "B-10-Masked.raw");
    printf("  A. %3d lines in V-9, %3d lines in V-10\n", s_nLinesV09, s_nLinesV10);
    RadiometricResponse();
    printf("  C. Isophote repair of isolated pixels\n");
    IsophoteRepair(s_rgfV09, s_rgnV09Mask, 0);          // set bDiffuse=1 if not calling InPaintNoise
    IsophoteRepair(s_rgfV10, s_rgnV10Mask, 0);
    InPaintNoise();
    printf("  E. Aperture correction and time-base correction.\n");
    ApertureCorrection(s_rgfV09, s_rgnV09Mask);
    WriteFloatImage(s_rgfTranspose[0], "..\\ProjectResults\\Venera9Inpainted.bmp", 256, 1024);
    MaskTelemetry(s_rgfV09, s_rgnV09Mask);
    WriteFloatImage(s_rgfTranspose[0], "..\\ProjectResults\\Venera9.bmp", 256, 1024);
    ApertureCorrection(s_rgfV10, s_rgnV10Mask);
    WriteFloatImage(s_rgfTranspose[0], "..\\ProjectResults\\Venera10Inpainted.bmp", 256, 1024);
    MaskTelemetry(s_rgfV10, s_rgnV10Mask);
    WriteFloatImage(s_rgfTranspose[0], "..\\ProjectResults\\Venera10.bmp", 256, 1024);
}
