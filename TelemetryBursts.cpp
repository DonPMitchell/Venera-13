#include "stdafx.h"
//
//  Find and flag periodic telemetry bursts
//  D. P. Mitchell  05/15/2003.
//
#include "Venus.h"
#include "ImageProcessing.h"

//
//  Mask out the periodic telemetry bursts, automatically and manually.
//  Carefully salvage partial scan line in first and last line of a burst.
//
#define SHORT_BURST 244
#define LONG_BURST  477

void
MaskBurst(int jFirst, int jLimit, int iFirst, int iLimit, int nRotate)
{
    int i, j, k, iTM, jTM;

    iTM = 0;
    jTM = g_jTelemetry;
    for (k = 0; k < 4; k++) {
        for (j = jFirst; j < jLimit; j++) {
            if (j == jFirst) {
                for (i = iFirst; i < 252; i++) {
                    if (k == 0 && iTM < 4096)
                        g_rgnTelemetry[jTM][iTM++] = g_rgnRaw[k][j][i];
                    SetBurst(g_rgnRaw[k][j][i]);
                }
            } else if (j == jLimit - 1) {
                for (i = 0; i < iLimit; i++) {
                    if (k == 0 && iTM < 4096)
                        g_rgnTelemetry[jTM][iTM++] = g_rgnRaw[k][j][i];
                    SetBurst(g_rgnRaw[k][j][i]);
                }
            } else {
                for (i = 0; i < 252; i++) {
                    if (k == 0 && iTM < 4096)
                        g_rgnTelemetry[jTM][iTM++] = g_rgnRaw[k][j][i];
                    SetBurst(g_rgnRaw[k][j][i]);
                }
            }
        }
        if (nRotate)
            RotateColumns(g_rgnRaw[k], jLimit-1, jLimit, nRotate, ROT_ROTATE);
    }
    g_jTelemetry++;
}

inline static int
Sign(int n)
{
    if (n < 0)
        return -1;
    if (n > 0)
        return 1;
    return 0;
}
//
//  Find telemetry bursts.  Rotate the last column if necessary, to re-align
//  image data.
//
void
AutoMaskBurst(unsigned short rgnData[MAXLINES][252], int jFirstBurstCenter,
              int bShort, int jBurstLimit)
{
    int iFirst, iLimit, i, j, jBurst, jLastBurst, jFirst, jLimit, bShift, iMaxCorr;
    int iDelta, nRotate, nLastStart, nShort, nLong;
    double fCorr, fRemainCorr, fMaxCorr;

    nShort = SHORT_BURST;
    nLong = LONG_BURST;
    jLastBurst = 0;
    nLastStart = 0;
    printf("    AutoMaskBurst\n");
    for (jBurst = jFirstBurstCenter; jBurst < jBurstLimit; bShort = !bShort, jBurst += bShort ? nShort : nLong) {
        if (rgnData[jBurst][0] & BADPIXEL)
            continue;
        //
        //  Bursts are about 244 and 477 lines apart.  Search for exact bounds.
        //
        Assert(jBurst < jBurstLimit);
        for (j = jBurst - 8; j < jBurst; j++) {
            fCorr = Correlation(&rgnData[j][0], &rgnData[j-1][0], 42);
            if (fCorr < 0.90) {
                jFirst = j;
                break;
            }
        }
        for (j = jBurst + 8; j > jBurst; --j) {
            fCorr = Correlation(&rgnData[j][0], &rgnData[j+1][0], 42);
            if (fCorr < 0.90) {
                jLimit = j+1;
                break;
            }
        }
        //
        //  Beginning of burst is a zero pixel.
        //
        iFirst = 0;
        for (i = 1; i < 252; i++) {
            if (DATA(rgnData[jFirst - 1][i]) == 0) {
                iFirst = i;
                jFirst = jFirst -1;
                break;
            }
        }
        //
        //  Zero pixel before end of bursts also
        //
        iLimit = 252;
        bShift = 0;
        for (i = 251; i > 1; --i) {
            if (rgnData[jLimit-1][i] == 0) {
                iLimit = i + 5;
                break;
            }
        }
        nRotate = 0;
        if (iLimit > 252)
            iLimit = 252;
        if (iLimit < 252) {
            //
            //  Is the remainder of the line strongly correlated with neighbor?
            //
            fMaxCorr = -1.0;
            iMaxCorr = -1;
            fRemainCorr = Correlation(&rgnData[jLimit-1][iLimit], &rgnData[jLimit][iLimit], 252 - iLimit);
            //
            //  Does burst end with mis-positioned image data?  Find max correlation
            //  position.  Watch out for duplicated copy of pixels from next line.
            //
            for (iDelta = iLimit; iDelta > 1; --iDelta) {
                fCorr = Correlation(&rgnData[jLimit-1][iLimit], &rgnData[jLimit][iDelta], 252 - iLimit);
                if (fCorr > fMaxCorr) {
                    fMaxCorr = fCorr;
                    iMaxCorr = iDelta;
                }
            }
            if (fMaxCorr < 0.5 || fMaxCorr > 0.99999)
                iLimit = 252;
            else if (fRemainCorr < 0.5)
                nRotate = iLimit - iMaxCorr;
        }
        if (jLimit - jFirst <= 12) {
            printf("      Burst at %7d %6d %4d ", iFirst + 252*jFirst, iFirst+252*jFirst-nLastStart, jLimit-1);
            nLastStart = iFirst+252*jFirst;
            printf("Video Corr= %6.3f ", fRemainCorr);
            if (iMaxCorr > -1)
                printf("Max Corr= %6.3f ", fMaxCorr);
            if (nRotate)
                printf("Rot %3d", nRotate);
            if (iLimit == 252)
                printf("Kill");
            else if (!nRotate)
                printf("Keep");
            printf("\n");
            MaskBurst(jFirst, jLimit, iFirst, iLimit, nRotate);
        } else {
            printf("      BadBurst %7d\n", iFirst + 252*jFirst);
            continue;
        }
        jLastBurst = jBurst;
        jBurst = (jLimit + jFirst - 1)/2;           // correct the burst center
        if (bShort)
            nShort += Sign(jBurst - jLastBurst);    // correct inter-burst distances
        else
            nLong += Sign(jBurst - jLastBurst);
    }
}
