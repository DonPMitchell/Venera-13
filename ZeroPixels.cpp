#include "stdafx.h"
//
//  Remove zero pixels that are inserted into data stream
//  D. P. Mitchell  05/15/2003.
//
#include "Venus.h"
#include "ImageProcessing.h"
//
//  These procedures solves a special problem in V-13-II and other transmissions.
//  Zeros are occasionally inserted, shifting the rest of a scan line down.  Sometimes
//  a zero heralds more bizarre damage, like 62 missing lines in V-13-I.  In regions of
//  strong salt & pepper noise, zeros occur naturally, without shifting.
//
void
ShiftSubColumn(unsigned short rgnData[MAXLINES][252], int jFirst,
               int iFirst, int iLimit)
{
    int i;

    //
    //  iFirst is the pixel below the zero pixel.  That was dumb.
    //
    for (i = iFirst; i < iLimit; i++)
        rgnData[jFirst][i-1] = rgnData[jFirst][i];  // iFirst moved to iFirst-1
    SetBad(rgnData[jFirst][iLimit - 1]);
}

void
NegShiftSubColumn(unsigned short rgnData[MAXLINES][252], int jFirst,
               int iFirst, int iLimit)
{
    int i;

    //
    //  iFirst is the pixel below the zero pixel.  That was dumb.
    //
    for (i = iLimit-1; i > iFirst; --i)
        rgnData[jFirst][i] = rgnData[jFirst][i-1];  // iFirst moved to iFirst+1
    //SetBad(rgnData[jFirst][iFirst]);
    rgnData[jFirst][iFirst] = BADPIXEL|BURSTPIXEL;
}

static int
ScanningDirection(unsigned short rgnData[MAXLINES][252], int j)
{
    int i, nGoodUpper, nGoodLower;
    double fUpper, fLower;

    fUpper = fLower = 0.0;
    nGoodUpper = nGoodLower = 0;
    for (i = 11; i <= 21; i++) {
        fUpper += UVALUE(rgnData[j][i], 9);
        nGoodUpper += IsGood(rgnData[j][i]);
    }
    for (; i <= 31; i++) {
        fLower += UVALUE(rgnData[j][i], 9);
        nGoodLower += IsGood(rgnData[j][i]);
    }
    if (nGoodUpper < 5 || nGoodLower < 5 || fUpper == fLower)
        return 0;
    else if (fUpper > fLower)
        return 1;
    else
        return -1;
}

static double
SquareError(unsigned short rgnData[MAXLINES][252], int j, int iShift,
            int iFirst, int iLimit)
{
    double rgnEstimate[252];
    int i;
    double fError, fPixel, fEst, fSum, fN;

    for (i = 0; i < 252; i++) {
        if (IsGood(rgnData[j-1][i]) && IsGood(rgnData[j+1][i]))
            rgnEstimate[i] = 0.5*(UVALUE(rgnData[j-1][i], 9) + UVALUE(rgnData[j+1][i], 9));
        else if (IsGood(rgnData[j-1][i]))
            rgnEstimate[i] = UVALUE(rgnData[j-1][i], 9);
        else if (IsGood(rgnData[j+1][i]))
            rgnEstimate[i] = UVALUE(rgnData[j+1][i], 9);
        else
            rgnEstimate[i] = -1.0;
    }
    fSum = 0.0;
    fN = 1.0e-6;
    for (i = iFirst; i < iLimit; i++) {
        if (IsBad(rgnData[j][i]) || rgnEstimate[i] == -1.0)
            continue;
        fEst = rgnEstimate[i-iShift];
        fPixel = UVALUE(rgnData[j][i], 9);
        fError = (fEst - fPixel);
        fError = fError*fError;
        fSum += fError*fError;
        fN += 1.0;
    }
    return sqrt(fSum/fN);
}

#define VERBOSE (kVersion == g_kVerbose)

void
FixZeroPixels(unsigned short rgnFix[MAXLINES][252],
                   unsigned short rgnCheck[MAXLINES][252],
                   int kVersion)
{
    int i, j, jFirst, jLimit, iLook, iShift, nFixed, nZeros, bMult, nCount, nBrown, bCheck;
    double rgfSim[3], rgfSelfSim[3], rgfError[3];
    double fStat0, fStat1, fStat2, fStat3, fStat4;

    nFixed = nBrown = 0;
    for (jFirst = 2; ; jFirst++)
        if (IsGood(rgnFix[jFirst][11]) || IsGood(rgnFix[jFirst][53]))
            break;
    for (jLimit = MAXLINES-2; ; --jLimit)
        if (IsGood(rgnFix[jLimit-1][11]) || IsGood(rgnFix[jLimit-1][53]))
            break;
    for (j = jFirst; j < jLimit; j++) {
        //
        //  rgnData[j][0] is always zero, rgnData[j][251] can't be repaired.
        //
        if (IsBurst(rgnFix[j][53]))
            continue;
        if (ScanningDirection(rgnFix, j) == -ScanningDirection(rgnFix, j+2) ||
            ScanningDirection(rgnFix, j-1) == -ScanningDirection(rgnFix, j+1))
            continue;   // always some zeros when the camera reverses
        for (i = 250; i >= 1; --i) {
            if (ISZERO(rgnFix[j][i])) {
                //
                //  Is there data to compare with in the check version?
                //
                bCheck = 0;
                for (iLook = i; iLook < 252; iLook++) {
                    if (IsGood(rgnCheck[j][iLook])) {
                        bCheck = 1;
                        break;
                    }
                }
                //
                //  Human inspection if situation too complex for routine
                //
                if (ISZERO(rgnCheck[j][i]) && !IsBurst(rgnCheck[j][i])) {
                    if (VERBOSE)
                        printf("INSPECT %4d %3d, zeros in both versions\n", j, i);
                    goto NextColumn;
                }
                for (nZeros = 1; i - nZeros > 0 && ISZERO(rgnFix[j][i-nZeros]); nZeros++)
                    ;
                bMult = 0;
                for (iLook = 1; iLook < i - nZeros; iLook++) {
                    if (ISZERO(rgnFix[j][iLook])) {
                        bMult++;
                    }
                }
                //
                //  Statistics comparing possibly shifted segment with the other
                //  version (similitude) or with its neighborhood (rms error).
                //
                nCount = 252 - i - 1;
                for (iShift = 0; iShift <= (1 + (nZeros>1)); iShift++) {
                    rgfSim[iShift]     = Similitude(rgnFix[j]+i+1, rgnCheck[j]+i+1-iShift, nCount);
                    rgfSelfSim[iShift] = Similitude(rgnFix[j]+i+1+iShift, rgnFix[j]+i+1, nCount-iShift);
                    rgfError[iShift] = SquareError(rgnFix, j, iShift, i + 1, 252);
                }
                fStat0 = (fabs(rgfSim[1] - rgfSelfSim[1]) +
                               fabs(rgfSim[0] - 1.0)) * sqrt(double(nCount));
                fStat1 = (fabs(rgfSim[0] - rgfSelfSim[1]) +
                               fabs(rgfSim[1] - 1.0)) * sqrt(double(nCount));
                fStat2 = (fabs(rgfSim[0] - rgfSelfSim[2]) +
                               fabs(rgfSim[2] - 1.0)) * sqrt(double(nCount));
                fStat3 = rgfError[0]/(rgfError[0] + rgfError[1] + 1.0e-6);
                fStat4 = rgfError[2]/(rgfError[0] + rgfError[2] + 1.0e-6);
                if (nZeros == 1) {
                    //
                    //  The similitude two-sided hypothesis test, very reliable
                    //
                    if (bCheck && fStat0 > 2.5 && fStat1 < 0.8) {
                        if (bMult) {
                            if (VERBOSE)
                                printf("INSPECT %4d %3d, multiple zeros above\n", j, i);
                            goto NextColumn;
                        }
                        if (VERBOSE)
                            printf("shifting %4d %3d\n", j, i);
                        if (IsPorch(&rgnFix[j][i])) printf("Front Porch\n");
                        ShiftSubColumn(rgnFix, j, i+1);
                        if (g_rgnFlags[kVersion] & FLAGS_ANOMALOUS)
                            ShiftSubColumn(g_rgnBackup[kVersion], j, i+1);
                        nFixed++;
                        //
                        //  In the Brown U. data, zero pixels cause the other half
                        //  of the line to shift up, so we can replace the missing
                        //  pixel with an 8-bit version (if there is no 9-bit data
                        //  in the check version).
                        //
                        if (IsBad(rgnCheck[j][251]))
                            nBrown += BrownReplace(kVersion, j, 251);
                        goto NextColumn;
                    } else if (bCheck && fStat0 < 0.8 && fStat1 > 2.5) {
                        goto NextColumn;
                    } else if (fStat3 > g_fShiftLimit) {
                        //
                        //  With multiple versions, VisualizeDifference makes it easy
                        //  to spot any auto-shift mistakes.  In the last part of
                        //  Venera 13 data, we only have one version, so we print
                        //  a warning and visually inspect what was done.
                        //
                        if (bMult) {
                            if (VERBOSE)
                                printf("INSPECT %4d %3d, multiple zeros above\n", j, i);
                            goto NextColumn;
                        }
                        if (!bCheck)
                            if (VERBOSE)
                                printf("SHIFTED: %4d %3d, %f %f\n",
                                        j, i, fStat3, fStat4);
                        ShiftSubColumn(rgnFix, j, i+1);
                        if (g_rgnFlags[kVersion] & FLAGS_ANOMALOUS)
                            ShiftSubColumn(g_rgnBackup[kVersion], j, i+1);
                        nFixed++;
                        if (IsBad(rgnCheck[j][251]))
                            nBrown += BrownReplace(kVersion, j, 251);
                        goto NextColumn;
                    } else {
                        if (!bCheck && fStat3 > 0.45)
                            if (VERBOSE)
                                printf("INSPECT %4d %3d, no confident hypothesis %f\n",
                                            j, i, fStat3);
                        goto NextColumn;
                    }
                } else {
                    //
                    //  The case of 2 zeros, usually means a shift by 2.  Never see
                    //  3-zero cases that should be shifted.
                    //
                    if (bCheck && fStat0 > 2.5 && fStat1 > 2.5 && fStat2 < 0.8) {
                        if (nZeros > 2) {
                            if (VERBOSE)
                                printf("INSPECT %4d %3d, %2d zeros\n", j, i, nZeros);
                            goto NextColumn;
                        }
                        if (VERBOSE)
                            printf("Shift Twice %4d %3d\n", j, i);
                        ShiftSubColumn(rgnFix, j, i);   // not i+1 for double zeros
                        ShiftSubColumn(rgnFix, j, i);
                        if (g_rgnFlags[kVersion] & FLAGS_ANOMALOUS) {
                            ShiftSubColumn(g_rgnBackup[kVersion], j, i+1);
                            ShiftSubColumn(g_rgnBackup[kVersion], j, i+1);
                        }
                        if (IsBad(rgnCheck[j][250]))
                            nBrown += BrownReplace(kVersion, j, 250);
                        if (IsBad(rgnCheck[j][251]))
                            nBrown += BrownReplace(kVersion, j, 251);
                        nFixed++;
                        goto NextColumn;
                    } else if (bCheck && fStat0 < 0.8 && fStat1 > 2.5 && fStat2 > 2.5) {
                        goto NextColumn;
                    } else {
                        if (nZeros < 5)
                            if (VERBOSE)
                                printf("INSPECT %4d %3d, no hypothesis (%d zeros) %f %f\n",
                                            j, i, nZeros, fStat3, fStat4);
                        goto NextColumn;
                    }
                }
                goto NextColumn;
            }
        }
NextColumn:
        ;
    }
    printf("    AutoShiftSubColumn: %d lines fixed, %d brown\n", nFixed, nBrown);
}
