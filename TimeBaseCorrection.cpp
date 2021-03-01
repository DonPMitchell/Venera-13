#include "stdafx.h"
//
//  Correct video syncronization errors in Venera telemetry
//  D. P. Mitchell  05/15/2003.
//
#include "Venus.h"
#include "ImageProcessing.h"
#include "ImageFile.h"
//
//  Video Time-Base Correction
//
//  1. Find mis-synced lines by maximum-correlation with 41-bit video front porch.
//  2. Classify:    Early   - top equals bottom of previous line (buffer output early)
//                  Late    - bottom equals top of next line
//                  Both    - data matching previous and next lines
//                  Wrong   - no matches with previous or next line
//  3. Match patterns and correct:  NEW*N   - shift data helically backward
//                                  NW*LN   - shift data helically forward
//                                  NBN     - delete line B (never been seen)
//                                  NW*N    - insert line after W and helically shift forward
//REVIEW: NW*N forward or backward? check rotation sign
#define SYNC_NORMAL '-'
#define SYNC_EARLY  'E'
#define SYNC_LATE   'L'
#define SYNC_BOTH   'B'
#define SYNC_WRONG  'W'

#define TBC_FULLREPAIR  0
#define TBC_NOSHIFTING  1
#define TBC_ROTATEONLY  2

#define NEGATE(N)   ((N) ^ 0xFF00)
#define UNNEGATE(N) (N)

static char             s_rgnSync[MAXLINES];
static unsigned char    s_rgnSyncShift[MAXLINES];
static unsigned char    s_rgnNegated[MAXLINES];

//REVIEW: build new version in a second buffer, will lose fewer bits

void
TimeBaseCorrection(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
                   int jRamp, int nType, double fCorrLimit, int jAuxRamp)
{
    unsigned short rgnLine[252+41], rgnRamp[41];
    unsigned short rgnAuxRamp[41], rgnRampC[41], rgnAuxRampC[41];
    int i, j, iMaxCorr, k;
    int i2, nPrevMatch, nPrevDiffer, nNextMatch, nNextDiffer, nZero;
    double fCorr, fMaxCorr, fEarly, fLate;

    printf("    Time-base correction %4d to %4d\n      ", jFirst, jLimit);
    for (j = jFirst; j < jLimit; j++)
        for (i = 0; i < 252; i++)
            g_rgnTest2[j][i] = BADPIXEL;            // assemble rotated lines in buffer
    for (j = 0; j < MAXLINES; j++) {
        s_rgnSync[j] = 0;
        s_rgnSyncShift[j] = 0;
    }
    if (jRamp < 0)
        jRamp = jFirst;
    for (k = 0; k < 4; k++)
        if (!(g_rgnRaw[k][jRamp][20] & BADPIXEL))  // rgnData not defined at jRamp
            break;
    if (k == 4) {
        printf("No calibration ramp at %d\n", jRamp);
        return;
    }
    for (i = 0; i < 41; i++) {
        rgnRamp[i] = g_rgnRaw[k][jRamp][i];
    }
    if (jAuxRamp > 0) {
        for (i = 0; i < 41; i++) {
            rgnAuxRamp[i] = g_rgnRaw[k][jAuxRamp][i];
            rgnRampC[i] = (~g_rgnRaw[k][jRamp][i]) & NINEBITS;
            rgnAuxRampC[i] = (~g_rgnRaw[k][jAuxRamp][i]) & NINEBITS;
        }
    }
    for (j = jFirst; j < jLimit; j++) {
        for (i = 0; i < 252; i++)
            rgnLine[i] = rgnData[j][i];
        for (; i < 252+41; i++)
            rgnLine[i] = rgnData[j][i-252];
        iMaxCorr = 0;
        fMaxCorr = -1.0;
        for (i = 0; i < 252; i++) {
            fCorr = Correlation(rgnRamp, rgnLine+i, 41);
            if (fCorr > fMaxCorr) {
                iMaxCorr = i;
                fMaxCorr = fCorr;
            }
            if (jAuxRamp > 0) {
                fCorr = Correlation(rgnAuxRamp, rgnLine+i, 41);
                if (fCorr > fMaxCorr) {
                    iMaxCorr = i;
                    fMaxCorr = fCorr;
                }
                /*
                fCorr = Correlation(rgnRampC, rgnLine+i, 41);
                if (fCorr > fMaxCorr) {
                    iMaxCorr = i;
                    fMaxCorr = fCorr;
                }
                fCorr = Correlation(rgnAuxRampC, rgnLine+i, 41);
                if (fCorr > fMaxCorr) {
                    iMaxCorr = i;
                    fMaxCorr = fCorr;
                }
                */
            }
        }
        if (iMaxCorr != 0 && fMaxCorr > fCorrLimit) {
            nPrevMatch = nPrevDiffer = nNextMatch = nNextDiffer = 0;
            fEarly = fLate = 0.0;
            nZero = rgnData[j][iMaxCorr] >> 7;   // supposed to be zero
            //
            //  If output is early, some data from previous line is still in buffer.
            //  If late, maybe some data will appear again in the next line.
            //
            for (i = 0; i < iMaxCorr; i++) {
                i2 = 252 - iMaxCorr + i;
                if (IsGood(rgnData[j-1][i2])) {
                    if (rgnData[j-1][i2] == rgnData[j][i])
                        nPrevMatch++;
                    else
                        nPrevDiffer++;
                }
            }
            for (; i < 252; i++) {
                i2 = i - iMaxCorr;
                if (IsGood(rgnData[j+1][i2])) {
                    if (rgnData[j+1][i2] == rgnData[j][i])
                        nNextMatch++;
                    else
                        nNextDiffer++;
                }
            }
            s_rgnSyncShift[j] = iMaxCorr;
            fEarly = double(nPrevMatch)/double(iMaxCorr);
            fLate = double(nNextMatch)/double(252-iMaxCorr-7);
            if (fEarly > 0.5 && fLate > 0.5) {
                s_rgnSync[j] = SYNC_BOTH;      // both late and early
            } else if (fEarly > 0.5) {
                s_rgnSync[j] = SYNC_EARLY;
            } else if (fLate > 0.5) {
                s_rgnSync[j] = SYNC_LATE;
            } else {
                s_rgnSync[j] = SYNC_WRONG;     // neither late nor early, just wrong
            }
            if (s_rgnSync[j] == SYNC_WRONG && s_rgnSync[j-1] == SYNC_NORMAL &&
                            (iMaxCorr == 1 || iMaxCorr == 251 || fMaxCorr < 0.6)) {
                s_rgnSync[j] = SYNC_NORMAL;     // just caused by noise, not video sync
                s_rgnSyncShift[j] = 0;
            }
            if (s_rgnSync[j] == SYNC_WRONG && nType == TBC_NOSHIFTING) {
                if (s_rgnSyncShift[j] < 126)
                    s_rgnSync[j] = SYNC_EARLY;
                else
                    s_rgnSync[j] = SYNC_LATE;
            }
        } else {
            s_rgnSync[j] = SYNC_NORMAL;
        }
        printf("%c", s_rgnSync[j]);
        if (j%70 == 69)
            printf("\n      ");
    }
    printf("\n");
    //
    //  This mode is used for study of noisy regions like last red 14-I panorama
    //
    if (nType == TBC_ROTATEONLY) {
        for (j = jFirst; j < jLimit; j++) {
            if (s_rgnSync[j] != SYNC_NORMAL)
                RotateColumns(rgnData, j, j+1, s_rgnSyncShift[j], ROT_ROTATE);
            if (s_rgnNegated[j])
                for (i = 0; i < 252; i++)
                    rgnData[j][i] = UNNEGATE(rgnData[j][i]);
        }
        return;
    }
    //
    //  Scan forward and repair EW* errors.  nRots are usually equal for EW+ patterns.
    //
    for (j = jFirst; j < jLimit; ) {
        if (s_rgnSync[j] == SYNC_EARLY) {
            do {
                RotateColumns(rgnData, j, j+1, s_rgnSyncShift[j], ROT_HELICAL);
                s_rgnSync[j] = SYNC_NORMAL;
                j++;
            } while (s_rgnSync[j] == SYNC_WRONG);
        } else
            j++;
    }
    //
    //  Scan backward and repair W*L errors.
    //
    for (j = jLimit-1; j >= jFirst; ) {
        if (s_rgnSync[j] == SYNC_LATE) {
            do {
                RotateColumns(rgnData, j, j+1, s_rgnSyncShift[j] - 252, ROT_HELICAL);
                s_rgnSync[j] = SYNC_NORMAL;
                --j;
            } while (s_rgnSync[j] == SYNC_WRONG);
        } else
            --j;
    }
    //
    //  Scan backward and repair N*WN* errors.  These are believed to be caused
    //  by very very late trigger of buffer output, causing a scanline to be skipped.
    //  SYNC_BOTH events not seen.  If they are, column is repeat, delete it.
    //
    for (j = jLimit-1; j >= jFirst; ) {
        if (s_rgnSync[j] == SYNC_WRONG) {
            if (nType == TBC_FULLREPAIR)
                ShiftColumns(rgnData, j+1, 1);
            printf("      NW*N event at %4d\n", j+1);
            do {
                RotateColumns(rgnData, j, j+1, s_rgnSyncShift[j] - 252, ROT_HELICAL);
                s_rgnSync[j] = SYNC_NORMAL;
                --j;
            } while (s_rgnSync[j] == SYNC_WRONG);
        } else
            --j;
    }
}
//
//  An ad hoc scheme for repairing the last red panorama in 14-I
//
void
TimeBaseRecode(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit, int jRamp)
{
    unsigned short rgnLine[252+41], rgnRamp[41], rgnNegRamp[41];
    int i, j, iMaxCorr, bNegate, k, jLastRot, nLastRot;
    double fCorr, fMaxCorr;

    for (j = 0; j < MAXLINES; j++)
        s_rgnNegated[j] = 0;
    if (jRamp < 0)
        jRamp = jFirst;
    for (k = 0; k < 4; k++)
        if (!(g_rgnRaw[k][jRamp][20] & BADPIXEL))  // rgnData might not be defined at jRamp
            break;
    for (i = 0; i < 41; i++)
        rgnRamp[i] = g_rgnRaw[k][jRamp][i];
    for (i = 0; i < 41; i++)
        rgnNegRamp[i] = NEGATE(rgnRamp[i]);     // to fix V-14-I last red panorama
    jLastRot = nLastRot = 0;
    for (j = jFirst; j < jLimit; j++) {
        for (i = 0; i < 252; i++)
            rgnLine[i] = rgnData[j][i];
        for (; i < 252+41; i++)
            rgnLine[i] = rgnData[j][i-252];
        iMaxCorr = 0;
        fMaxCorr = -1.0;
        bNegate = 0;
        for (i = 0; i < 252; i++) {
            fCorr = Correlation(rgnRamp, rgnLine+i, 41);
            if (fCorr > fMaxCorr) {
                iMaxCorr = i;
                fMaxCorr = fCorr;
                bNegate = 0;
            }
            fCorr = Correlation(rgnNegRamp, rgnLine+i+2, 41);
            if (fCorr > fMaxCorr) {
                iMaxCorr = i;
                fMaxCorr = fCorr;
                bNegate = 1;
                s_rgnNegated[j] = 1;
            }
        }
        if (bNegate)
            for (i = 0; i < 252; i++)
                rgnData[j][i] =  NEGATE(rgnData[j][i]);
    }
}
//
//  New Time Base Correction.
//
//  The problem is better understood now.  Data is copied into a buffer, which
//  saves more pixels than helical rotations in place.  The old way of measuring
//  "early" and "late" buffer dumps did not handle the relative shifting of very
//  badly damaged data.
//
#define MAXTBC  2000

static unsigned short s_rgnBuffer[MAXTBC][252];
static unsigned short s_rgnAlternate[MAXTBC][252];
//
//  Match with the 41-line calibration front porch, including the photometric ramp
//
static double
MatchCalibration(unsigned short *pnRamp, unsigned short *pnData)
{
    int i, nPorch;
    double fCorr;

    nPorch = 0;
    for (i = 0; i < 5; i++)
        nPorch += (pnData[i] == g_rgnPorch[i]);
    fCorr = Correlation(pnRamp, pnData, 41);
    if (nPorch >= 3 && fCorr > 0.5)             // bonus for digital porch
        return double(nPorch)/20.0 + fCorr;
    else
        return fCorr;
}
//
//  Shifted data above calibration should equal some of the previous line, if
//  previous pixels are there.  If not, a whole line might be missing, which is rare.
//
static int
Similarity(unsigned short *pnData, unsigned short *pnPrev, int n,
           int *pnGood = 0, int *pnEqual = 0)
{
    int i, nGood, nEqual;

    nGood = nEqual = 0;
    for (i = 0; i < n; i++) {
        if (IsGood(pnPrev[i])) {
            nGood++;
            nEqual += (pnPrev[i] == pnData[i]);
        }
    }
    if (pnGood)
        *pnGood = nGood;
    if (pnEqual)
        *pnEqual = nEqual;
    return n - nGood + nEqual;
}

void
NewTimeBaseCorrection(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
                      int jRamp)
{
    int i, j, jBuf, jLate, iMaxCorr, iLastMax, nSim, nGood, nEqual, nShift;
    unsigned short rgnRamp[41], rgnLine[252+41];
    double fCorr, fMaxCorr, rgfMaxCorr[MAXTBC];
    int rgiMaxCorr[MAXTBC];

    for (i = 0; i < 41; i++) {
        rgnRamp[i] = rgnData[jRamp][i];
        rgnLine[252+i] = BADPIXEL;
    }
    for (j = 0; j < MAXTBC; j++) {
        for (i = 0; i < 252; i++) {
            s_rgnBuffer[j][i] = BADPIXEL;
            s_rgnAlternate[j][i] = BADPIXEL;
        }
        rgiMaxCorr[j] = 0;
        rgfMaxCorr[j] = -1.0;
    }
    iLastMax = 0;
    for (j = jFirst; j < jLimit; j++) {
        iMaxCorr = 0;
        fMaxCorr = -1.0;
        for (i = 0; i < 252; i++)
            rgnLine[i] = rgnData[j][i];
        for (i = 0; i < 252-20; i++) {
            fCorr = MatchCalibration(rgnRamp, rgnLine+i);
            if (fCorr > fMaxCorr) {
                iMaxCorr = i;
                fMaxCorr = fCorr;
            }
        }
        rgiMaxCorr[j-jFirst+1] = iMaxCorr;
        rgfMaxCorr[j-jFirst+1] = fMaxCorr;
    }
    jBuf = 1;
    for (j = jFirst; j < jLimit; j++) {
        iMaxCorr = rgiMaxCorr[j-jFirst+1];
        fMaxCorr = rgfMaxCorr[j-jFirst+1];
        if (fMaxCorr > 0.33 && iMaxCorr == 0) {
            //
            //  Normal scanline, front porch at top
            //
            nSim = Similarity(&rgnData[j][0], &s_rgnBuffer[jBuf][0], 252);
            if (nSim < 252/2) //REVIEW: do a proper confidence limit here
                jBuf++;
            for (i = 0; i < 252; i++) {
                if (IsBad(s_rgnBuffer[jBuf][i]))
                    s_rgnBuffer[jBuf][i] = rgnData[j][i];
                else if (s_rgnBuffer[jBuf][i] != rgnData[j][i])
                    s_rgnAlternate[jBuf][i] = rgnData[j][i];
            }
            jBuf++;
        } else {
            //  Out-of-sync line
            //
            //      j: __________________ABCDEFGH
            //    j+1: IJKLMNOPQRSTUVWXYZ________
            //
            //  If its early, ABCDEFGH should match previous line.  If late,
            //  move it over a line, and IJKL... should match the next line.
            //
            if (fMaxCorr < 0.33)
                iMaxCorr = iLastMax;    // junky line
            nSim = Similarity(&rgnData[j][0], &s_rgnBuffer[jBuf-1][252-iMaxCorr],
                              iMaxCorr, &nGood, &nEqual);
            jLate = (nSim < iMaxCorr/2);
            for (i = 0; i < iMaxCorr; i++) {
                if (IsBad(s_rgnBuffer[jBuf-1+jLate][252 - iMaxCorr + i]))
                    s_rgnBuffer[jBuf-1+jLate][252 - iMaxCorr + i] = rgnData[j][i];
                else if (s_rgnBuffer[jBuf-1+jLate][252 - iMaxCorr + i] != rgnData[j][i])
                    s_rgnAlternate[jBuf-1+jLate][252 - iMaxCorr + i] = rgnData[j][i];
            }
            for (i = iMaxCorr; i < 252; i++)
                s_rgnBuffer[jBuf+jLate][i - iMaxCorr] = rgnData[j][i];
            printf("%4d %3d %8.4f %8.4f %d %d\n", j, iMaxCorr, fMaxCorr,
                    double(nSim)/double(iMaxCorr), nGood, nEqual);
            jBuf += 1;
        }
        iLastMax = iMaxCorr;
    }
    printf("TBC shift: %4d vs. %4d\n", jFirst+jBuf, jLimit);
    WriteShortImage(s_rgnBuffer[0], "tbc.bmp", jBuf, 252);
    WriteShortImage(s_rgnAlternate[0], "tbcalt.bmp", jBuf, 252);
    nShift = jFirst + jBuf - jLimit - 1;    // jBuf started out +1
    Assert(nShift >= 0);
    if (nShift)
        ShiftColumns(rgnData, jLimit, nShift);
    for (j = jFirst; j < jFirst+jBuf-1; j++)
        for (i = 0; i < 252; i++)
            rgnData[j][i] = s_rgnBuffer[j-jFirst+1][i];
}




