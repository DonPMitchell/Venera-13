#include "stdafx.h"
//
//  Assemble multiple transmissions into one high-quality master version
//  D. P. Mitchell  05/15/2003.
//
#include "Venus.h"
#include "Recoding.h"
#include "ImageProcessing.h"

#define NSEG    21

static void
VisualizeNoise2(unsigned short rgnError[MAXLINES][252],
               unsigned short rgn[MAXLINES][252], double fSigma)
{
    int i, j, nPixel, nSigma, rgnLow[32], rgnHigh[32], rgnSigmas[32];
    double fError, fPixel;

    for (i = 0; i < 32; i++)
        rgnLow[i] = rgnHigh[i] = rgnSigmas[i] = 0;
    for (j = 1; j < MAXLINES - 1; j++) {
        for (i = 1; i < 252-1; i++) {
            if (!(rgn[j][i] & BADPIXEL)) {
                nPixel = rgn[j][i] >> 7;
                if (nPixel < 32)
                    rgnLow[nPixel]++;
                if (nPixel - 480 >= 0)
                    rgnHigh[nPixel-480]++;
                fError = PixelVariance(rgn, j, i);
                nSigma = int(2.0*fError);
                if (nSigma > 31)
                    nSigma = 31;
                rgnSigmas[nSigma]++;
                fPixel = fabs(fError)/fSigma;
                if (fPixel > 1.0)
                    fPixel = 1.0;
                rgnError[j][i] = int(fPixel * 0xFFF0);
            }
        }
    }
    for (i = 0; i < 32; i++)
        printf("%3d: %7d     %3d: %7d     %0.3f: %7d\n",
                    i, rgnLow[i], i+480, rgnHigh[i], 0.5*double(i), rgnSigmas[i]);
}

void
MaskNoise(unsigned short rgnData[MAXLINES][252], double fSigma)
{
    int i, j, iPorch, nPixel, nMasked;

    printf("      MaskNoise - ");
    //
    //  Some transmissions contain spurious appearences of the digital front porch.
    //
    nMasked = 0;
    for (j = 0; j < MAXLINES; j++) {
        for (i = 2; i < 252-5; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            if (!IsPorch(&rgnData[j][i]))
                continue;
            for (iPorch = 0; iPorch < 5; iPorch++)
                SetBurst(rgnData[j][i+iPorch]);
            nMasked++;
        }
    }
    printf("%dv ", nMasked);
    //
    //  Mask obvious bad pixel values
    //
    nMasked = 0;
    for (j = 0; j < MAXLINES; j++) {
        for (i = 1; i < 252; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            nPixel = rgnData[j][i] >> 7;
            if (nPixel == 0 || nPixel >= 509) {
                SetBad(rgnData[j][i]);
                nMasked++;
            }
        }
    }
    printf("%dz ", nMasked);
    //
    //  Make a high-noise cut-off first, to robustify the statistics of next pass
    //
    nMasked = 0;
    for (j = 1; j < MAXLINES-1; j++) {
        for (i = 0; i < 252; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            nPixel = rgnData[j][i] >> 7;
            if (PixelVariance(rgnData, j, i) > 2.0*fSigma) {
                SetBad(rgnData[j][i]);
                nMasked++;
            }
        }
    }
    printf("%dbb ", nMasked);
    nMasked = 0;
    for (j = 1; j < MAXLINES-1; j++) {
        for (i = 0; i < 252; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            nPixel = rgnData[j][i] >> 7;
            if (PixelVariance(rgnData, j, i) > fSigma) {
                SetBad(rgnData[j][i]);
                nMasked++;
            }
        }
    }
    printf("%db ", nMasked);
    nMasked = 0;
    //
    //  Masking out some noise allows more accurate pixel-deviation recalculation.
    //
    for (j = 1; j < MAXLINES-1; j++) {
        for (i = 0; i < 252; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            nPixel = rgnData[j][i] >> 7;
            if (PixelVariance(rgnData, j, i) > fSigma)
                SetBad(rgnData[j][i]);
        }
    }
    printf("%db2\n", nMasked);
    nMasked = 0;
}
//
//  Assemble a rough master version my judging noisiness of short segments
//
static void
AssembleSegments(unsigned short rgnMaster[MAXLINES][252], int nCopies)
{
    int i, iSeg, j, k, kBest, rgnSrc[4];
    double fChiSqr, fMinChiSqr, fN, fMaxN;
    //
    //  First measure Chi-Square total deviation for short segments of each version
    //  and assemble the best segments.  Don't use anomalous transmissions yet.
    //
    for (k = 0; k < 4; k++)
        rgnSrc[k] = 0;
    for (j = 1; j < MAXLINES - 1; j++) {
        for (iSeg = 0; iSeg < 252; iSeg += NSEG) {
            fMinChiSqr = 1.0e10;
            fMaxN = 0.0;
            kBest = -1;
            for (k = 0; k < nCopies; k++) {
                if (g_rgnFlags[k] & FLAGS_ANOMALOUS)
                    continue;
                fChiSqr = 0.0;
                fN = 0.0;
                for (i = iSeg; i < iSeg+NSEG; i++) {
                    if (i == 0 || i == 251)
                        continue;
                    if (g_rgnRaw[k][j][i] & BADPIXEL)
                        continue;
                    fChiSqr += PixelVariance(g_rgnRaw[k], j, i);
                    fN += 1.0;
                }
                if (fN == 0.0)      // segment all bad?
                    continue;
                fChiSqr /= fN;
                if (fChiSqr/fMinChiSqr < 0.95 || (fChiSqr/fMinChiSqr < 1.05 && fN > fMaxN)) {
                    fMinChiSqr = fChiSqr;
                    fMaxN = fN;
                    kBest = k;      // better data or more data
                }
            }
            if (kBest >= 0) {
                rgnSrc[kBest]++;
                for (i = iSeg; i < iSeg+NSEG; i++) {
                    rgnMaster[j][i] = g_rgnRaw[kBest][j][i];
                    g_rgnSource[j][i] = kBest;
                }
            } else {
                for (i = iSeg; i < iSeg+NSEG; i++) {
                    rgnMaster[j][i] = BADPIXEL;
                    g_rgnSource[j][i] = -1;
                }
            }
        }
    }
    for (k = 0; k < nCopies; k++)
        printf("      %5d segments from %d\n", rgnSrc[k], k);
}

void
MasterVersion(unsigned short rgnMaster[MAXLINES][252], int nCopies, void (*pSpecial)(void),
              int jLimitCleanData, int kPrefered)
{
    double fChiSqr, fMinChiSqr, fChiSqr2, fProbCopy, fProb, fN;
    int i, iCode, iVersion, j, jFirst, jLimit, k, kBest, kCheck, n, nVersions;
    int nIterations, rgnSrc[7], nBurst, nBad, nNotAnomalous, nValue, bAllBurst;
    unsigned short nPixel;
    PixelCount rgcMultiValues[2*NVALUES+2];


    printf("    Master Version\n");
    //
    //  Step 1.  Reject pixels more than 8 sigmas from the mean.
    //
    nNotAnomalous = nCopies;
    for (k = 0; k < nCopies; k++)
        if (g_rgnFlags[k] & FLAGS_ANOMALOUS)
            --nNotAnomalous;
    //for (k = 0; k < nNotAnomalous; k++)
    //    MaskNoise(g_rgnRaw[k], SIGMA_REJECT);
    //
    //  Step 2. Make a rough version of the master image
    //
    AssembleSegments(rgnMaster, nCopies);
    //
    //  Step 3. Recode anomalous transmissions
    //
    for (k = nNotAnomalous; k < nCopies; k++) {
        InitRecoding();
        RecordTransmission(g_rgnRaw[k], rgnMaster, 0, jLimitCleanData);
        BuildRecodingMap(s_rgrmMap[k]);
        Recode(g_rgnRaw[k], rgnMaster, g_rgnBackup[k], s_rgrmMap[k]);
    }
 //ExamineErrors(g_rgnRaw[0], g_rgnRaw[2], 3724, 3727);
    //
    //  Step 4. Fix lines where a zero pixel has been inserted.
    //
    for (k = 0; k < nCopies; k++) {
        if (nCopies == 2)                   // experiments
            kCheck = k ^ 1;
        else if (nCopies == 4)
            kCheck = (k < 2) ? 2 : 0;       // V-14 I, II
        else if (nCopies - nNotAnomalous)
            kCheck = (k == 0) ? 2 : 0;      // V-13 II
        else
            kCheck = (k == 0) ? 1 : 0;      // V-13 I
        //printf("Autoshift comparing %d to %d\n", k, kCheck);
        FixZeroPixels(g_rgnRaw[k], g_rgnRaw[kCheck], k);
    }
    //
    //  Do special fixes and mask out noise
    //
    pSpecial();
    for (k = 0; k < nNotAnomalous; k++)
        MaskNoise(g_rgnRaw[k], SIGMA_REJECT);
    //
    //  Step 5. Redo the rough master after recoding and zero pixel fixes.
    //
    AssembleSegments(rgnMaster, nCopies);
    //
    //  The remaining steps interatively improve the master version
    //
    for (nIterations = 0; nIterations < 4; nIterations++) {
        //
        //  Step 6. Recode anomalous transmissions against latest master
        //
        for (k = nNotAnomalous; k < nCopies; k++) {
            for (j = 0; j < MAXLINES; j++) {
                for (i = 0; i < 252; i++) {
                    g_rgnRaw[k][j][i] = g_rgnBackup[k][j][i];
                }
            }
        }
        for (k = nNotAnomalous; k < nCopies; k++) {
            InitRecoding();
            RecordTransmission(g_rgnRaw[k], rgnMaster, 0, jLimitCleanData);
            BuildRecodingMap(s_rgrmMap[k]);
            Recode(g_rgnRaw[k], rgnMaster, g_rgnBackup[k], s_rgrmMap[k]);
        }
        //
        //  Step 7. Assign probabilities to likely values for each pixel
        //
        for (j = 0; j < MAXLINES; j++) {
            for (i = 0; i < 252; i++) {
                nVersions = 0;
                //
                //  Don't count values masked out as telemetry-burst pixels, unless
                //  all of the version are bursts.
                //
                fN = 0.0;
                bAllBurst = 0;
                for (k = 0; k < nCopies; k++)
                    fN += !IsBurst(g_rgnRaw[k][j][i]);
                if (fN)
                    fProbCopy = 1.0/fN;
                else {
                    fProbCopy = 1.0/double(nCopies);
                    bAllBurst = 1;
                }
                //
                //  All unique possible values and their total probabilities
                //
                for (k = 0; k < nCopies; k++) {
                    if (!bAllBurst && IsBurst(g_rgnRaw[k][j][i]))
                        continue;
                    if (g_rgnFlags[k] & FLAGS_ANOMALOUS) {
                        for (iCode = 0; iCode < NVALUES; iCode++) {
                            nPixel = g_rgnBackup[k][j][i] >> 7;
                            fProb = s_rgrmMap[k][nPixel].rgfProb[iCode];
                            nValue = s_rgrmMap[k][nPixel].rgnValue[iCode];
                            for (iVersion = 0; iVersion < nVersions; iVersion++) {
                                if (nValue == rgcMultiValues[iVersion].nValue) {
                                    rgcMultiValues[iVersion].fScore += float(fProb*fProbCopy);
                                    break;
                                }
                            }
                            if (iVersion == nVersions) {
                                nVersions++;
                                rgcMultiValues[iVersion].fScore = float(fProb*fProbCopy);
                                rgcMultiValues[iVersion].nValue = nValue;
                            }
                        }
                    } else {
                        nValue = g_rgnRaw[k][j][i] & NINEBITS;
                        for (iVersion = 0; iVersion < nVersions; iVersion++) {
                            if (nValue == rgcMultiValues[iVersion].nValue) {
                                rgcMultiValues[iVersion].fScore += float(fProbCopy);
                                break;
                            }
                        }
                        if (iVersion == nVersions) {
                            nVersions++;
                            rgcMultiValues[iVersion].fScore = float(fProbCopy);
                            rgcMultiValues[iVersion].nValue = nValue;
                        }
                    }
                }
                QuickSort(rgcMultiValues, nVersions);
                nPixel = rgcMultiValues[nVersions-1].nValue;
                if (bAllBurst)
                    SetBurst(nPixel);
                g_rgnTest2[j][i] = nPixel;
            }
        }
        //
        //  Fill in bad pixels if any versions of a good pixel
        //
        n = 0;
        for (j = 0; j < MAXLINES; j++) {
            for (i = 0; i < 252; i++) {
                if (rgnMaster[j][i] & BADPIXEL) {
                    for (k = 0; k < nCopies; k++)
                        if (IsGood(g_rgnRaw[k][j][i])) {
                            rgnMaster[j][i] = g_rgnRaw[k][j][i];
                            g_rgnSource[j][i] = k;
                            n++;
                            break;
                        }
                }
            }
        }
        printf("      %d good pixels filled in\n", n);
        //
        //  Sometimes one version is known to be much better
        //
        if (kPrefered >= 0) {
            for (j = 0; j < MAXLINES; j++) {
                for (i = 0; i < 252; i++) {
                    if (!(g_rgnRaw[kPrefered][j][i] & BADPIXEL)) {
                        rgnMaster[j][i] = g_rgnRaw[kPrefered][j][i];
                        g_rgnSource[j][i] = kPrefered;
                    }
                }
            }
            printf("      Moving prefered versions\n");
        }
        //
        //  Now try to replace junky master-version pixels with "better" pixels
        //
        n = 0;
        for (j = 1; j < MAXLINES-1; j++) {
            for (i = 1; i < 252-1; i++) {
                if (!(rgnMaster[j][i] & BADPIXEL)) {
                    fChiSqr = PixelVariance(rgnMaster, j, i);
                    if (fChiSqr > SIGMA_JUNKY) {
                        fMinChiSqr = 100.0;
                        kBest = -1;
                        nPixel = rgnMaster[j][i];
                        for (k = 0; k < nCopies; k++) {
                            if (!(g_rgnRaw[k][j][i] & BADPIXEL)) {
                                rgnMaster[j][i] = g_rgnRaw[k][j][i];
                                fChiSqr2 = PixelVariance(rgnMaster, j, i);
                                if (fChiSqr2 < fMinChiSqr) {
                                    fMinChiSqr = fChiSqr2;
                                    kBest = k;
                                }
                            }
                        }
                        if (kBest >= 0 && fMinChiSqr < fChiSqr) {
                            n++;
                            rgnMaster[j][i] = g_rgnRaw[kBest][j][i];
                            g_rgnSource[j][i] = kBest;
                        } else {
                            rgnMaster[j][i] = nPixel;
                        }
                    }
                }
            }
        }
        printf("      %d pixels improved\n", n);
        //
        //  Mask noise left over in master copy.  Feed mask back to sources.
        //
        MaskNoise(rgnMaster, SIGMA_REJECT);
        for (j = 0; j < MAXLINES; j++)
            for (i = 0; i < 252; i++)
                if (IsBad(rgnMaster[j][i])) {
                    k = g_rgnSource[j][i];
                    if (k >= 0 && k < nCopies)
                        SetBad(g_rgnRaw[k][j][i]);
                }
    }
    for (k = 0; k < 7; k++)
        rgnSrc[k] = 0;
    nBurst = nBad = 0;
    for (jFirst = 2; ; jFirst++)
        if (g_rgnSource[jFirst][10] != -1 || g_rgnSource[jFirst][50] != -1)
            break;
    for (jLimit = MAXLINES-2; ; --jLimit)
        if (g_rgnSource[jLimit-1][10] != -1 || g_rgnSource[jLimit-1][50] != -1)
            break;
    for (j = jFirst; j < jLimit; j++) {
        for (i = 0; i < 252; i++) {
            if (IsBurst(rgnMaster[j][i]))
                nBurst++;
            else if (IsBad(rgnMaster[j][i]))
                nBad++;
            if (g_rgnSource[j][i] == -1)
                rgnSrc[5]++;
            else if (g_rgnSource[j][i] == -4)
                rgnSrc[4]++;
            else
                rgnSrc[g_rgnSource[j][i]]++;
        }
    }
    printf("Master Panorama: 252 x %d\n", jLimit - jFirst);
    for (k = 0; k < 4; k++)
        printf("  %7d pixels from transmission %d\n", rgnSrc[k], k);
    printf("  %7d pixels from Brown University data\n", rgnSrc[SOURCE_BROWN]);
    printf("  %7d pixels are unanimous\n", rgnSrc[SOURCE_UNANIMOUS]);
    printf("  %7d pixels telemetry bursts\n", nBurst);
    printf("  %7d pixels are bad data\n", nBad);
    VisualizeDifferences(nCopies);
}

