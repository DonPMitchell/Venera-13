#include "stdafx.h"
//
//  Recode pixel values in anomalous Venera images
//  D. P. Mitchell  05/15/2003.
//
#include "Venus.h"
#include "Recoding.h"
#include "ImageProcessing.h"
//
//  Some values are unchanged, some map to a unique different value, some
//  map probabilistically to just two or three values.
//
RecodingMap         s_rgrmMap[4][512];
static RecodingMap  s_rmZero;
static unsigned     s_rgnHistogram[512][512];   // Histogram[nBad][nGood]

void
InitRecoding()
{
    int i, j;

    for (i = 0; i < 512; i++) {
        for (j = 0; j < 512; j++)
            s_rgnHistogram[j][i] = 0;
    }
}

static void
RecordSample(unsigned short nBad, unsigned short nGood)
{
    if ((nGood & BADPIXEL) || (nBad & BADPIXEL))
        return;
    s_rgnHistogram[nBad >> 7][nGood >> 7]++;
}

void
RecordTransmission(unsigned short rgnBad[MAXLINES][252],
                   unsigned short rgnGood[MAXLINES][252], int jFirst, int jLimit)
{
    int i, j;

    for (j = jFirst; j < jLimit; j++) {
        //
        //  Zero pixels herald shifted scanlines or other trouble, so exclude
        //  those scan lines from the recoding statistics.
        //
        for (i = 1; i < 252; i++)
            if ((rgnBad[j][i] & NINEBITS) == 0 || (rgnGood[j][i] & NINEBITS) == 0)
                goto NextColumn;
        for (i = 0; i < 252; i++)
            RecordSample(rgnBad[j][i], rgnGood[j][i]);
NextColumn:
        ;
    }
}

inline static IsLess(PixelCount &c1, PixelCount &c2)
{
    return c1.fScore < c2.fScore;
}

void
QuickSort(PixelCount rgtItem[], int nItems)
{
    int i, j;
    PixelCount tTemp, tPivot;

    nItems--;   // iLast
    while (nItems >= 30) {
        i = nItems/2;
        //
        //  Sort first, middle and last elements.  This provides median-of-3
        //  partitioning, limits partitioning to only N - 3 remaining items,
        //  and creates sentinals to simplify the inner loop.
        //
        if (IsLess(rgtItem[i], rgtItem[0])) {      // 2.48 compares on average
            tTemp = rgtItem[0];
            rgtItem[0] = rgtItem[i];
            rgtItem[i] = tTemp;
        }
        if (IsLess(rgtItem[nItems], rgtItem[i])) {
            tTemp = rgtItem[nItems];
            rgtItem[nItems] = rgtItem[i];
            if (IsLess(tTemp, rgtItem[0])) {
                rgtItem[i] = rgtItem[0];
                rgtItem[0] = tTemp;
            } else
                rgtItem[i] = tTemp;
        }
        j = nItems - 1;
        tPivot = rgtItem[i];
        rgtItem[i] = rgtItem[j];
        rgtItem[j] = tPivot;
        i = 0;
        //
        //  Partition, using Sedgewick's "j < i" suggestion.  Oddly, it is
        //  faster to loop on i before looping on j (on the Pentium 4).
        //
        for(;;) {
            while(IsLess(rgtItem[++i], tPivot))
                ;
            while(IsLess(tPivot, rgtItem[--j]))
                ;
            if (j < i)
                break;
            tTemp = rgtItem[i];
            rgtItem[i] = rgtItem[j];
            rgtItem[j] = tTemp;
        }
        tTemp = rgtItem[nItems - 1];
        rgtItem[nItems - 1] = rgtItem[i];
        rgtItem[i] = tTemp;
        //
        //  Recursing on smaller partition yields O(log N) stack growth.
        //
        if (j < nItems - i - 1) {
            QuickSort(rgtItem ,j + 1);
            rgtItem += i + 1;
            nItems -= i + 1;
        } else {
            QuickSort(rgtItem + i + 1, nItems - i);
            nItems = j;
        }
    }
    //
    //  Small partitions are insertion sorted.  Distribution is wedge
    //  shaped, with only about 3.8 comparisons done on average, and
    //  benefit gained from structuring the loop for quick out.
    //
    for (i = 1; i <= nItems; i++) {
        j = i;
        tTemp = rgtItem[j];
        if (IsLess(tTemp,rgtItem[j - 1])) {
            do {
                rgtItem[j] = rgtItem[j - 1];
                j = j - 1;
            } while (j > 0 && IsLess(tTemp, rgtItem[j - 1]));
            rgtItem[j] = tTemp;
        }
    }
}

void
BuildRecodingMap(RecodingMap rgrm[512])
{
    unsigned nTotal;
    int jBad, iGood, i;
    PixelCount rgcGoodValues[512];
    //
    //  Sort historgram and calculate bad-pixel recoding probabilites.
    //
    for (jBad = 0; jBad < 512; jBad++) {
        nTotal = 0;
        for (iGood = 0; iGood < 512; iGood++) {
            rgcGoodValues[iGood].fScore = float(s_rgnHistogram[jBad][iGood]);
            rgcGoodValues[iGood].nValue = iGood << 7;
            nTotal += s_rgnHistogram[jBad][iGood];
        }
        QuickSort(rgcGoodValues, 512);
        rgrm[jBad].nCount = nTotal;
        for (i = 0; i < NVALUES; i++) {
            iGood = 511 - i;
            if (nTotal) {
                rgrm[jBad].rgfProb[i] = float(rgcGoodValues[iGood].fScore)/float(nTotal);
                rgrm[jBad].rgnValue[i] = rgcGoodValues[iGood].nValue;
            } else {
                rgrm[jBad].rgfProb[i] = 0.0;
                rgrm[jBad].rgnValue[i] = BADPIXEL;
            }
        }
    }
}

static void
DisplayMap(RecodingMap rgrm[512])
{
    int i, j, nCount;

    for (j = 0; j < 512; j++) {
        printf("%3d->", j);
        for (i = 0; i < NVALUES; i++)
            printf("%3d/%4.2f ", rgrm[j].rgnValue[i]>>7, rgrm[j].rgfProb[i]);
        nCount = rgrm[j].nCount;
        if (nCount > 9999)
            nCount = 9999;
        printf("(%4d) ", nCount);
        if (j % 2 == 1)
            printf("\n");
    }
}

static void
DisplayPixels(int jFirst, int jLimit, int iFirst, int iLimit, int kCopy, int kReference)
{
    int i, j, nOut, nReal, nPixel;

    nOut = 0;
    for (j = jFirst; j < jLimit; j++) {
        for (i = iFirst; i < iLimit; i++) {
            if (IsGood(g_rgnRaw[kReference][j][i])) {
                nReal  = g_rgnRaw[kReference][j][i] >> 7;
                nPixel = g_rgnRaw[kCopy][j][i] >> 7;
                printf("%3d %3d %4d %4d     ", nReal, nPixel, nReal-nPixel, nReal/2 - nPixel);
            } else {
                printf("--- %3d ---- ----     ", nPixel);
            }
            if (nOut % 3 == 2)
                printf("\n");
            nOut++;
        }
    }
    printf("\n");
}

void
Recode(unsigned short rgnData[MAXLINES][252], unsigned short rgnMaster[MAXLINES][252],
       unsigned short rgnBackup[MAXLINES][252], RecodingMap rgrm[512])
{
    int i, j, iValue, nPixel, nNew, nBad, nGood, nRemapped;
    int nSampleSize, bRecoded, nOld, nMaster, nBad2, nRemapped2;
    double fProb, fDeviate, fLimit;
    PixelCount rgcMultiValues[NVALUES];
    //
    //  Recode pixels where there is one confident value
    //
    nBad = nGood = nRemapped = 0;
    for (j = 0; j < MAXLINES; j++) {
        for (i = 0; i < 252; i++) {
            bRecoded = 0;
            rgnBackup[j][i] = rgnData[j][i];    // save so recoding can be repeated
            //
            //  MaskBurst and ShiftSubColumn set BADPIXEL on anomalous images,
            //  but MaskNoise must never be run on them.
            //
            if (IsBad(rgnData[j][i])) {     //REVIEW: all wrong, recode burst pixels and preserve their flags
                SetBurst(rgnData[j][i]);
                continue;
            }
            nPixel = rgnData[j][i] >> 7;
            nOld = nPixel << 7;
            fProb = rgrm[nPixel].rgfProb[0];
            nSampleSize = rgrm[nPixel].nCount;
            //
            //  Does pixel decode into likely value equal to a clean master value?
            //REVIEW: set confidence interval conservatively here
            if (nSampleSize > 25 && PixelVariance(rgnMaster, j, i) < SIGMA_JUNKY) {
                nMaster = rgnMaster[j][i];
                for (iValue = 0; iValue < NVALUES; iValue++) {
                    if (rgrm[nPixel].rgfProb[iValue] >= 0.25 &&
                        rgrm[nPixel].rgnValue[iValue] == nMaster) {
                            bRecoded = 1;
                            rgnData[j][i] = nNew = nMaster;
                            break;
                        }
                }
            }
            //
            //  Is there one highly probably decode value?
            //
            if (!bRecoded) {
                if (nSampleSize < 25)
                    fLimit = 1.0 - 0.05*sqrt(double(nSampleSize))/5.0;
                else
                    fLimit = 0.95;
                if (fProb >= fLimit) {
                    nNew = rgrm[nPixel].rgnValue[0];
                    bRecoded = 1;
                    rgnData[j][i] = nNew;
                }
            }
            if (bRecoded) {
                if (nNew != nOld)
                    nRemapped++;
                nGood++;
            } else {
                nBad++;
                SetRecode(rgnData[j][i]);      // bad pixels are recoding rejections
            }
        }
    }
    //
    //  Initial recoding rejects may be resolvable with local image statistics
    //REVIEW: weight with fDeviate and fProb both
    //REVIEW: confidence limit on rgfProb[NVALUE-1]
    nRemapped2 = nBad2 = 0;
    /*
    for (j = 1; j < MAXLINES-1; j++) {
        for (i = 0; i < 252; i++) {
            nPixel = rgnData[j][i];
            if (!IsRecode(nPixel))
                continue;
            if (rgrm[nPixel>>7].rgfProb[NVALUES-1] > 0.05)  // too many multiple values
                continue;
            nBad2++;
            for (iValue = 0; iValue < NVALUES; iValue++) {
                nNew = rgrm[nPixel>>7].rgnValue[iValue];
                fProb = rgrm[nPixel>>7].rgfProb[iValue];
                rgnData[j][i] = nNew;
                fDeviate = PixelVariance(rgnData, j, i);
                rgnData[j][i] = nPixel;
                rgcMultiValues[iValue].fScore = float(fDeviate);   // weight with fProb somehow
                rgcMultiValues[iValue].nValue = nNew;
            }
            QuickSort(rgcMultiValues, NVALUES);
            if (rgcMultiValues[0].fScore == 0.0)        // not enough local statistics
                continue;
            if (rgcMultiValues[0].fScore < 1.5 && rgcMultiValues[1].fScore > 1.5) {
                //
                //  Only one recode value is plausible in image context
                //
                nRemapped2++;
                nNew = rgcMultiValues[0].nValue;
                rgnData[j][i] = nNew;
            }
        }
    }
    */
    printf("      Recode: %d bad, %d good, %d+%d remapped  %5.2f percent salvaged\n",
        nBad, nGood, nRemapped, nRemapped2, 100*double(nGood+nRemapped2)/double(nGood+nBad));
}
