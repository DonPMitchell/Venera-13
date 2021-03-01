#include "stdafx.h"
//
//  Convert unencoded 9-bit data into linear radiance values
//  D. P. Mitchell  03/03/2003.
//
#include "Venus.h"
#include "ImageProcessing.h"
#include "ImageFile.h"
//
//  Chartless radiometric camera calibration [
//
#define NWIDE 840
#define NHIGH 256
static unsigned s_rgnPlot[NHIGH][NWIDE];
static unsigned s_rgnT[512][512];

static unsigned short s_rgnRed[42];         // calibration from constant-gain mode
static unsigned short s_rgnGreen[42];
static unsigned short s_rgnBlue[42];
static unsigned short s_rgnClear[42];
static unsigned short s_rgnAutoGain[42];    // calibration in automatic-gain mode

static void
PlotWedge(unsigned short rgnWedge[42], unsigned short rgnBase[42], unsigned nColor)
{
    int i, iDelta, nPixel, nPixel0, nPixel1, nBase0, nBase1, nBase, nDiff;

    for (i = 5; i < 42 - 1; i++) {
        nPixel0 = rgnWedge[i];
        nPixel1 = rgnWedge[i+1];
        if (IsBad(nPixel0) || IsBad(nPixel1))
            continue;
        nPixel0 >>= 7;
        nPixel1 >>= 7;
        nBase0 = rgnBase[i] >> 7;
        nBase1 = rgnBase[i+1] >> 7;
        for (iDelta = 0; iDelta < 10; iDelta++) {
            nPixel = (nPixel0*(10 - iDelta) + nPixel1*iDelta) / 10;
            nBase =  ( nBase0*(10 - iDelta) +  nBase1*iDelta) / 10;
            nDiff = (nPixel - nBase) + 256;
            if (nDiff < 0) nDiff = 0;
            if (nDiff > 511) nDiff = 511;
            s_rgnPlot[255 - nPixel/2][10*i + iDelta] = nColor;
            s_rgnPlot[255 - nDiff/2][10*i + iDelta + 420] = nColor;
        }
    }
}

static void
PlotWedgeInteval(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit, unsigned nColor)
{
    int j;

    for (j = jFirst; j < jLimit; j++)
        PlotWedge(rgnData[j], s_rgnAutoGain, nColor);
}

static void
AverageWedge(unsigned short rgnAve[42], unsigned short rgnData[MAXLINES][252],
             int jFirst, int jLimit)
{
    int i, j, n, nSum;

    for (i = 0; i < 42; i++)
        rgnAve[i] = 0;
    for (i = 5; i < 42; i++) {
        n = nSum = 0;
        for (j = jFirst; j < jLimit; j++) {
            if (IsBad(rgnData[j][i]))
                continue;
            n++;
            nSum += rgnData[j][i] >> 7;
        }
        rgnAve[i] = (nSum/n) << 7;
    }
}
//
//  Plot same input, same gain, different exposures (glass filter densities)
//
static void
PlotVaryingLight()
{
    int i, j;
    unsigned nClear, nR, nG, nB;

    for (j = 0; j < NHIGH; j++)
        for (i = 0; i < NWIDE; i++)
            s_rgnPlot[j][i] = 0;
    nClear = ColorPixel(255, 255, 255);
    nR = ColorPixel(255, 100, 100);
    nG = ColorPixel(100, 255, 100);
    nB = ColorPixel(100, 100, 255);
    PlotWedge(s_rgnClear, s_rgnRed, nClear);      // clear glass filter
    PlotWedge(s_rgnRed, s_rgnRed, nR);            // red glass filter
    PlotWedge(s_rgnGreen, s_rgnRed, nG);          // green glass filter
    PlotWedge(s_rgnBlue, s_rgnRed, nB);           // blue glass filter
    WriteLongImage(s_rgnPlot[0], "VaryExposure.bmp", NHIGH, NWIDE);
}
//
//  Plot same input, same exposure, difference phototube gains
//
static void
PlotVaryingGain(unsigned short rgnData[MAXLINES][252], int jFirstBurn, int jLimitBurn,
                int jFirstDec, int jLimitDec)
{
    int i, j;
    unsigned nClear;
    unsigned short rgn[42];

    for (j = 0; j < NHIGH; j++)
        for (i = 0; i < NWIDE; i++)
            s_rgnPlot[j][i] = 0;
    nClear = ColorPixel(255, 255, 100);
    PlotWedgeInteval(rgnData, jFirstBurn, jLimitBurn, nClear);        // clear glass filter, automatic gain
    nClear = ColorPixel(255, 150, 100);
    for (j = jFirstDec; j < jLimitDec; j += 20) {
        AverageWedge(rgn, rgnData, j - 1, j + 2);
        PlotWedge(rgn, s_rgnAutoGain, nClear);
    }
    WriteLongImage(s_rgnPlot[0], "VaryGain.bmp", NHIGH, NWIDE);
}

static void
PlotTransfer(unsigned short rgnWedgeX[42], unsigned short rgnWedgeY[42], unsigned nColor)
{
    int i, iDelta, nX, nX0, nX1, nY, nY0, nY1;

    for (i = 9; i < 33 - 1; i++) {
        nX0 = rgnWedgeX[i] >> 7;
        nX1 = rgnWedgeX[i+1] >> 7;
        nY0 = rgnWedgeY[i] >> 7;
        nY1 = rgnWedgeY[i+1] >> 7;

        for (iDelta = 0; iDelta < 10; iDelta++) {
            nX = (nX0*(10 - iDelta) + nX1*iDelta) / 10;
            nY =  ( nY0*(10 - iDelta) +  nY1*iDelta) / 10;
            s_rgnT[511-nY][nX] = nColor;
        }
    }
}

static void
BrightnessTransfer()
{
    int nX, nY;

    for (nX = 0; nX < 512; nX++)
        for (nY = 0; nY < 512; nY++)
            s_rgnT[nX][nY] = 0;
    PlotTransfer(s_rgnGreen, s_rgnBlue, ColorPixel(100, 255, 255));
    PlotTransfer(s_rgnGreen, s_rgnRed, ColorPixel(255, 255, 100));
    PlotTransfer(s_rgnGreen, s_rgnClear, ColorPixel(100, 255, 100));
    WriteLongImage(s_rgnT[0], "BrightTransfer.bmp", 512, 512);
}

static void
AnalyzeV13C1()
{
    printf("  A. Analyze Venera 13, Camera I\n");
    //
    //  Gather averaged data on calibration signal in reverse-direction scanning.
    //
    AverageWedge(s_rgnClear, g_rgnMasterV13C1, 1198, 1203);     // clear glass filter
    AverageWedge(s_rgnRed, g_rgnMasterV13C1, 1205, 1520);       // red glass filter
    AverageWedge(s_rgnGreen, g_rgnMasterV13C1, 1859, 1864);     // green glass filter
    AverageWedge(s_rgnBlue, g_rgnMasterV13C1, 1865, 2185);      // blue glass filter
    AverageWedge(s_rgnAutoGain, g_rgnMasterV13C1, 2566, 2666);  // clear, auto-gain
    PlotVaryingGain(g_rgnMasterV13C1, 2217, 2227, 2768, 2870);
    PlotVaryingLight();
    BrightnessTransfer();
}

extern void TestLSE();
extern void TestQRCP();
extern void TestFit();

void
Radiometric()
{
    printf("II. Linearize images\n");
    TestFit();
    AnalyzeV13C1();
}