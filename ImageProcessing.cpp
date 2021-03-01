#include "stdafx.h"
#pragma intrinsic(sin, cos, fabs, sqrt, pow)
//
//  Basic signal-processing and numerical procedures for images
//  D. P. Mitchell  04/16.2003.
//
#include "Venus.h"
#include "ImageProcessing.h"

//
//  Resample signals with Kaiser-windowed sinc filter
//REVIEW: fix this up to interpolate and decimate 1D, or 2D by integer factors
//

#define D_PI 3.14159265358979323846264

static double
Sinc(double x)
{
    if (x < 0.01 && x > -0.01) {
        x *= D_PI;
        return 1.0 + x*x*(-1.0/6.0 + x*x*(1.0/120.0 + x*x*(-1.0/5040.0 + x*x/362880.0)));
    } else
        return sin(D_PI*x)/(D_PI*x);
}

static double
DSinc(double x)
{
    if (x < 0.001 && x > -0.001) {
        x *= D_PI;
        return D_PI * x*(-2.0/6.0 + x*x*(4.0/120.0 + x*x*(-6.0/5040.0)));
    }
    return (cos(D_PI*x) - Sinc(x))/x;
}

static double
DDSinc(double x)
{
    if (x < 0.001 && x > -0.001) {
        x *= D_PI;
        return D_PI*D_PI * (-1.0/3.0 + x*x*(1.0/10.0 - x*x*(1.0/168.0)));
    }
    return -D_PI*D_PI*Sinc(x) - 2.0*DSinc(x)/x;
}

static float gs_fI0Alpha = 1.0f/11.3019219521363f;
static float gs_fAlphaSquared = 4.0f*4.0f;

static float
Kaiser4Filter(float x, float fLaplacian = 0.0)
{
    float x1, x2, xTerm, fN, fSinc, fKaiser;

    x1 = x * 1.0f/4.0f;  // 1.0/window-half-width
    x2 = x1*x1;
    if (x2 > 1.0)
        return 0.0;
    x2 = 0.25f * gs_fAlphaSquared * (1.0f - x2);    // Alpha set to make sidebands just invisible
    fKaiser = 1.0f + x2;
    fN = 2.0f;
    xTerm = x2;
    while (xTerm > 1.0e-4) {    // Window based on modified Bessel function I0
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
    }
    fSinc = float(Sinc(x));
    //
    //  Don't use this.  Separable laplacian hyper-emphasizes diagonal high frequency.
    //
    if (fLaplacian)
        fSinc -= fLaplacian * float(DDSinc(x));             // add a -omega**2 frequency emphasis
    return fSinc * fKaiser * gs_fI0Alpha;
}
//
//  Interpolate by 2, src data is in upper half of Dst.  This is for Venera 9/10 images.
//
#define PHASE_CORRECTION -0.25f     // 1/2 - 1/2MAGNIFICATION

static float s_rgfSplat[17];

static void
Interpolate1D(float rgfDst[], float rgfSrc[], int n)
{
    int i, k, iSrc, iDst, n2;
    float f;

    for (i = 0; i < n; i++) {
        rgfSrc[i] = rgfDst[i + n];
        rgfDst[i + n] = 0.0;
    }
    n2 = 2*n;
    for (i = -4; i < n+4; i++) {
        iSrc = i;
        if (iSrc < 0)
            iSrc = -iSrc;
        if (iSrc >= n)
            iSrc = n + n - iSrc - 2;    //fixed -1 bug
        f = rgfSrc[iSrc];
        for (k = -8; k <= 8; k++) {
            iDst = 2*i + k;
            if (iDst >= 0 && iDst < n2)
                rgfDst[iDst] += f*s_rgfSplat[k+8];
        }
    }
}

void
Interpolate(float *pfDst, float *pfWorkSpace, int nHigh, int nWide, float fLaplacian)
{
    int j, k;

    for (k = -8; k <= 8; k++)
        s_rgfSplat[k+8] = Kaiser4Filter(float(k)/2.0f + PHASE_CORRECTION, fLaplacian);
    for (j = 0; j < nHigh; j++)
        Interpolate1D(pfDst+2*j*nWide, pfWorkSpace, nWide);
}
//
//  Magnify and Minify are better, more general interpolation routines.  Written later
//  for Venera 13/14 superresolution calculation.
//
#define KHALF   4
#define MAXFAC  32
static int      s_nFactor = 0;
static int      s_nKernel;
static float    s_rgfKernel[2*KHALF*MAXFAC];

static int
InitializeKernel(int nFactor)
{
    float x;
    int k;

    if (nFactor > MAXFAC || nFactor <= 0)
        return 0;
    s_nFactor = nFactor;
    if (nFactor%2 == 1) {
        x = -float(KHALF*nFactor) + 1.0f;   // careful not to introduce a phase shift
        s_nKernel = 2*KHALF*nFactor - 1;
    } else {
        x = -float(KHALF*nFactor) + 0.5f;
        s_nKernel = 2*KHALF*nFactor;
    }
    for (k = 0; x < float(KHALF*nFactor); x += 1.0, k++) {
        s_rgfKernel[k] = Kaiser4Filter(x/float(nFactor), 0.0f);
    }
    s_rgfKernel[k] = 0.0;
    return 1;
}

int
Magnify(float rgfDst[], const float rgfSrc[], int nSrc, int nMagnification, int nStride)
{
    int iSrc, iDst, i, k, nDst;
    float fSrc;

    if (s_nFactor != nMagnification) {
        if (InitializeKernel(nMagnification) == 0)
            return 0;
    }
    nDst = nSrc*nMagnification;
    for (iDst = 0; iDst < nDst; iDst++)
        rgfDst[iDst*nStride] = 0.0f;
    iDst = -KHALF*nMagnification + (nMagnification + 1)/2 - KHALF*nMagnification;
    for (iSrc = -KHALF; iSrc < nSrc+KHALF; iSrc++, iDst += nMagnification) {
        if (iSrc < 0)
            i = -iSrc;
        else if (iSrc > nSrc-1)
            i = 2*nSrc - iSrc - 2;
        else
            i = iSrc;
        fSrc = rgfSrc[i*nStride];
        for (k = 0; k < s_nKernel; k++) {
            i = iDst + k;
            if (i >= 0 && i < nDst)
                rgfDst[i*nStride] += fSrc * s_rgfKernel[k]; // splat
        }
    }
    return 1;
}

int
Minify(float rgfDst[], const float rgfSrc[], int nSrc, int nMinification, int nStride)
{
    int iSrc, iDst, i, k, nDst;
    float fDst;

    if (s_nFactor != nMinification) {
        if (InitializeKernel(nMinification) == 0)
            return 0;
    }
    nDst = nSrc/nMinification;
    for (iDst = 0; iDst < nDst; iDst++)
        rgfDst[iDst*nStride] = 0.0f;
    iSrc = -KHALF*nMinification + (nMinification + 1)/2;
    for (iDst = 0; iDst < nDst; iDst++, iSrc += nMinification) {
        fDst = 0.0;
        for (k = 0; k < s_nKernel; k++) {
            i = iSrc + k;
            if (i < 0)
                i = -i;
            if (i > nSrc-1)
                i = 2*nSrc - i - 2;
            fDst += rgfSrc[i*nStride] * s_rgfKernel[k];     // slurp
        }
        rgfDst[iDst*nStride] += fDst/float(nMinification);
    }
    return 1;
}
//
//  Laplacian sharpening operation (for Venera 9 aperture correction)
//  Not used anymore, laplacian included in windowed filter now.
//
void
Laplacian(float *pfDst, float fA, float *pfSrc, int nHigh, int nWide)
{
    int i, j, ij;
    float fLaplace;

    for (j = 0; j < nHigh; j++) {
        for (i = 0; i < nWide; i++) {
            ij = i + j*nWide;
            if (j == 0 || i == 0 || j == nHigh-1 || i == nWide-1) {
                pfDst[ij] = pfSrc[ij];
                continue;
            }
            fLaplace = 4.0f*pfSrc[ij] - pfSrc[ij+1] - pfSrc[ij-1] - pfSrc[ij+nWide] - pfSrc[ij-nWide];
            pfDst[ij] = pfSrc[ij] + fA*fLaplace;
        }
    }
}

//
//  Repair pixel by isophote interpolation of a pixel from its neighborhood
//
inline float
Lerp(double f1, double f2, double a)
{
    return float(f1 + a*(f2 - f1));
}
//
//  Works with NE, SE, NW, SW neighbors, N S E W can be missing
//
int
IsoPhote(float *pfImage, char *pnMask, int ij, int nWide)
{
    double fDx, fDy, fPx, fPy, fAve, fTan;

    //
    //  Do nothing if not enough good pixels in neighborhood
    //
    if (pnMask[ij-1-nWide] || pnMask[ij+1-nWide] || pnMask[ij-1+nWide] || pnMask[ij+1+nWide])
            return 0;
    //
    //  If working on a column of bad pixels, interpolate across bad pixels, but
    //  don't unset their mask.  They will get replaced by isophote interpolation.
    //
    if (pnMask[ij-1])
        pfImage[ij-1] = Lerp(pfImage[ij-1-nWide], pfImage[ij-1+nWide], 0.5);
    if (pnMask[ij+1])
        pfImage[ij+1] = Lerp(pfImage[ij+1-nWide], pfImage[ij+1+nWide], 0.5);
    if (pnMask[ij-nWide])
        pfImage[ij-nWide] = Lerp(pfImage[ij+1-nWide], pfImage[ij-1-nWide], 0.5);
    if (pnMask[ij+nWide])
        pfImage[ij+nWide] = Lerp(pfImage[ij+1+nWide], pfImage[ij-1+nWide], 0.5);
    //
    //  Isophote direction is just 90 degree rotation of gradient vector
    //
    fDx =  (pfImage[ij+1-nWide]   - pfImage[ij-1-nWide]
      + 2.0*pfImage[ij+1]     - 2.0*pfImage[ij-1]
          + pfImage[ij+1+nWide]   - pfImage[ij-1+nWide])/8.0;
    fDy = (pfImage[ij-1+nWide] + 2.0*pfImage[ij+nWide] + pfImage[ij+1+nWide] 
         - pfImage[ij-1-nWide] - 2.0*pfImage[ij-nWide] - pfImage[ij+1-nWide])/8.0;
    //
    //  If gradient is small, just average the surrounding neighborhood
    //
    if (fabs(fDx) + fabs(fDy) < 1.0e-2) {
        fAve = (pfImage[ij-1-nWide] + 2.0*pfImage[ij-nWide] + pfImage[ij+1-nWide]
          //+ 2.0*pfImage[ij-1]                           + 2.0*pfImage[ij+1]
          +     pfImage[ij-1+nWide] + 2.0*pfImage[ij+nWide] + pfImage[ij+1+nWide])/8.0;
    } else {
        fPx = -fDy;
        fPy = +fDx;
        if (fPx < 0.0) {
            fPx = -fPx;
            fPy = -fPy;
        }
        if (fabs(fPx) > fabs(fPy)) {
            fTan = fPy/fPx;
            if (fTan < 0.0) {
                fAve = Lerp(pfImage[ij+1], pfImage[ij+1-nWide], -fTan)
                     + Lerp(pfImage[ij-1], pfImage[ij-1+nWide], -fTan);
            } else {
                fAve = Lerp(pfImage[ij+1], pfImage[ij+1+nWide], +fTan)
                     + Lerp(pfImage[ij-1], pfImage[ij-1-nWide], +fTan);
            }
        } else {
            fTan = fPx/fPy; // co-tangent
            if (fTan < 0.0) {
                fAve = Lerp(pfImage[ij-nWide], pfImage[ij+1-nWide], -fTan)
                     + Lerp(pfImage[ij+nWide], pfImage[ij-1+nWide], -fTan);
            } else {
                fAve = Lerp(pfImage[ij-nWide], pfImage[ij-1-nWide], +fTan)
                     + Lerp(pfImage[ij+nWide], pfImage[ij+1+nWide], +fTan);
            }
        }
        fAve *= 0.5;
    }
    pfImage[ij] = float(fAve);
    pnMask[ij] = 0;
    return 1;
}
//
//  This version works with N S E W neighbors, and corners can be missing
//
int
IsoPhote2(float *pfImage, char *pnMask, int ij, int nWide)
{
    double fDx, fDy, fPx, fPy, fAve, fTan;

    //
    //  Do nothing if not enough good pixels in neighborhood
    //
    if (pnMask[ij-1] || pnMask[ij+1] || pnMask[ij+nWide] || pnMask[ij-nWide])
            return 0;
    //
    //  If working on a column of bad pixels, interpolate across bad pixels, but
    //  don't unset their mask.  They will get replaced by isophote interpolation.
    //
    if (pnMask[ij-1-nWide])
        pfImage[ij-1-nWide] = Lerp(pfImage[ij-1], pfImage[ij-nWide], 0.5);
    if (pnMask[ij+1-nWide])
        pfImage[ij+1-nWide] = Lerp(pfImage[ij+1], pfImage[ij-nWide], 0.5);
    if (pnMask[ij-1+nWide])
        pfImage[ij-1+nWide] = Lerp(pfImage[ij-1], pfImage[ij+nWide], 0.5);
    if (pnMask[ij+1+nWide])
        pfImage[ij+1+nWide] = Lerp(pfImage[ij+1], pfImage[ij+nWide], 0.5);
    //
    //  Isophote direction is just 90 degree rotation of gradient vector
    //
    fDx = (pfImage[ij+1] - pfImage[ij-1]) * 0.5;
    fDy = (pfImage[ij+nWide] - pfImage[ij-nWide]) * 0.5;
    //
    //  If gradient is small, just average the surrounding neighborhood
    //
    if (fabs(fDx) + fabs(fDy) < 1.0e-2) {
        fAve = (         pfImage[ij-nWide]
          + pfImage[ij-1]                +pfImage[ij+1]
                       + pfImage[ij+nWide])/4.0;
    } else {
        fPx = -fDy;
        fPy = +fDx;
        if (fPx < 0.0) {
            fPx = -fPx;
            fPy = -fPy;
        }
        if (fabs(fPx) > fabs(fPy)) {
            fTan = fPy/fPx;
            if (fTan < 0.0) {
                fAve = Lerp(pfImage[ij+1], pfImage[ij+1-nWide], -fTan)
                     + Lerp(pfImage[ij-1], pfImage[ij-1+nWide], -fTan);
            } else {
                fAve = Lerp(pfImage[ij+1], pfImage[ij+1+nWide], +fTan)
                     + Lerp(pfImage[ij-1], pfImage[ij-1-nWide], +fTan);
            }
        } else {
            fTan = fPx/fPy; // co-tangent
            if (fTan < 0.0) {
                fAve = Lerp(pfImage[ij-nWide], pfImage[ij+1-nWide], -fTan)
                     + Lerp(pfImage[ij+nWide], pfImage[ij-1+nWide], -fTan);
            } else {
                fAve = Lerp(pfImage[ij-nWide], pfImage[ij-1-nWide], +fTan)
                     + Lerp(pfImage[ij+nWide], pfImage[ij+1+nWide], +fTan);
            }
        }
        fAve *= 0.5;
    }
    pfImage[ij] = float(fAve);
    pnMask[ij] = 0;
    return 1;
}
//
//  Try to exprapolate isophotes
//
int
IsoPhote3(float *pfImage, char *pnMask, int ij, int nWide)
{
    double fDx, fDy;

    //
    //  Try to compute the best possible gradient value
    //
    if (pnMask[ij-1] == 0 && pnMask[ij+1] == 0)
        fDx = 0.5*(pfImage[ij+1] - pfImage[ij-1]);
    else if (pnMask[ij] && pnMask[ij-1])
        fDx = pfImage[ij] - pfImage[ij-1];
    else if (pnMask[ij] && pnMask[ij+1])
        fDx = pfImage[ij+1] - pfImage[ij];
    else
        return 0;
    if (pnMask[ij-nWide] == 0 && pnMask[ij+nWide] == 0)
        fDy = 0.5*(pfImage[ij+nWide] - pfImage[ij-nWide]);
    else if (pnMask[ij] && pnMask[ij-nWide])
        fDy = pfImage[ij] - pfImage[ij-nWide];
    else if (pnMask[ij] && pnMask[ij+nWide])
        fDy = pfImage[ij+nWide] - pfImage[ij];
    else
        return 0;

    return 1;
}

int
DiffuseExtrapolate(float *pfImage, char *pnMask, int ij, int nWide)
{
    double fMean, fN;
    int iDelta, jDelta;

    if (pnMask[ij]) {
        fMean = fN = 0.0;
        for (jDelta = -nWide; jDelta <= nWide; jDelta += nWide) {
            for (iDelta = -1; iDelta <= 1; iDelta += 1) {
                if (!pnMask[ij+iDelta+jDelta]) {
                    fMean += pfImage[ij+iDelta+jDelta];
                    fN += 1.0;
                }
            }
        }
        if (fN) {
            fMean /= fN;
            pfImage[ij] = float(fMean);
            pnMask[ij] = 2;
            return 1;
        } else
            return 0;
    }
    return 0;
}

//
//  Normalized cross correlation (ushort version)
//
double
Correlation(unsigned short rgn1[], unsigned short rgn2[], int n)
{
    double fCorr, fN, fMean1, fMean2, fVar1, fVar2, f1, f2, fAccuracy;
    int i;

    fN = fMean1 = fMean2 = fVar1 = fVar2 = fCorr = 0.0;
    for (i = 0; i < n; i++) {
        if ((rgn1[i] & BADPIXEL) || (rgn2[i] & BADPIXEL))
            continue;
        fN += 1.0;
        f1 = UVALUE(rgn1[i], 9);
        f2 = UVALUE(rgn2[i], 9);
        fMean1 += f1;
        fMean2 += f2;
        fVar1 += f1*f1;
        fVar2 += f2*f2;
        fCorr += f1*f2;
    }
    if (fN < 2.0)
        return -1.0;
    fVar1 = (fN*fVar1 - fMean1*fMean1);
    fVar2 = (fN*fVar2 - fMean2*fMean2);
    fCorr = (fN*fCorr - fMean1*fMean2);
    //
    //  Maxing in the instrument error (fAccuracy) prevents dark, low-variance
    //  regions from giving strong bogus correlations.
    //
    fAccuracy = 0.03 - 0.02*fMean2/fN;      // sqrt(Instrument Variance)
    fAccuracy = fAccuracy*fAccuracy*fN*fN;  // "fN*fVar"
    if (fVar1 < fAccuracy)
        fVar1 = fAccuracy;
    if (fVar2 < fAccuracy)
        fVar2 = fAccuracy;
    return fCorr/sqrt(fVar1*fVar2);
}
//
//  Measure the deviation from exact identity of pixels that may be the same
//
double
Similitude(unsigned short rgnData[252], unsigned short rgnRefer[252], int iCount)
{
    int i, n, nEqual;

    nEqual = n = 0;
    for (i = 0; i < iCount; i++) {
        if (IsGood(rgnData[i]) && IsGood(rgnRefer[i])) {
            n++;
            if (rgnData[i] == rgnRefer[i])
                nEqual++;
        }
    }
    if (n)
        return double(nEqual)/double(n) + 1.0e-6;
    else
        return 0.0;
}
//
//  Local pixel deviations from the mean
//
#define HNHOOD  3
#define VNHOOD  1
//
//  Simple estimate of instrument sigma, measured from photoramp in constant-gain
//
static double
Accuracy(double fU)
{
    if (fU < 0.09)
        return 0.009;
    if (fU < 0.2)
        return 0.1*fU;
    if (fU > 0.5)
        return 0.005;
    return 0.02 - 0.015*(fU - 0.2)/(0.5 - 0.2);
}

double
PixelVariance(unsigned short rgnData[MAXLINES][252], int jPixel, int iPixel)
{
    double fN, fU, fMean, fVar;
    int i, j, jLo, jHi, iLo, iHi;

    Assert(iPixel >= 0 && iPixel < 252);
    if (IsBad(rgnData[jPixel][iPixel]))
        return 2.5*SIGMA_REJECT;
    //
    //  The first five pixels of the front porch are fixed digital values.
    //
    if (iPixel < 5) {
        if (rgnData[jPixel][iPixel] == g_rgnPorch[iPixel])
            return 0.0;
        else
            return 2.5*SIGMA_REJECT;
    }
    //
    //  Rows 5 to 40 of the front porch are photometer data from a test pattern,
    //  and the 211 rows below that are image data.  
    //
    if (iPixel >= 5 && iPixel < 41 ) {
        iLo = iHi = iPixel;
        jLo = jPixel - HNHOOD;
        jHi = jPixel + HNHOOD;
    } else {
        iLo = iPixel - VNHOOD;
        iHi = iPixel + VNHOOD;
        if (iHi > 251)
            iHi = 251;
        jLo = jPixel - VNHOOD;
        jHi = jPixel + VNHOOD;
    }
    if (jLo < 0)
        jLo = 0;
    if (jHi >= MAXLINES)
        jHi = MAXLINES - 1;
    //
    //  Return chi = (pixel - mean)/standard_deviation
    //
    fN = fMean = fVar = 0.0;
    for (j = jLo; j <= jHi; j++) {
        if (j == jPixel)
            continue;
        for (i = iLo; i <= iHi; i++) {
            if (IsBad(rgnData[j][i]))
                continue;
            fU = UVALUE(rgnData[j][i], 9);
            fN += 1.0;
            fMean += fU;
            fVar += fU*fU;
        }
    }
    if (fN < 2.0)
        return 0.0;         // not enough statistics to judge a pixel
    fMean /= fN;
    fVar = (fVar - fN*fMean*fMean)/(fN - 1.0);
    if (fVar < 0.0)
        fVar = 0.0;
    fVar = sqrt(fVar);
    if (fVar < 0.02 && fVar < Accuracy(fMean))
        fVar = Accuracy(fMean);
    fU = UVALUE(rgnData[jPixel][iPixel], 9);
    return fabs(fU - fMean)/fVar;
}
//
//  Return a best estimation of a missing pixel value.  bBorder indicates the
//  pixel is on the top or bottom row of the image.
//
unsigned short
PixelEstimate(unsigned short *pn, int nWide, int bBorder)
{
    float fPixel, rgfNHood[3][3];
    char rgnMask[3][3];
    int i, j, jDelta;
    unsigned short nPixel;

    jDelta = bBorder ? 0 : 1;
    for (j = -jDelta; j <= jDelta; j++) {
        for (i = -1; i <= 1; i++) {
            nPixel = pn[i + j*nWide];
            rgnMask[j+1][i+1] = nPixel & BADPIXEL;
            rgfNHood[j+1][i+1] = float(UVALUE(nPixel, 9));
        }
    }
    if (!bBorder && IsoPhote(&rgfNHood[1][1], &rgnMask[1][1], 0, 3))        // isophote interpolation
        fPixel = rgfNHood[1][1];
    else if (!bBorder && IsoPhote2(&rgfNHood[1][1], &rgnMask[1][1], 0, 3))
        fPixel = rgfNHood[1][1];
    else if (rgnMask[1][0] == 0 && rgnMask[1][2] == 0)
        fPixel = Lerp(rgfNHood[1][0], rgfNHood[1][2], 0.5);     // linear interpolation
    else if (rgnMask[1][0] == 0)
        fPixel = rgfNHood[1][0];                                // nearest neighbor
    else if (rgnMask[1][2] == 0)
        fPixel = rgfNHood[1][2];
    else
        return BADPIXEL;
    return NVALUE(fPixel, 9);
}
//
//  Estimate mean and standard deviation of an image region
//
int
GatherStatistics(double &fMean, double &fSigma, double &fRange, unsigned short rgnData[MAXLINES][252],
                 int jFirst, int jLimit, int iFirst, int iLimit)
{
    int i, j;
    double fSum, fVar, fN, fU, fMin, fMax;

    fN = fSum = 0.0;
    fMin = 1.1;
    fMax = -0.1;
    for (j = jFirst; j < jLimit; j++)
        for (i = iFirst; i < iLimit; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            fN += 1.0;
            fU = UVALUE(rgnData[j][i], 9);
            fSum += fU;
            if (fU < fMin) fMin = fU;
            if (fU > fMax) fMax = fU;
        }
    if (fN < 2.0)
        return 0;
    fMean = fSum/fN;
    fRange = fMax - fMin;
    fVar = 0.0;
    for (j = jFirst; j < jLimit; j++)
        for (i = iFirst; i < iLimit; i++) {
            if (rgnData[j][i] & BADPIXEL)
                continue;
            fU = UVALUE(rgnData[j][i], 9);
            fVar += (fU - fMean)*(fU - fMean);
        }
    fSigma = sqrt(fVar/(fN - 1.0));
    return 1;
}
//
//  Convert 9-bit data into linear floating-point format with Russian calibration
//  result.  This is not accurate for darker pixels, and may have the wrong gamma
//  for brighter pixels.
//
float
RussianLuminance(float fU)
{
    float fDensity, fLuminance;

    fDensity = 2.0f*(1.0f - fU);    // Russian calibration (incorrect at low brightness)
    fLuminance = float(pow(10.0, double(-fDensity)));
    return fLuminance;
}

void
ContrastExpansion(float rgf[], int n)
{
    float fMin, fMax;
    int i;

    fMin = fMax = rgf[0];
    for (i = 1; i < n; i++) {
        if (rgf[i] > fMax)
            fMax = rgf[i];
        if (rgf[i] < fMin)
            fMin = rgf[i];
    }
    for (i = 0; i < n; i++)
        rgf[i] = (rgf[i] - fMin)/(fMax - fMin);
}