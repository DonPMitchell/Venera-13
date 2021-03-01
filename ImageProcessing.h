#pragma once
//
//  Basic Image Processing routines
//  D. P. Mitchell  04/16/2003
//
#define UVALUE(N, B) (double((N) >> (16-B))/double((1 << B) - 1))
#define NVALUE(U, B) (unsigned short(U * double((1 << B) - 1) + 0.5) << (16-B))
#define NHOOD8  1
#define NHOOD6  2
#define NHOOD2  3

struct ConfidenceInterval {
    float   fLower;
    float   fUpper;
};

extern void     Interpolate(float *pfDst, float *pfWorkSpace, int nHigh, int nWide, float fLaplacian = 0.0f);
extern void     Laplacian(float *pfDst, float fSharpness, float *pfSrc, int nHigh, int nWide);
extern float    Lerp(double f1, double f2, double a);
extern int      IsoPhote(float *pfImage, char *pnMask, int ij, int nWide);
extern int      IsoPhote2(float *pfImage, char *pnMask, int ij, int nWide);
extern int      DiffuseExtrapolate(float *pfImage, char *pnMask, int ij, int nWide);
extern double   Correlation(unsigned short rgn1[], unsigned short rgn2[], int n);
extern double   Similitude(unsigned short rgnData[252], unsigned short rgnRefer[252], int iCount);
extern double   PixelVariance(unsigned short rgnData[MAXLINES][252], int jPixel, int iPixel);
extern unsigned short PixelEstimate(unsigned short *pn, int nWide, int bBorder = 0);
extern int      GatherStatistics(double &fMean, double &fSigma, double &fRange, unsigned short rgn[MAXLINES][252],
                                 int jFirst, int jLimit, int iFirst = 0, int iLimit = 252);
extern int      Magnify(float rgfDst[], const float rgfSrc[], int nSrc, int nMagnification, int nStride=1);
extern int      Minify(float rgfDst[], const float rgfSrc[], int nSrc, int nMinification, int nStride=1);
extern void     ImageFFT(float *pfPower, float *pfPhase, const float *pfImage, int nHigh, int nWide);
extern void     ImageInverseFFT(float *pfImage, float *pfPower, float *pfPhase, int nHigh, int nWide);
extern float    RussianLuminance(float fU);
extern void     ContrastExpansion(float rgf[], int n);

