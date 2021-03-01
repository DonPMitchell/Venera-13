#include "stdafx.h"
#pragma intrinsic(sqrt, fabs)
//
//  InPainting - Bertalmio's algorithm for filling in missing image regions
//  D. P. Mitchell  04/12/2003.
//

//#define NINPAINT    15
//#define NDIFFUSE    2
//#define DELTAT      0.1f
#define NINPAINT    25
#define NDIFFUSE    2
#define INITISO     100
#define ROUNDS      2000
#define DELTAT      0.01f

#define MISSING_REGION  0
#define IMAGE_REGION    1
#define ALL_REGIONS     2
#define ISOTROPIC_DIF   0
#define ANISOTROPIC_DIF 1

static float
Random()
{
    static unsigned nSeed = 1234567;

    nSeed = 1099087573*nSeed + 715136305;
    return float(double(nSeed)/4294967296.0);
}
//
//  Shear the image in a wavey pattern to prevent long straight edges.  In-painting
//  algorithms go crazy on long straight borders, a spurious jet-stream flow forms there.
//
static void
Wavey(float *pfImage, char *pnMask, int nHigh, int nWide, int bInvert, float fPhase = 0.5f, float fAmp=2.0)
{
    int i, j, jDelta;
    float *pf;
    char *pn;

    pf = new float[nHigh];
    pn = new char[nHigh];
    for (i = 0; i < nWide; i++) {
        jDelta = int(fAmp*sin(0.5*double(i) + 2.0f*3.1415926f*fPhase));
        if (bInvert)
            jDelta = -jDelta;
        for (j = 0; j < nHigh; j++) {
            pf[j] = pfImage[i + j*nWide];
            pn[j] = pnMask[i + j*nWide];
        }
        for (j = 0; j < nHigh; j++) {
            pfImage[i + j*nWide] = pf[(j + jDelta + nHigh) % nHigh];
            pnMask[i + j*nWide] = pn[(j + jDelta + nHigh) % nHigh];
        }
    }
    delete [] pf;
    delete [] pn;
}

static void
Laplacian(float *pfLaplacian, float *pfImage, int nHigh, int nWide)
{
    double fDxx, fDyy;
    int i, j, ij;

    for (i = 1; i < nWide-1; i++) {
        for (j = 1; j < nHigh-1; j++) {
            ij = i + j*nWide;
            fDxx = pfImage[ij + 1] - 2.0*pfImage[ij] + pfImage[ij - 1];
            fDyy = pfImage[ij + nWide] - 2.0*pfImage[ij] + pfImage[ij - nWide];
            pfLaplacian[ij] = float(fDxx + fDyy);
        }
    }
}

inline static void
Gradient(double &fDx, double &fDy, float *pfImage, char *pnMask, int ij, int nWide)
{
    fDx = 0.5*(pfImage[ij+1] - pfImage[ij-1]);
    fDy = 0.5*(pfImage[ij+nWide] - pfImage[ij-nWide]);
    return;
    //
    //  Avoid the artificial gradient step at the missing-region boundary
    //
    if (pnMask[ij-1] == pnMask[ij+1])
        fDx = 0.5*(pfImage[ij+1] - pfImage[ij-1]);
    else if (pnMask[ij] == pnMask[ij-1])
        fDx = pfImage[ij] - pfImage[ij-1];
    else
        fDx = pfImage[ij+1] - pfImage[ij];

    if (pnMask[ij-nWide] == pnMask[ij+nWide])
        fDy = 0.5*(pfImage[ij+nWide] - pfImage[ij-nWide]);
    else if (pnMask[ij] == pnMask[ij-nWide])
        fDy = pfImage[ij] - pfImage[ij-nWide];
    else
        fDy = pfImage[ij+nWide] - pfImage[ij];
}
//
//  One step of inpainting, propogate image Laplacian along direction of isophotes
//
static double
InPaint(float *pfDelta, float *pfLaplacian, float *pfImage, char *pnMask, int nHigh, int nWide)
{
    int i, j, ij;
    double fDx, fDy, fNorm, fNx, fNy, fLx, fLy, fDxBack, fDxFore, fDyBack, fDyFore;
    double fMinDxBack, fMaxDxBack, fMinDxFore, fMaxDxFore;
    double fMinDyBack, fMaxDyBack, fMinDyFore, fMaxDyFore;
    double fBeta, fModGradient, fDelta, fTotalDelta;

    for (ij = 0; ij < nWide*nHigh; ij++)
        pfDelta[ij] = 0.0;
    fTotalDelta = 0.0;
    for (i = 1; i < nWide-1; i++) {
        for (j = 1; j < nHigh-1; j++) {
            ij = i + j*nWide;
            if (pnMask[ij] == 0)
                continue;
            //
            //  Compute isophote direction, 90 degrees to image gradient
            //
            Gradient(fDx, fDy, pfImage, pnMask, ij, nWide);
            fNorm = sqrt(fDx*fDx + fDy*fDy + 1.0e-6);
            if (fNorm > 1.0e-4) {
                fNx = fDx/fNorm;
                fNy = fDy/fNorm;
            } else {
                fNx = 0.0;  // no gradient, grab any local pixel?
                fNy = 1.0;
            }
            //
            //  slope-limited gradient
            //
            fDxBack = pfImage[ij] - pfImage[ij - 1];        // forward difference
            fDxFore = pfImage[ij + 1] - pfImage[ij];        // backward difference
            fDyBack = pfImage[ij] - pfImage[ij - nWide];
            fDyFore = pfImage[ij + nWide] - pfImage[ij];
            fMinDxBack = fDxBack < 0.0 ? fDxBack : 0.0;     // min(fDxBack, 0.0)
            fMaxDxBack = fDxBack > 0.0 ? fDxBack : 0.0;     // max(fDxBack, 0.0)
            fMinDxFore = fDxFore < 0.0 ? fDxFore : 0.0;
            fMaxDxFore = fDxFore > 0.0 ? fDxFore : 0.0;
            fMinDyBack = fDyBack < 0.0 ? fDyBack : 0.0;
            fMaxDyBack = fDyBack > 0.0 ? fDyBack : 0.0;
            fMinDyFore = fDyFore < 0.0 ? fDyFore : 0.0;
            fMaxDyFore = fDyFore > 0.0 ? fDyFore : 0.0;
            //
            //  grad(Laplacian) . Isophote_Direction
            //
            Gradient(fLx, fLy, pfLaplacian, pnMask, ij, nWide);
            fBeta = -fNy*fLx + fNx*fLy;
            if (fBeta > 0.0)
                fModGradient = sqrt(fMinDxBack*fMinDxBack + fMaxDxFore*fMaxDxFore +
                                    fMinDyBack*fMinDyBack + fMaxDyFore*fMaxDyFore);
            else
                fModGradient = sqrt(fMaxDxBack*fMaxDxBack + fMinDxFore*fMinDxFore +
                                    fMaxDyBack*fMaxDyBack + fMinDyFore*fMinDyFore);
            fBeta = fBeta * fModGradient;
            fDelta = sqrt(sqrt(fabs(fBeta)));
            if (fBeta < 0.0)
                pfDelta[ij] = -float(fDelta);
            else
                pfDelta[ij] = +float(fDelta);
            fTotalDelta += pfDelta[ij];
        }
    }
    return fTotalDelta;
}
//
//  Anisotropic diffusion along isophote directions.  Bertalmio used a very simple diffusion
//  equation, Laplacian with an inhibition in the gradient direction, no "edge enhancement".
//
#define LAMBDA  0.2
#define THRESH  5.0e-4

static void
Diffuse(float *pfImage, float *pfTemp, char *pnMask, int nHigh, int nWide,
        int nDomain, int nStyle=ANISOTROPIC_DIF, double fAniso=1.00)
{
    int i, j, ij;
    double fDx, fDy, fDxx, fDyy, fDxy, fNorm2, fLaplacian, fAntiGradient;

    for (ij = 0; ij < nWide*nHigh; ij++)
        pfTemp[ij] = pfImage[ij];
    for (i = 1; i < nWide-1; i++) {
        for (j = 1; j < nHigh-1; j++) {
            ij = i + j*nWide;
            if (nDomain == MISSING_REGION && pnMask[ij] == 0)
                continue;
            if (nDomain == IMAGE_REGION && pnMask[ij] != 0)
                continue;
            Gradient(fDx, fDy, pfImage, pnMask, ij, nWide);
            fDxx = pfImage[ij + 1]     - 2.0*pfImage[ij] + pfImage[ij - 1];
            fDyy = pfImage[ij + nWide] - 2.0*pfImage[ij] + pfImage[ij - nWide];
            fDxy = 0.25*( pfImage[ij + 1 + nWide] - pfImage[ij - 1 + nWide]
                        - pfImage[ij + 1 - nWide] + pfImage[ij - i - nWide]);
            fLaplacian = fDxx + fDyy;
            fAntiGradient = fDxx*fDx*fDx + fDyy*fDy*fDy + 2.0*fDxy*fDx*fDy;
            fNorm2 = fDx*fDx + fDy*fDy;
            //
            //  If numerically valid, add gradient-direction-inhibition term
            //
            if (nStyle == ANISOTROPIC_DIF && fNorm2 > THRESH)
                fLaplacian -= fAniso*fAntiGradient/fNorm2;
            pfTemp[ij] = float(LAMBDA*fLaplacian + pfTemp[ij]);
        }
    }
    for (ij = 0; ij < nWide*nHigh; ij++)
        pfImage[ij] = pfTemp[ij];
}
//
//  Bertalmio's isophotic flow in-painting.  These algorithms are touchy and unstable.
//
void
InPainting(float *pfImage, char *pnMask, int nHigh, int nWide, float fPhase)
{
    float *pfLaplacian, *pfDelta;
    double fMaxDelta;
    int ij, nTime, nIteration;

    pfLaplacian = new float[nWide*nHigh];
    pfDelta = new float[nWide*nHigh];
    //
    //  Initialize the missing region with noise, diffuse it around to minimize strong
    //  gradients normal to the region boundary.
    //
    for (ij = 0; ij < nWide*nHigh; ij++) {
        pfLaplacian[ij] = 0.0;
        if (pnMask[ij])
            pfImage[ij] = 0.4f*Random();
    }
    //
    //  Wavey shears the image to remove long straight edges.  Strong spurious flows
    //  can form along the boundary gradient of long straight sections.
    //
    Wavey(pfImage, pnMask, nHigh, nWide, 0, fPhase, 2.0);
    for (nIteration = 0; nIteration < INITISO; nIteration++)
        Diffuse(pfImage, pfDelta, pnMask, nHigh, nWide, MISSING_REGION, ISOTROPIC_DIF);
    for (nTime = 0; ROUNDS; nTime++) {
        for (nIteration = 1; nIteration < NINPAINT; nIteration++) {
            Laplacian(pfLaplacian, pfImage, nHigh, nWide);
            fMaxDelta = InPaint(pfDelta, pfLaplacian, pfImage, pnMask, nHigh, nWide);
            for (ij = 0; ij < nWide*nHigh; ij++)
                pfImage[ij] += DELTAT * pfDelta[ij];
        }
        if (nTime > ROUNDS/NINPAINT)
            break;
        for (nIteration = 0; nIteration < NDIFFUSE; nIteration++)
            Diffuse(pfImage, pfDelta, pnMask, nHigh, nWide, MISSING_REGION, ANISOTROPIC_DIF);
        if (nTime % 10 == 0) {
            printf("%3d delta = %f\n", nTime, fMaxDelta);
        }
    }
    Wavey(pfImage, pnMask, nHigh, nWide, 1, fPhase, 2.0);
    delete [] pfDelta;
    delete [] pfLaplacian;
}
