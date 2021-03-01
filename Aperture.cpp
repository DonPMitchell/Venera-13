#include "stdafx.h"
//
//  Calculation of camera aperture modulation
//  D. P. Mitchell  06/11/2003
//
#include "Venus.h"
#include "ImageProcessing.h"
#include "ImageFile.h"

#define TPI 6.283185307179586476925286766559

static void
HammingWindow(float *pf, int nHigh, int nWide)
{
    int i, j;
    double x, y;
    
    for (j = 0; j < nHigh; j++) {
        y = 1.0 - cos(TPI*double(j)/double(nHigh-1));
        for (i = 0; i < nWide; i++) {
            x = 1.0 - cos(TPI*double(i)/double(nWide-1));
            pf[i + j*nWide] *= 0.25f*float(x*y);
        }
    }
}

static void
RadialSpectrum(float rgfSpec[], float *pf, int nWide, float fScale)
{
    int i, j, ir, nHalf, rgn[1000];
    double r, xDelta, yDelta;
    
    for (i = 0; i < nWide; i++) {
        rgfSpec[i] = 0.0;
        rgn[i] = 0;
    }
    nHalf = nWide/2;
    for (j = 0; j < nWide; j++) {
        for (i = 0; i < nWide; i++) {
            xDelta = (i - nHalf)/fScale;
            yDelta = j - nHalf;
            r = sqrt(xDelta*xDelta + yDelta*yDelta);
            ir = int(r + 0.5);
            if (ir < nWide) {
                rgn[ir]++;
                rgfSpec[ir] += pf[i + j*nWide];
            }
        }
    }
    for (i = 0; i < nWide; i++)
        rgfSpec[i] /= float(rgn[i]);
}

#define N 128

static void
TestFFT(int jFirst, int iFirst, float rgfImage[N][N], float rgfPower[N][N])
{
    float rgfPhase[N][N];
    float rgfCombo[N][3*N];
    float rgfSpec[N];
    int i, j;
    unsigned short nU;
    float f;

    for (j = 0; j < N; j++) {
        for (i = 0; i < N; i++) {
            nU = g_rgnMasterV13C1[j+jFirst][i+iFirst];
            if (IsBad(nU))
                nU = 0;
            f = float(UVALUE(nU, 9));
            f = RussianLuminance(f);
            rgfImage[i][j] = f;
        }
    }
    //ContrastExpansion(rgfImage[0], N*N);
    HammingWindow(rgfImage[0], N, N);
    ImageFFT(rgfPower[0], rgfPhase[0], rgfImage[0], N, N);
    RadialSpectrum(rgfSpec, rgfPower[0], N, 0.6667f);
    for (j = 0; j < N; j++)  {
        for (i = 0; i < N; i++) {
            rgfCombo[i][j] = rgfImage[i][j];
            rgfCombo[i][j+N] = rgfPower[i][j];
            rgfCombo[i][j+N+N] = 0.0f;
        }
    }
    for (j = 0; j < N; j++) {
        f = rgfPower[j][N/2];
        if (f > 1.0)
            f = 1.0;
        for (i = int(66.0*f); i >= 0; --i)
            rgfCombo[N-i-1][j+N+N] = 1.0f;
        f = rgfPower[N/2][j];
        if (f > 1.0)
            f = 1.0;
        for (i = int(f*N/3); i >= 0; --i)
            rgfCombo[N-i-N/3][j+N+N] = 0.8f;
        f = rgfSpec[j/2];
        if (f > 1.0)
            f = 1.0;
        for (i = int(f*N/3); i >= 0; --i) {
            rgfCombo[N-i-2*N/3][j/2+N+N+N/2] = 0.6f;
            rgfCombo[N-i-2*N/3][N+N+N/2-j/2] = 0.6f;
        }
    }
    WriteFloatImage(rgfCombo[0], "FFT.bmp", N, 3*N);
}
//
//  Generate images of sections and their power spectrum, average many spectra together
//
static void
Periodogram()
{
    float rgfImage[N][N], rgfPower[N][N], rgfTotal[N][N];
    static float rgfSpectrum[2*N][2*N], rgf[2*N][7*N];
    int i, j, k, nSpectra;

    for (i = 0; i < 7*N; i++)
        for (j = 0; j < 2*N; j++)
            rgf[j][i] = 0.0f;
    //
    //  Image spectra of 6 typical areas
    //
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            rgfTotal[j][i] = 0.0f;
    for (k = -N; k <= N; k += N) {
        TestFFT(3398-N/2+k, 252-N, rgfImage, rgfPower);
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i] / 6.0f;
                rgf[j][i+N+k] = rgfImage[j][i];
                rgf[j+N][i+N+k] = rgfPower[j][i];
            }
        }
    }
    for (k = -N; k <= N; k += N) {
        TestFFT(3398-N/2+k, 252-200, rgfImage, rgfPower);
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i] / 6.0f;
                rgf[j][i+4*N+k] = rgfImage[j][i];
                rgf[j+N][i+4*N+k] = rgfPower[j][i];
            }
        }
    }
    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++)
            rgf[j][i+6*N] = rgfTotal[j][i];
    //
    //  Average a periodogram over a large number of samples
    //
    nSpectra = 0;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            rgfTotal[j][i] = 0.0f;
    for (k = 490; k <= 706-N; k += N/4) {
        TestFFT(k, 252-N, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 718; k <= 950-N; k += N/4) {
        TestFFT(k, 252-N, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 3130; k <= 3564-N; k += N/4) {
        TestFFT(k, 252-N, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 5800; k <= 6004-N; k += N/4) {
        TestFFT(k, 252-N, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 6030; k <= 6330-N; k += N/4) {
        TestFFT(k, 252-N, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 252-N-N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 2220; k <= 2394-N; k += N/4) {
        TestFFT(k, 45, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 45+N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 45+N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 2408; k <= 2606-N; k += N/4) {
        TestFFT(k, 45, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 45+N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 45+N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (k = 4900; k <= 5038-N; k += N/4) {
        TestFFT(k, 45, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 45+N/4, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
        TestFFT(k, 45+N/2, rgfImage, rgfPower);
        nSpectra++;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                rgfTotal[j][i] += rgfPower[j][i];
            }
        }
    }
    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++)
            rgf[j+N][i+6*N] = rgfTotal[j][i]/float(nSpectra);
    WriteFloatImage(rgf[0], "Aperture.bmp", 2*N, 7*N);
}

void
ApertureCorrection()
{
    printf("Aperture Correction\n");
    Periodogram();
}
