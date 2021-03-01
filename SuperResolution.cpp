#include "stdafx.h"
//
//  Super-resolution calculation
//  D. P. Mitchell  06/18/2003.
//
#include "Venus.h"
#include "ImageFile.h"
#include "ImageProcessing.h"
static void
TestResampling()
{
    static float f, fAve, rgf[209][460], rgf2[209][2*460], rgf3[2*209][2*460], rgfGrad[209];
    static float rgf4[209][230], rgf5[104][230], rgfLap[209][460];
    int i, j;
    unsigned short nU;

    for (j = 0; j < 460; j++) {
        for (i = 0; i < 209; i++) {
            nU = g_rgnMasterV13C1[j+2869-460][i+252-209];
            if (IsBad(nU))
                nU = PixelEstimate(&g_rgnMasterV13C1[j+2869-460][i+252-209], 252, i == 251);
            f = float(UVALUE(nU, 9));
            f = RussianLuminance(f);
            rgf[i][j] = f;
        }
    }
    ContrastExpansion(rgf[0], 209*460);
    for (i = 0; i < 209; i++) {
        rgfGrad[i] = 0.0;
        for (j = 0; j < 460; j++)
            rgfGrad[i] += rgf[i][j];
    }
    for (i = 2; i < 209-2; i++) {
        fAve = (rgfGrad[i-2]+rgfGrad[i-1]+rgfGrad[i]+rgfGrad[i+1]+rgfGrad[i+2])/5.0f;
        f = rgfGrad[i]/fAve;
        for (j = 0; j < 460; j++)
            //rgf[i][j] -= f/460.0; // Doesn't seem to be addative
            rgf[i][j] /= f;
    }
    Laplacian(rgfLap[0], 0.4f, rgf[0], 209, 460);
    for (j = 0; j < 209; j++)
        Magnify(rgf2[j], rgfLap[j], 460, 2);
    for (i = 0; i < 2*460; i++)
        Magnify(&rgf3[0][i], &rgf2[0][i], 209, 2, 2*460);
    WriteFloatImage(rgf3[0], "Resample.bmp", 2*209, 2*460);
    for (j = 0; j < 209; j++)
        Minify(rgf4[j], rgf[j], 460, 2);
    for (i = 0; i < 2*460; i++)
        Minify(&rgf5[0][i], &rgf4[0][i], 209, 2, 460/2);
    WriteFloatImage(rgf5[0], "Resample2.bmp", 209/2, 460/2);
}