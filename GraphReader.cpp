#include "stdafx.h"
//
//  Simple routine to read data off images of a graph.  Scale it to exactly 100 times
//  the desired number of data values, keep image wide but flat so slopes are small.
//  Write out in raw mode, 8-bit, greyscale.  Image must be inverted, the graph being
//  white on black.
//
#include "ImageFile.h"

static int
GraphReader(float rgyData[], char *szFile, int nWidth, double yLow, double yHigh, int nRate = 100)
{
    HANDLE hFile;
    DWORD nBytesRead;
    int i, j;
    unsigned char *pch;
    double yMin, yMax, *pfWeight;

    hFile = CreateFile(szFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                       FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("%s not found by GraphReader\n", szFile);
        return 0;
    }
    pch = new unsigned char[nWidth];
    pfWeight = new double[nWidth/nRate];
    for (i = 0; i < nWidth/nRate; i++) {
        rgyData[i] = 0.0;
        pfWeight[i] = 0.0;
    }
    for (j = 0; ; --j) {
        ReadFile(hFile, pch, nWidth, &nBytesRead, NULL);
        if (nBytesRead != nWidth)
            break;
        for (i = nRate/2; i < nWidth; i += nRate) {
            rgyData[i/nRate] += float(double(j)*double(pch[i-1] + pch[i] + pch[i+1]));
            pfWeight[i/nRate] += double(pch[i-1] + pch[i] + pch[i+1]);
        }
    }
    CloseHandle(hFile);
    delete [] pch;
    yMin = yMax = rgyData[0]/pfWeight[0];
    for (i = 0; i < nWidth/nRate; i++) {
        rgyData[i] = float(rgyData[i]/pfWeight[i]);
        if (rgyData[i] < yMin)
            yMin = rgyData[i];
        if (rgyData[i] > yMax)
            yMax = rgyData[i];
    }
    delete [] pfWeight;
    for (i = 0; i < nWidth/nRate; i++)
        rgyData[i] = float(yLow + (yHigh - yLow)*(rgyData[i] - yMin)/(yMax - yMin));
    return 1;
}

static void
PrintData(float rgf[], int n, char *szName)
{
    int i;

    printf("double %s[%d] = {\n", szName, n);
    for (i = 0; i < n; i++) {
        printf("%10f,", rgf[i]);
        if (i%6 == 5)
            printf("\n");
    }
    printf("\n};\n");
}

#define SAMPLERATE 44100
#define PERSAMPLE ((SAMPLERATE/3)/2)
#define FILTER 10

static float
Random()
{
    static unsigned nSeed = 1234567;

    nSeed = 1099087573*nSeed + 715136305;
    return float(double(nSeed)/4294967296.0);
}

static short rgnSamples[PERSAMPLE];
static short rgnFiltered[PERSAMPLE/FILTER];

static void
AudioGraphs()
{
    float rgfLevel[600];
    HANDLE hFile;
    DWORD nBytesWritten;
    int i, j, k, nSample;

    hFile = CreateFile("..//ProjectResults//Venera14.wav", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    WriteHeaderWAVE(hFile, PERSAMPLE*600/FILTER, SAMPLERATE);
    GraphReader(rgfLevel, "VeneraMicrophone.raw", 30000, 1.5-1.0, 6.5-1.0, 50);
    for (i = 0; i < 600; i++) {
        for (j = 0; j < PERSAMPLE; j++)
            rgnSamples[j] = int(1500.0 * rgfLevel[i] * (Random() - 0.5));
        for (j = 0; j < PERSAMPLE; j += FILTER) {
            nSample = 0;
            for (k = 0; k < FILTER; k++)    // filter to about 2200 Hz, the microphone limit
                nSample += rgnSamples[j + k];
            rgnFiltered[j/FILTER] = nSample;
        }
        WriteFile(hFile, rgnSamples, sizeof(short)*PERSAMPLE/FILTER, &nBytesWritten, NULL);
    }
    CloseHandle(hFile);
}

void
ReadGraphs()
{
    static float rgf[512];

    printf("0. Read Graphical Data\n");
    printf("  A. Venera 9 and 10 radiometric response.\n");
    GraphReader(rgf, "Venera10_UtoD.raw", 6400, 0.0, 1.146);
    PrintData(rgf, 64, "V10Response");
    GraphReader(rgf, "Venera9_UtoD.raw", 6400, 0.0, 1.440);
    PrintData(rgf, 64, "V9Response");
    printf("  B. Venera 14 audio.\n");
    AudioGraphs();
}