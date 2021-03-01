#include "stdafx.h"
#include "ImageFile.h"

void
Sounds()
{
    HANDLE hFile;
    DWORD nBytesWritten;
    short *pn;
    int i;

    hFile = CreateFile("test.wav", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    WriteHeaderWAVE(hFile, 320000, 44100);
    pn = new short[320000];
    for (i = 0; i < 320000; i++)
        pn[i] = int(4000.0*sin(double(i)*0.05));
    WriteFile(hFile, pn, 320000*2, &nBytesWritten, NULL);
    CloseHandle(hFile);
}

