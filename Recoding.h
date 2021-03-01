#pragma once
//
//  Recode anomalous pixel values
//
#define NVALUES 4

struct RecodingMap {
    unsigned        nCount;
    float           rgfProb[NVALUES];
    unsigned short  rgnValue[NVALUES];
};

struct PixelCount {
    float           fScore;
    unsigned short  nValue;
};

extern RecodingMap  s_rgrmMap[4][512];