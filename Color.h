#pragma once
//
//  Color - Managing spectral and trichromatic color data
//  D. P. Mitchell  03/30/2003.
//
#define SPECTRUM_START  360
#define SPECTRUM_END    830
#define SPECTRUM_STEP     5
#define SPECTRUM_NDATA  ((SPECTRUM_END - SPECTRUM_START)/SPECTRUM_STEP + 1)
#define SPECTRUM_CIE    (SPECTRUM_END - SPECTRUM_START + 1)

extern float sdVenusSky[SPECTRUM_NDATA];

