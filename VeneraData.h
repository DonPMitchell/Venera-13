#pragma once
//
//  Global data tables for Venera project
//  D. P. Mitchell  05/30/2003.
//

//
//  Spectral response of cameras with colored glass filters
//
#define NGLASS 37

extern float    g_rgfGlassWavelengths[NGLASS];
extern float    g_rgfBlue_13I[NGLASS];
extern float    g_rgfGreen_13I[NGLASS];
extern float    g_rgfRed_13I[NGLASS];
extern float    g_rgfBlue_13II[NGLASS];
extern float    g_rgfGreen_13II[NGLASS];
extern float    g_rgfRed_13II[NGLASS];
extern float    g_rgfBlue_14I[NGLASS];
extern float    g_rgfGreen_14I[NGLASS];
extern float    g_rgfRed_14I[NGLASS];
extern float    g_rgfBlue_14II[NGLASS];
extern float    g_rgfGreen_14II[NGLASS];
extern float    g_rgfRed_14II[NGLASS];

#define NPAINT 33

extern float    g_rgfPaintWavelengths[NPAINT];
extern float    g_rgfBluePaint[NPAINT];
extern float    g_rgfGreenPaint[NPAINT];
extern float    g_rgfRedPaint[NPAINT];
extern float    g_rgfWhitePaint[NPAINT];
extern float    g_rgfGrayPaint[NPAINT];
extern float    g_rgfBlackPaint[NPAINT];
extern float    g_rgfTeeth[NPAINT];         // titanium teeth on landing ring

#define NCAMERA 21

extern float    g_rgfOpticalDensity[NCAMERA];
extern float    g_rgfSignal_13I[NCAMERA];
extern float    g_rgfSignal_13II[NCAMERA];
extern float    g_rgfSignal_14I[NCAMERA];
extern float    g_rgfSignal_14II[NCAMERA];
