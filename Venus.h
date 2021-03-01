#pragma once
//
//  Analysis of Venera 13 and Venera 14 lander images
//  D. P. Mitchell  03/01/2003.
//
//  Formats: Raw 9-bit data is stored in high-order end of unsigned short,
//              with the low-order 7 bits holding flags.
//           Processed pixels are stored in floats with an auxiliary array
//              of chars storing invalid-pixel flag.
//
#define CLEARGLASS  0
#define REDGLASS    1   
#define GREENGLASS  2
#define BLUEGLASS   3

struct Panorama {
    int nGlassFilter;
    int nScanDirection;
    int nWidth;
    float   rgfVGain[1024];
    float   rgfHGain[252];
    float   rgfImage[1024][252];
    char    rgnFlags[1024][252];
};

#define MAXLINES        7500
#define BADPIXEL        1
#define BURSTPIXEL      2
#define RECODEPIXEL     4
#define NINEBITS        0xFF80
#define DATA(N)         ((N) & NINEBITS)
#define ISZERO(N)       (((N) & NINEBITS) == 0 && !IsBurst(N))
#define SIGMA_JUNKY     4.0
#define SIGMA_REJECT    8.0
#define ROT_ROTATE      0
#define ROT_SHIFT       1
#define ROT_HELICAL     2
#define SOURCE_TRANS0    0
#define SOURCE_TRANS1    1
#define SOURCE_TRANS2    2
#define SOURCE_TRANS3    3
#define SOURCE_BROWN     4
#define SOURCE_UNANIMOUS 5
#define SOURCE_UNKNOWN   6
#define FLAGS_NODATA    0
#define FLAGS_TELEMETRY 1
#define FLAGS_ANOMALOUS 2
#define FLAGS_BROWN     4
#define TBC_FULLREPAIR  0
#define TBC_NOSHIFTING  1
#define TBC_ROTATEONLY  2

extern unsigned short   g_rgnRaw[4][MAXLINES][252];         // Up to 4 copies of raw telemetry
extern unsigned short   g_rgnBackup[4][MAXLINES][252];      // Un-recoded anomalous images
extern char             g_rgnSource[MAXLINES][252];         // Origin of master pixels
extern unsigned short   g_rgnMasterV13C1[MAXLINES][252];    // Venera 13, Camera 1
extern unsigned short   g_rgnMasterV13C2[MAXLINES][252];    // Venera 13, Camera 2
extern unsigned short   g_rgnMasterV14C1[MAXLINES][252];    // Venera 14, Camera 1
extern unsigned short   g_rgnMasterV14C2[MAXLINES][252];    // Venera 14, Camera 2
extern unsigned short   g_rgnTest1[MAXLINES][252];          // Visualizeations
extern unsigned short   g_rgnTest2[MAXLINES][252];
extern unsigned short   g_rgnTelemetry[128][4096];          // Telemetry-burst data
extern int              g_jTelemetry;
extern unsigned short   g_rgnBrown[MAXLINES][252];          // 8-bit Brown-University image
extern int              g_rgkBrown[2];                      // Which versions are Brown
extern unsigned short   g_rgnFlags[4];                      // version info
extern int              g_kVerbose;
extern double           g_fShiftLimit;                      // Repair zero shifts without check
extern unsigned short   g_rgnPorch[5];                      // Digital video front porch


extern void SetBad(unsigned short &nPixel);
extern void SetBurst(unsigned short &nPixel);
extern void SetRecode(unsigned short &nPixel);
extern int  IsBad(unsigned short nPixel);
extern int  IsBurst(unsigned short nPixel);
extern int  IsRecode(unsigned short nPixel);
extern int  IsGood(unsigned short nPixel);
extern int  IsPorch(unsigned short rgnData[]);
extern void VisualizeDifferences(int nCopies);
extern void BadPixels(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
                      int iFirst = 0, int iLimit = 252);
extern void RotateColumns(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
                          int nRot, int bShift = ROT_ROTATE, int kBrown = -1);
extern void ShiftColumns(unsigned short rgnData[MAXLINES][252], int jFirst, int nShift,
                         int iFirst = 0);
extern void MaskBurst(int jFirst, int jLimit, int iFirst = 0, int iLimit = 252,
                      int nRotate = 0);
extern void AutoMaskBurst(unsigned short rgnData[MAXLINES][252], int jFirstBurstCenter,
                          int bShort, int jBurstLimit=MAXLINES);
extern int  LoadBrownImage(unsigned short rgnRaw[MAXLINES][252], char *szFile, int jFirst);
extern void FixBrown();
extern unsigned short   BrownPixel(int jRaw, int iRaw, int nMode = 0);
extern int  BrownReplace(int k, int j, int i);
extern int  BrownRepair(int k, int j, int iFirst, int iLimit, int bBadOnly=0);
extern void VisualizeBrown(int jFirst, int jLimit, int nMode = 0);
extern void TimeBaseCorrection(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
                               int jRamp = -1, int nType = TBC_FULLREPAIR,
                               double fCorrLimit = 0.4, int jAuxRamp = -1);
extern void NewTimeBaseCorrection(unsigned short rgnData[MAXLINES][252], int jFirst, int jLimit,
                      int jRamp);
extern void TimeBaseRecode(unsigned short rgnData[MAXLINES][252], int jFirst,
                           int jLimit, int jRamp = -1);
extern void FixZeroPixels(unsigned short rgnFix[MAXLINES][252],
                          unsigned short rgnCheck[MAXLINES][252], int kVersion);
extern void ShiftSubColumn(unsigned short rgnData[MAXLINES][252], int jFirst,
                           int iFirst, int iLimit = 252);
extern void NegShiftSubColumn(unsigned short rgnData[MAXLINES][252], int jFirst,
               int iFirst, int iLimit = 252);
extern void InitRecoding();
extern void RecordTransmission(unsigned short rgnBad[MAXLINES][252],
                               unsigned short rgnGood[MAXLINES][252], int jFirst, int jLimit);
struct PixelCount;
struct RecodingMap;
extern void QuickSort(PixelCount rgtItem[], int nItems);
extern void BuildRecodingMap(RecodingMap rgrm[512]);
extern void Recode(unsigned short rgnData[MAXLINES][252],
                   unsigned short rgnMaster[MAXLINES][252],
                   unsigned short rgnBackup[MAXLINES][252], RecodingMap rgrm[512]);
extern void MaskNoise(unsigned short rgnData[MAXLINES][252], double fSigma);
extern void MasterVersion(unsigned short rgnMaster[MAXLINES][252], int nCopies,
                          void (*pSpecial)(void), int jLimitCleanData, int kPrefered = -1);
extern void RepairHorizon14I();
extern void AlignHorizon14II();
extern void RepairHorizon14II();
extern void BitShift(unsigned short rgnData[], int nShift, int nLength);
extern void BitShiftPixels(unsigned short rgnData[MAXLINES][252], int nShift,
               int jFirst, int jLimit, int iFirst, int iLimit);
extern void BitShiftVenera9Pixels(unsigned short rgnData[512][128], int nShift,
               int jFirst, int jLimit, int iFirst, int iLimit);
extern void SpecialV13C1();
extern void SpecialV13C2();
extern void SpecialV14C1();
extern void SpecialV14C2();
extern void ExamineErrors(unsigned short rgnData[MAXLINES][252], unsigned short rgnCheck[MAXLINES][252],
              int jFirst, int jLimit);
extern void CopyBrown(unsigned short rgnData[MAXLINES][252], int jData, int jBrown,
          int iFirst, int iLimit);


#define Assert(__expr)                                      \
    do {                                                    \
    if (!(__expr)) {                                        \
        printf("Assert failed in %s\n",__FILE__);           \
        printf(" line %d: (%s)\n", __LINE__, #__expr);      \
        __asm { int 3 };                                    \
    }                                                       \
    } while (0)

//
//  Microsoft Programming Conventions:
//
//  g_  global data
//  s_  static data
//  m_  class-member data
//
//  i   index (elements in a line)
//  j   index (lines in a transmission)
//  k   index (transmissions)
//  n   integer
//  f   floating point
//  b   boolean
//  ch  character/byte
//  sz  string (zero terminated)
//  rg  array of
//  p   pointer to
//
//  First, Limit    Limit = Last+1, e.g. for (i = iFirst; i < iLimit; i++)
//  Max, Min
//