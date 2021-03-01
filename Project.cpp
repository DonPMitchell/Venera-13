// Project.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Venus.h"

extern void     LoadImages();               // Load and process undecoded 9-bit telemetry
extern void     Radiometric();              // Find radiometric response function and linearize images
extern void     ReadGraphs();               // Read images of graphs and measure function
extern void     ProcessVenera9and10();      // Process images from Venera 9 and Venera 10 panoramas
extern void     ApertureCorrection();       // Analyze and correct camera aperture modulation

int _tmain(int argc, _TCHAR* argv[])
{
    //ReadGraphs();
    //ProcessVenera9and10();
    LoadImages();
    Radiometric();
    //ApertureCorrection();
	return 0;
}

