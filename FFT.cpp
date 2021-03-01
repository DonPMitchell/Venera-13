#include "stdafx.h"
#pragma intrinsic(log, exp, sqrt, atan2, cos, sin)
//
//  Fourier transform of a 2D array (image)
//  D. P. Mitchell  06/09/2003.
//
#define TPI 6.283185307179586476925286766559

//
//  General-Radix FFT
//
typedef struct Complex {
	double  r, i;
} Complex;

static Complex zero = {0.0, 0.0};
static AB;
static rgnPrime[] = {
   4,   2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,
  43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103,
 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,   1
};
static Complex z[10000];

static void
SlowTransform(Complex x0[], Complex x1[], int A, int B, Complex root)
{
	int aB, alphaB;
	int B2, B3, B4;
	Complex w0, w1, t, t1, t2, t3, t4, t5, m0, m1, m2, m3, m4, m5;

	switch (A) {

	case 2:
		x1[0].r = x0[0].r + x0[B].r;
		x1[0].i = x0[0].i + x0[B].i;
		x1[B].r = x0[0].r - x0[B].r;
		x1[B].i = x0[0].i - x0[B].i;
		break;
	case 3:
		B2 = B + B;
		t1.r = x0[B].r + x0[B2].r;
		t1.i = x0[B].i + x0[B2].i;
		m0.r = x0[0].r + t1.r;
		m0.i = x0[0].i + t1.i;
		m1.r = (root.r - 1.0)*t1.r + m0.r;
		m1.i = (root.r - 1.0)*t1.i + m0.i;
		m2.r = root.i * (x0[B2].i - x0[B].i);
		m2.i = root.i * (x0[B].r - x0[B2].r);
		x1[0] = m0;
		x1[B].r = m1.r + m2.r;
		x1[B].i = m1.i + m2.i;
		x1[B2].r = m1.r - m2.r;
		x1[B2].i = m1.i - m2.i;
		break;
	case 4:
		B2 = B + B;
		B3 = B2 + B;
		t1.r = x0[0].r + x0[B2].r;
		t1.i = x0[0].i + x0[B2].i;
		t2.r = x0[B].r + x0[B3].r;
		t2.i = x0[B].i + x0[B3].i;
		x1[0].r = t1.r + t2.r;
		x1[0].i = t1.i + t2.i;
		x1[B2].r = t1.r - t2.r;
		x1[B2].i = t1.i - t2.i;
		m2.r = x0[0].r - x0[B2].r;
		m2.i = x0[0].i - x0[B2].i;
		m3.r = x0[B3].i - x0[B].i;
		m3.i = x0[B].r - x0[B3].r;
		x1[B].r = m2.r + m3.r;
		x1[B].i = m2.i + m3.i;
		x1[B3].r = m2.r - m3.r;
		x1[B3].i = m2.i - m3.i;
		break;
	case 5:
		B2 = B + B;
		B3 = B2 + B;
		B4 = B3 + B;
		w0.r = root.r * root.r - root.i * root.i;
		w0.i = root.r * root.i;
		w0.i += w0.i;
		t1.r = x0[B].r + x0[B4].r;
		t1.i = x0[B].i + x0[B4].i;
		t2.r = x0[B2].r + x0[B3].r;
		t2.i = x0[B2].i + x0[B3].i;
		t3.r = x0[B].r - x0[B4].r;
		t3.i = x0[B].i - x0[B4].i;
		t4.r = x0[B3].r - x0[B2].r;
		t4.i = x0[B3].i - x0[B2].i;
		t5.r = t1.r + t2.r;
		t5.i = t1.i + t2.i;
		m0.r = x0[0].r + t5.r;
		m0.i = x0[0].i + t5.i;
		m1.i = m1.r = 0.5 * (root.r + w0.r) - 1.0;
		m1.r *= t5.r;
		m1.i *= t5.i;
		m2.r = m2.i = 0.5 * (root.r - w0.r);
		m2.r *= (t1.r - t2.r);
		m2.i *= (t1.i - t2.i);
		m3.r = root.i * (t3.i + t4.i);
		m3.i = root.i * (-t3.r - t4.r);
		m4.r = (root.i + w0.i) * t4.i;
		m4.i = - (root.i + w0.i) * t4.r;
		m5.r = - (root.i - w0.i) * t3.i;
		m5.i = (root.i - w0.i) * t3.r;
		t3.r = m3.r - m4.r;
		t3.i = m3.i - m4.i;
		t5.r = m3.r + m5.r;
		t5.i = m3.i + m5.i;
		t1.r = m0.r + m1.r;
		t1.i = m0.i + m1.i;
		t2.r = t1.r + m2.r;
		t2.i = t1.i + m2.i;
		t4.r = t1.r - m2.r;
		t4.i = t1.i - m2.i;
		x1[0]   = m0;
		x1[B4].r = t2.r + t3.r;
		x1[B4].i = t2.i + t3.i;
		x1[B3].r = t4.r + t5.r;
		x1[B3].i = t4.i + t5.i;
		x1[B2].r = t4.r - t5.r;
		x1[B2].i = t4.i - t5.i;
		x1[B].r = t2.r - t3.r;
		x1[B].i = t2.i - t3.i;
		break;
	default:
		t = zero;
		for (aB = 0; aB < AB; aB += B) {
			t.r += x0[aB].r;
			t.i += x0[aB].i;
		}
		x1[0] = t;
		for (w0 = root, alphaB = B; alphaB < AB; alphaB += B) {
			x1[alphaB] = x0[0];
			for (w1 = w0, aB = B; aB < AB; aB += B) {
				x1[alphaB].r += w1.r*x0[aB].r - w1.i*x0[aB].i;
				x1[alphaB].i += w1.r*x0[aB].i + w1.i*x0[aB].r;
				t.r = w1.r*w0.r - w1.i*w0.i;
				t.i = w1.r*w0.i + w1.i*w0.r;
				w1 = t;
			}
			t.r = root.r*w0.r - root.i*w0.i;
			t.i = root.r*w0.i + root.i*w0.r;
			w0 = t;
		}
	}
}

static void
FastTransform(Complex x[], int factors[], int A, int C, Complex *pRoot, Complex *pSlowRoot)
{
	int aC, bC, aBC, bAC, BC, AC, B;
	Complex w0, w1, t;

	B = A / *factors;
	BC = B * C;
	A = *factors++;
	AC = A * C;
	for (bC = 0; bC < BC; bC += C)
		SlowTransform(x + bC, z + bC, A, BC, *pSlowRoot);
	for (bC = bAC = 0; bC < BC; bC += C, bAC += AC)
		x[bAC] = z[bC];
	for (w0 = *pRoot, aC = C, aBC = BC; aC < AC; aC += C, aBC += BC) {
		x[aC] = z[aBC];
		for (w1 = w0, bC = C, bAC = AC; bC < BC; bC += C, bAC += AC) {
			t = z[bC + aBC];
			x[aC + bAC].r = w1.r*t.r - w1.i*t.i;
			x[aC + bAC].i = w1.r*t.i + w1.i*t.r;
			t.r = w1.r*w0.r - w1.i*w0.i;
			t.i = w1.r*w0.i + w1.i*w0.r;
			w1 = t;
		}
		t.r = pRoot->r*w0.r - pRoot->i*w0.i;
		t.i = pRoot->r*w0.i + pRoot->i*w0.r;
		w0 = t;
	}
	if (B > 1) {
		pRoot++;
		pSlowRoot++;
		for (aC = 0; aC < AC; aC += C)
			FastTransform(x + aC, factors, B, AC, pRoot, pSlowRoot);
	}
}

static void
FFT(Complex x[], int n, int bShift)
{
	int i, j, k, lastfactor, nHalf;
	int factors[32];
	Complex cTmp, root[32], slowroot[32];
    static float re[1000], im[1000], reOut[1000], imOut[1000];

	k = n;
	AB = n;
    //
    //  factor n
    //
	lastfactor = 0;
	for (j = 0; k > 1; j++) {
		for (i = 0; k % rgnPrime[i]; i++)
			;
		factors[j] = rgnPrime[i];
		root[j].r = cos(TPI / (double)k);
		root[j].i = sin(TPI / (double)k);
		if (factors[j] == lastfactor)
			slowroot[j] = slowroot[j - 1];
		else {
			slowroot[j].r = cos(TPI / (double)factors[j]);
			slowroot[j].i = sin(TPI / (double)factors[j]);
		}
		lastfactor = factors[j];
		k = k / rgnPrime[i];
	}
	factors[j] = 1;
    //
    //  Cooley-Tukey FFT
    //
	FastTransform(x, factors, n, 1, root, slowroot);
    //
    //  Shift zero frequency to the center
    //
    if (bShift) {
        nHalf = n/2;
        for (i = 0; i < nHalf; i++) {
            cTmp = x[i];
            x[i] = x[i + nHalf];
            x[i + nHalf] = cTmp;
        }
    }
}

static void
InverseFFT(Complex x[], int n, int bShift)
{
    int i,nHalf;
    double p, q;
    Complex cTmp;
    //
    //  Shift zero frequency from center
    //
    if (bShift) {
        nHalf = n/2;
        for (i = 0; i < nHalf; i++) {
            cTmp = x[i];
            x[i] = x[i + nHalf];
            x[i + nHalf] = cTmp;
        }
    }
    for (i = 0; i < n; i++)
        x[i].i = -x[i].i;       // negate imaginary
    FFT(x, n, 0);                  // forward FFT
    p = 1.0/double(n);
    q = -p;
    for (i = 0; i < n; i++) {
        x[i].r *= p;
        x[i].i *= q;
    }
}

void
ImageFFT(float *pfPower, float *pfPhase, const float *pfImage, int nHigh, int nWide)
{
    Complex *pc;
    double fScale, fReal, fImag, fMag, fPhase;
    int i, j;

    pc = new Complex[nHigh > nWide ? nHigh : nWide];
    fScale = 1.0/double(nHigh*nWide);
    for (j = 0; j < nHigh; j++) {
        for (i = 0; i < nWide; i++) {
            pc[i].r = pfImage[i + j*nWide];
            pc[i].i = 0.0;
        }
        FFT(pc, nWide, 1);
        for (i = 0; i < nWide; i++) {
            pfPower[i + j*nWide] = float(pc[i].r);
            pfPhase[i + j*nWide] = float(pc[i].i);
        }
    }
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nHigh; j++) {
            pc[j].r = pfPower[i + j*nWide];
            pc[j].i = pfPhase[i + j*nWide];
        }
        FFT(pc, nHigh, 1);
        for (j = 0; j < nHigh; j++) {
            fReal = pc[j].r;
            fImag = pc[j].i;
            fPhase = atan2(fReal, fImag);
            fMag = 0.09*log(fScale*sqrt(fReal*fReal + fImag*fImag) + 1.0e-30) + 1.0;
            pfPower[i + j*nWide] = float(fMag);
            pfPhase[i + j*nWide] = float(fPhase);
        }
    }
    for (i = 0; i < nWide; i++) {
        pc[i].r = pfImage[i];
        pc[i].i = 0.0;
    }
    delete [] pc;
}

void
ImageInverseFFT(float *pfImage, float *pfPower, float *pfPhase, int nHigh, int nWide)
{
    Complex *pc;
    double fScale, fReal, fImag, fMag, fPhase;
    int i, j;

    pc = new Complex[nHigh > nWide ? nHigh : nWide];
    fScale = 1.0/double(nHigh*nWide);
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nHigh; j++) {
            fMag   = pfPower[i + j*nWide];
            fPhase = pfPhase[i + j*nWide];
            fMag = exp((fMag - 1.0)/0.09)/fScale;
            fReal = fMag * sin(fPhase);
            fImag = fMag * cos(fPhase);
            pc[j].r = fReal;
            pc[j].i = fImag;
        }
        InverseFFT(pc, nHigh, 1);
        for (j = 0; j < nHigh; j++) {
            pfPower[i + j*nWide] = float(pc[j].r);  // note: overwrites pfPower, pfPhase
            pfPhase[i + j*nWide] = float(pc[j].i);
        }
    }
    for (j = 0; j < nHigh; j++) {
        for (i = 0; i < nWide; i++) {
            pc[i].r = pfPower[i + j*nWide];
            pc[i].i = pfPhase[i + j*nWide];
        }
        InverseFFT(pc, nWide, 1);
        for (i = 0; i < nWide; i++)
            pfImage[i + j*nWide] = float(pc[i].r);
    }
    delete [] pc;
}