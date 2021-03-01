#include "stdafx.h"
//
//  Least-Square-Error routines
//  D. P. Mitchell  06/18/2003.
//
#include "Venus.h"
#pragma intrinsic(fabs, sqrt, log)
#define D_MACHINE_EPS   2.22044605e-016
//
//  Vector L2 Norm
//
double
VectorNorm(double *pv, int n, int nSpan)
{
    double f, fMax, fScale, fNorm;
    int i, nExp;

    fMax = 0;
    n = nSpan*(n - 1);
    for (i = n; i >= 0; i -= nSpan) {
        f = fabs(pv[i]);
        if (f > fMax)
            fMax = f;
    }
    //
    //  Scale to avoid overflow.  Use close power of two
    //  to avoid any roundoff error.
    //
    frexp(fMax, &nExp);
    fScale = ldexp(1.0, -nExp);
    fNorm = 0.0;
    for (i = n; i >= 0; i -= nSpan) {
        f = pv[i] * fScale;
        fNorm += f*f;
    }
    return sqrt(fNorm)/fScale;
}
//
//  A = B*C
//  where B: l x m,  C: m x n,  A: l x n
//
static int
MatrixMultiply(double *pmA, const double *pmB, const double *pmC, int l, int m, int n)
{
    int i, j, k;
    const double *pmSaveC;
    double f;

    if (pmA == pmB || pmA == pmC)
        return 0;
    pmSaveC = pmC;
    for (i = 0; i < l; i++) {
        for (j = 0; j < n; j++) {
            f = 0.0;
            for (k = 0; k < m; k++) {
                f += pmB[k] * pmC[j];
                pmC += n;
            }
            *pmA++ = f;
            pmC = pmSaveC;
        }
        pmB += m;
    }
    return 1;
}
//
//  QR Decomposition With Complete Pivoting: AP = QR
//  where: R is m x n.  V is n x m
//  j-th Householder transform: I - pfBeta[j]*Outer(V[m*j], V[m*j])
//
static double s_fMatrixCondition;

static int
HouseholderQR(double *V, double pfBeta[], double *R, int p[], int q[], int m, int n)
{
    int i, j, k, iMax, jMax, nExp;
    double fScale, fNorm, fAvj, fLastNorm, f, fEps, *pfNorms;

    if (m < n)
        return 0;
    pfNorms = V + m*(n - 1);
    for (k = 0; k < n; k++)
        q[k] = p[k] = k;
    fEps = 4.0*sqrt(double(n))*D_MACHINE_EPS;
    for (k = 0; k < n; k++) {
        //
        //  Downdate column norms, recalculate periodically or if bad cancellation
        //
        for (j = k; j < n; j++) {
            if ((k & 7) && fabs(f = R[n*(k - 1) + j]/pfNorms[j]) < 0.999)
                pfNorms[j] *= sqrt(1.0 - f*f);
            else
                pfNorms[j] = VectorNorm(R + n*k + j, m - k, n);
        }
        //
        //  Pivot best column
        //
        jMax = k;
        for (j = k + 1; j < n; j++)
            if (pfNorms[j] > pfNorms[jMax])
                jMax = j;
        p[k] = jMax;
        if (jMax != k) {
            f = pfNorms[k];
            pfNorms[k] = pfNorms[jMax];
            pfNorms[jMax] = f;
            for (i = 0; i < m; i++) {
                f = R[n*i + k];
                R[n*i + k] = R[n*i + jMax];
                R[n*i + jMax] = f;
            }
        }
        //
        //  Pivot best row
        //
        iMax = k;
        for (i = k + 1; i < m; i++)
            if (fabs(R[n*i + k]) > fabs(R[n*iMax + k]))
                iMax = i;
        q[k] = iMax;
        if (iMax != k) {
            for (j = 0; j < n; j++) {
                f = R[n*k + j];
                R[n*k + j] = R[n*iMax + j];
                R[n*iMax + j] = f;
            }
        }
        fNorm = pfNorms[k] = VectorNorm(R + n*k + k, m - k, n);
        //
        //  Rank determination
        //
        if (fNorm == 0.0) {
            s_fMatrixCondition = fNorm;
            return k;
        }
        if (k == 0) {
            fLastNorm = fNorm;
        } else if (fNorm < fLastNorm*fEps) {
            s_fMatrixCondition = fNorm/fLastNorm;
            return k;
        }
        //
        //  Householder: I - fBeta*Outer(V, V)
        //
        frexp(fNorm, &nExp);
        fScale = ldexp(1.0, -nExp);
        for (i = 0; i < k; i++)
            V[m*k + i] = 0.0;
        for (i = k; i < m; i++)
            V[m*k + i] = R[n*i + k] * fScale;   // overwrites pfNorms on last round
        if (R[n*k + k] < 0.0) {
            V[m*k + k] -= fNorm*fScale;
            R[n*k + k] = +fNorm;
        } else {
            V[m*k + k] += fNorm*fScale;
            R[n*k + k] = -fNorm;
        }
        pfBeta[k] = 1.0/fabs(fNorm * fScale * V[m*k + k]);
        for (i = k + 1; i < m; i++)
            R[n*i + k] = 0.0;
        for (j = k + 1; j < n; j++) {
            fAvj = 0.0;
            for (i = k; i < m; i++)
                fAvj += V[m*k + i] * R[n*i + j];
            fAvj *= pfBeta[k];
            for (i = k; i < m; i++)
                R[n*i + j] -= fAvj * V[m*k + i];
        }
    }
    s_fMatrixCondition = fNorm / fLastNorm;
    return k;
}

static int
MultiplyQTranspose(double *M, const double *V, double pfBeta[], int q[], int m, int n, int n2)
{
    int i, j, k, iMax, nResult;
    double fAvj, f, *W;

    nResult = 0;
    W = new double[n2];
    if (m < n || W == 0)
        goto Fail;
    for (k = 0; k < n; k++) {
        iMax = q[k];                            // Row permutation
        for (j = 0; j < n2; j++) {
            f = M[n2*k + j];
            M[n2*k + j] = M[n2*iMax + j];
            M[n2*iMax + j] = f;
        }
        if (pfBeta[k]) {
            for (j = 0; j < n2; j++) {          // Householder reflections
                fAvj = 0.0;
                for (i = k; i < m; i++)
                    fAvj += V[m*k + i] * M[n2*i + j];
                W[j] = fAvj * pfBeta[k];
            }
            for (j = 0; j < n2; j++) {
                for (i = k; i < m; i++)
                    M[n2*i + j] -= V[m*k + i] * W[j]; 
            }
        }

    }
    nResult = 1;
Fail:
    delete [] W;
    return nResult;
}
//
//  Linear Least-Square-Error via QR decomposition
//
int
LeastSquareError(double *A, double *b, double *x, int m, int n)
{
    int i, j, nRank, *p, *q, nResult;
    double f, *V, *pfBeta;

    nResult = 0;
    V = new double[m*n];
    pfBeta = new double[n];
    p = new int[n];
    q = new int[n];
    if (V == 0 || pfBeta == 0 || p == 0 || q == 0)
        goto Fail;
    nRank = HouseholderQR(V, pfBeta, A, p, q, m, n);
    if (nRank != n)
        goto Fail;
    MultiplyQTranspose(b, V, pfBeta, q, m, n, 1);
    for (i = 0; i < n; i++)
        x[i] = b[i];
    x[n - 1] /= A[n*(n - 1) + n - 1];
    for (i = n - 2; i >= 0; --i) {    // Back Substitution
        f = x[i];
        for (j = i + 1; j < n; j++)
            f -= A[n*i + j] * x[j];
        x[i] = f / A[n*i + i];
    }
    for (i = n - 1; i >= 0; --i) {
        f = x[i];
        x[i] = x[p[i]];
        x[p[i]] = f;
    }
    nResult = 1;
Fail:
    delete [] pfBeta;
    delete [] V;
    delete [] p;
    delete [] q;
    return nResult;
}

int
FitPolynomial(double a[], const double x[], const double y[], int nDegree, int nPoints)
{
    int i, j, M, N, nResult;
    double xi, Xj, *rgA, *rgB, *rgX;

    nResult = 0;
    M = nDegree + 1;
    N = nPoints;
    rgA = new double[M*N];
    rgB = new double[N];
    rgX = new double[M];
    for (i = 0; i < nPoints; i++) {
        rgB[i] = y[i];
        xi = x[i];
        Xj = 1.0;
        for (j = 0; j < M; j++) {
            rgA[M*i + j] = Xj;
            Xj *= xi;
        }
    }
    if (LeastSquareError(rgA, rgB, rgX, N, M) == 0)
        goto Fail;
    for (j = 0; j < M; j++)
        a[j] = rgX[j];
    nResult = 1;
Fail:
    delete [] rgX;
    delete [] rgB;
    delete [] rgA;
    return nResult;
}

static double
RandomDouble()
{
    static unsigned n = 123456789;

    n = 1099087573*n + 2654435761;
    return double(n)/4294967296.0;
}

void
TestFit()
{
    double x, y, a[5], b[5], rgX[100], rgY[100];
    int i, nTrials;

    for (nTrials = 0; nTrials < 10; nTrials++) {
        for (i = 0; i < 5; i++)
            b[i] = RandomDouble() - 0.5;
        for (i = 0; i < 100; i++) {
            x = RandomDouble();
            y = b[0] + x*(b[1] + x*(b[2] + x*(b[3] + x*b[4])));
            rgX[i] = x;
            rgY[i] = y + 0.01*(RandomDouble() - 0.5);
        }
        FitPolynomial(a, rgX, rgY, 4, 100);
        printf("Poly: %f %f %f %f %f\n", b[0], b[1], b[2], b[3], b[4]);
        printf("Fit:  %f %f %f %f %f\n", a[0], a[1], a[2], a[3], a[4]);
    }
}


/*
static void
RandomMatrix(double *A, int m, int n, int nRank)
{
    int i, j;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            A[n*i + j] = RandomDouble() - 0.5;
}

static void
GenerateQ(double *Q, double *V, double pfBeta[], int q[], int m, int n)
{
    int i, j, k, iMax;
    double fAvj, f;

    for (i = 0; i < m; i++)                     // Identity
        for (j = 0; j < m; j++)
            Q[m*i + j] = double(i == j);
    for (k = n - 1; k >= 0; --k) {                // Householder transforms
        if (pfBeta[k] == 0.0)
            continue;
        for (j = k; j < m; j++) {
            fAvj = 0.0;
            for (i = k; i < m; i++)
                fAvj += V[m*k + i] * Q[m*i + j];
            fAvj *= pfBeta[k];
            for (i = k; i < m; i++)
                Q[m*i + j] -= fAvj * V[m*k + i]; 
        }
        iMax = q[k];                            // Row permutation
        for (j = 0; j < m; j++) {
            f = Q[m*k + j];
            Q[m*k + j] = Q[m*iMax + j];
            Q[m*iMax + j] = f;
        }
    }
}

static int
MultiplyQ(double *M, double *V, double pfBeta[], int q[], int m, int n, int n2)
{
    int i, j, k, iMax, nResult;
    double fAvj, f, *W;

    nResult = 0;
    W = new double[n2];
    if (m < n || W == 0)
        goto Fail;
    for (k = n - 1; k >= 0; --k) {                // Householder transforms
        if (pfBeta[k]) {
            for (j = 0; j < n2; j++) {
                fAvj = 0.0;
                for (i = k; i < m; i++)
                    fAvj += V[m*k + i] * M[n2*i + j];
                W[j] = fAvj * pfBeta[k];
            }
            for (j = 0; j < n2; j++) {
                for (i = k; i < m; i++)
                    M[n2*i + j] -= V[m*k + i] * W[j]; 
            }
        }
        iMax = q[k];                            // Row permutation
        for (j = 0; j < n2; j++) {
            f = M[n2*k + j];
            M[n2*k + j] = M[n2*iMax + j];
            M[n2*iMax + j] = f;
        }
    }
    nResult = 1;
Fail:
    delete [] W;
    return nResult;
}

void
TestQRCP()
{
    double A[10*10], B[10*10], V[10*10], Q[10*10], R[10*10], QR[10*10], QT[10*10], pfNorms[10];
    double fOrthoError, fQRError, fQError, fQTError, f;
    int i, j, m, n, n2, nTrials, nRank, p[100], q[100];

    printf("    3. QRCP:\n");
    //
    //  Norm, sort increasing: 2.056, 9.365 errors
    //        sort decreasing: 2.104, 9.519
    //        no sorting:      2.073, 9.417
    //
    for (i = 0; i < 10; i++)
        pfNorms[i] = 1.0;
    f = VectorNorm(pfNorms, 10, 1);
    f = fabs(f - sqrt(10.0));
    Assert(f < 1.0e-8);
    //
    //  Test m x n cases, 10 >= m >= n >= 1
    //
    fOrthoError = fQRError = 0.0;
    for (nTrials = 0; nTrials < 1000; nTrials++) {
        for (n = 1; n <= 10; n++) {
            for (m = n; m <= 10; m++) {
                RandomMatrix(A, m, n, n);
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        R[n*i + j] = A[n*i + j];
                nRank = HouseholderQR(V, pfNorms, R, p, q, m, n);
                Assert(nRank == n);
                GenerateQ(Q, V, pfNorms, q, m, n);
                for (i = 0; i < m; i++)
                    for (j = 0; j < m; j++)
                        QT[m*j + i] = Q[m*i + j];
                MatrixMultiply((double *)QR, (double *)Q, (double *)QT, m, m, m);
                for (i = 0; i < m; i++)
                    for (j = 0; j < m; j++)
                        fOrthoError += fabs(QR[m*i + j] - double(i == j));
                MatrixMultiply((double *)QR, (double *)Q, (double *)R, m, m, n);
                for (j = 0; j < n - 1; j++) {
                    for (i = 0; i < m; i++) {
                        f = A[n*i + j];
                        A[n*i + j] = A[n*i + p[j]];
                        A[n*i + p[j]] = f;
                    }
                }
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        fQRError += fabs(QR[n*i + j] - A[n*i + j]);
            }
        }
    }
    fOrthoError /= double(nTrials);
    fQRError /= double(nTrials);
    printf("      a. Orthogonality Error: %g\n", fOrthoError);
    printf("      b. A - QR Error:        %g\n", fQRError);
    //
    //  Test Multiplication by Q transpose
    //
    printf("      e. Test MultiplyQ and MultiplyQTranspose:\n");
    fQError = fQTError = 0.0;
    for (nTrials = 0; nTrials < 100; nTrials++) {
        for (n = 1; n <= 10; n++) {
            for (m = n; m <= 10; m++) {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        R[n*i + j] = RandomDouble() - 0.5;
                nRank = HouseholderQR(V, pfNorms, R, p, q, m, n);
                Assert(nRank == n);
                GenerateQ(Q, V, pfNorms, q, m, n);
                for (i = 0; i < m; i++)
                    for (j = 0; j < m; j++)
                        QT[m*j + i] = Q[m*i + j];
                for (n2 = 1; n2 <= 10; n2++) {
                    for (i = 0; i < m; i++)
                        for (j = 0; j < n2; j++)
                            A[n2*i + j] = RandomDouble() - 0.5;                
                    MatrixMultiply((double *)B, (double *)Q, (double *)A, m, m, n2);
                    MultiplyQ((double *)A, (double *)V, pfNorms, q, m, n, n2);
                    for (i = 0; i < m; i++)
                        for (j = 0; j < n2; j++)
                            fQError += fabs(A[i*n2 + j] - B[i*n2 + j]);
                }
                for (n2 = 1; n2 <= 10; n2++) {
                    for (i = 0; i < m; i++)
                        for (j = 0; j < n2; j++)
                            A[n2*i + j] = RandomDouble() - 0.5;                
                    MatrixMultiply((double *)B, (double *)QT, (double *)A, m, m, n2);
                    MultiplyQTranspose((double *)A, (double *)V, pfNorms, q, m, n, n2);
                    for (i = 0; i < m; i++)
                        for (j = 0; j < n2; j++)
                            fQTError += fabs(A[i*n2 + j] - B[i*n2 + j]);
                }

            }
        }
    }
    printf("        Q Error:    %g\n", fQError/double(nTrials));
    printf("        Q**T Error: %g\n", fQTError/double(nTrials));
    Assert(fQError/double(nTrials) < 1.0e-10);
    Assert(fQTError/double(nTrials) < 1.0e-10);
}

void
TestLSE()
{
    double A[10*20], B[10*20], a[20], b[20], x[10], b2[20], x2[10];
    double fMinQR, fMin, fMinRandom, fRandom, fQR;
    int i, j, m, n, nProblems, nRandom;

    printf("    4. Testing LSE:\n");
    fRandom = fQR = 0.0;
    Assert(LeastSquareError(A, b, x, 3, 4) == 0);
    for (nProblems = 0; nProblems < 1; nProblems++) {
        for (n = 1; n < 10; n++) {                      // n, m == Numerical Recipes M, N
            for (m = n; m < 20; m++) {
                for (i = 0; i < m; i++) {
                    b[i] = a[i] = RandomDouble() - 0.5;
                    for (j = 0; j < n; j++)
                        B[n*i + j] = A[n*i + j] = RandomDouble() - 0.5;
                }
                Assert( LeastSquareError(A, b, x, m, n) );
                MatrixMultiply(b2, (double *)B, x, m, n, 1);
                fMinQR = 0.0;
                for (i = 0; i < m; i++)
                    fMinQR += (b2[i] - a[i]) * (b2[i] - a[i]);
                fQR += fMinQR;
                //
                //  LSE solution should be better than random solutions
                //
                fMinRandom = 1.0e+50;
                for (nRandom = 0; nRandom < 10000; nRandom++) {
                    for (j = 0; j < n; j++)
                        x2[j] = RandomDouble() - 0.5;
                    MatrixMultiply(b2, (double *)B, x2, m, n, 1);
                    fMin = 0.0;
                    for (i = 0; i < m; i++)
                        fMin += (b2[i] - a[i]) * (b2[i] - a[i]);
                    if (fMin < fMinRandom)
                        fMinRandom = fMin;
                }
                fRandom += fMinRandom;
                Assert(fMinRandom >= fMinQR);
            }
        }
    }
    printf("      QR residuals:   %g\n", fQR);
    printf("      RAND residuals: %g\n", fRandom);
}
*/
