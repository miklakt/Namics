#include "tools.h"
#include "namics.h"
#include "stdio.h"
#include <limits>
#include <cfloat>
#define MAX_BLOCK_SZ 512
#define HALF_MAX_BLOCK_SZ 256

Real* SUM_RESULT;


inline bool safe_mask_compare(const Real& mask_value, const int& query_value) {
	// if mask_value == query_value, i.e., if they're within the numeric limit
    return std::abs(mask_value - static_cast<Real>(query_value)) < std::numeric_limits<Real>::epsilon();
}

void bx(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy, int by1=0, int bz1=0, bool corner=false)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	int jy_by1=jy*by1;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
		if (corner) {
			if (y==0 && z==0) {
				P[i]=P[bx1_jx+i+jy_by1+bz1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1+bz1];
			}
			if (y==My-1 && z==0) {
				P[i]=P[bx1_jx+i-jy_by1+bz1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1+bz1];
			}
			if (y==0 && z==Mz-1) {
				P[i]=P[bx1_jx+i+jy_by1-bz1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1-bz1];
			}
			if (y==My-1 && z==Mz-1) {
				P[i]=P[bx1_jx+i-jy_by1-bz1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1-bz1];
			}
			if (y==0) {
				P[i]=P[bx1_jx+i+jy_by1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1];
			}
			if (y==My-1) {
				P[i]=P[bx1_jx+i-jy_by1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1];
			}
			if (z==0) {
				P[i]=P[bx1_jx+i+bz1];
				P[jx_mmx+i]=P[jx_bxm+i+bz1];
			}
			if (z==Mz-1) {
				P[i]=P[bx1_jx+i-bz1];
				P[jx_mmx+i]=P[jx_bxm+i-bz1];
			}
		}
	}
}
void b_x(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	(void)bxm;
	(void)bx1;
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy, int bz1=0, bool corner=false)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
		if (corner) {
			if (z==0) {
				P[i]=P[jy_by1+i+bz1];
				P[jy_mmy+i]=P[jy_bym+i+bz1];
			}
			if (z==Mz-1) {
				P[i]=P[jy_by1+i-bz1];
				P[jy_mmy+i]=P[jy_bym+i-bz1];
			}
		}
	}
}
void b_y(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	(void)bym;
	(void)by1;
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy, int bx1=0, bool corner=false)   {
	int i;
	int bx1_jx=bx1*jx;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
		if (corner) {
			if (x==0) {
				P[i]=P[i+bz1+bx1_jx];
				P[i+mmz]=P[i+bzm+bx1_jx];
			}
			if (x==Mx-1) {
				P[i]=P[i+bz1-bx1_jx];
				P[i+mmz]=P[i+bzm-bx1_jx];
			}
		}
	}
}
void b_z(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	(void)bzm;
	(void)bz1;
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}
void bx(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy, int by1=0, int bz1=0, bool corner=false)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	int jy_by1=jy*by1;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
		if (corner) {
			if (y==0 && z==0) {
				P[i]=P[bx1_jx+i+jy_by1+bz1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1+bz1];
			}
			if (y==My-1 && z==0) {
				P[i]=P[bx1_jx+i-jy_by1+bz1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1+bz1];
			}
			if (y==0 && z==Mz-1) {
				P[i]=P[bx1_jx+i+jy_by1-bz1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1-bz1];
			}
			if (y==My-1 && z==Mz-1) {
				P[i]=P[bx1_jx+i-jy_by1-bz1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1-bz1];
			}
			if (y==0) {
				P[i]=P[bx1_jx+i+jy_by1];
				P[jx_mmx+i]=P[jx_bxm+i+jy_by1];
			}
			if (y==My-1) {
				P[i]=P[bx1_jx+i-jy_by1];
				P[jx_mmx+i]=P[jx_bxm+i-jy_by1];
			}
			if (z==0) {
				P[i]=P[bx1_jx+i+bz1];
				P[jx_mmx+i]=P[jx_bxm+i+bz1];
			}
			if (z==Mz-1) {
				P[i]=P[bx1_jx+i-bz1];
				P[jx_mmx+i]=P[jx_bxm+i-bz1];
			}
		}
	}
}
void b_x(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	(void)bxm;
	(void)bx1;
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy, int bz1=0, bool corner=false)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
		if (corner) {
			if (z==0) {
				P[i]=P[jy_by1+i+bz1];
				P[jy_mmy+i]=P[jy_bym+i+bz1];
			}
			if (z==Mz-1) {
				P[i]=P[jy_by1+i-bz1];
				P[jy_mmy+i]=P[jy_bym+i-bz1];
			}
		}
	}
}
void b_y(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	(void)bym;
	(void)by1;
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy, int bx1=0, bool corner=false)   {
	int i;
	int bx1_jx=bx1*jx;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
		if (corner) {
			if (x==0) {
				P[i]=P[i+bz1+bx1_jx];
				P[i+mmz]=P[i+bzm+bx1_jx];
			}
			if (x==Mx-1) {
				P[i]=P[i+bz1-bx1_jx];
				P[i+mmz]=P[i+bzm-bx1_jx];
			}
		}
	}
}
void b_z(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	(void)bzm;
	(void)bz1;
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}


#define MAX_ITER 100
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x,y) ((x)>(y)?(x):(y))

Real PYTHAG(Real a, Real b)
{
    Real at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

int svdcmp(Real** a, int m, int n, Real *w, Real** v)
{
    int flag, i, its, j, jj, k, l, nm = 0;
    Real c, f, h, s, x, y, z;
    Real anorm = 0.0, g = 0.0, scale = 0.0;
    Real *rv1;

    if (m < n)
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }

    rv1 = (Real *)malloc((unsigned int) n*sizeof(Real));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs((Real)a[i][k]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a[i][k] = (Real)((Real)a[i][k]/scale);
                    s += ((Real)a[i][k] * (Real)a[i][k]);
                }
                f = (Real)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (Real)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((Real)a[i][k] * (Real)a[j][k]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[j][k] += (Real)(f * (Real)a[i][k]);
                    }
                }
                for (k = i; k < m; k++)
                    a[i][k] = (Real)((Real)a[i][k]*scale);
            }
        }
        w[i] = (Real)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += fabs((Real)a[k][i]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a[k][i] = (Real)((Real)a[k][i]/scale);
                    s += ((Real)a[k][i] * (Real)a[k][i]);
                }
                f = (Real)a[l][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[l][i] = (Real)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (Real)a[k][i] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((Real)a[k][j] * (Real)a[k][i]);
                        for (k = l; k < n; k++)
                            a[k][j] += (Real)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[k][i] = (Real)((Real)a[k][i]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((Real)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    v[j][i] = (Real)(((Real)a[j][i] / (Real)a[l][i]) / g);
                    /* Real division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((Real)a[k][i] * (Real)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (Real)(s * (Real)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (Real)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[j][i] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((Real)a[i][k] * (Real)a[j][k]);
                    f = (s / (Real)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[j][k] += (Real)(f * (Real)a[i][k]);
                }
            }
            for (j = i; j < m; j++)
                a[i][j] = (Real)((Real)a[i][j]*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a[i][j] = 0.0;
        }
        ++a[i][i];
    }

    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        if (a[i][j] != a[i][j])
            throw -3;

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((Real)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (Real)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (Real)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (Real)a[nm][j];
                            z = (Real)a[i][j];
                            a[nm][j] = (Real)(y * c + z * s);
                            a[i][j] = (Real)(z * c - y * s);
                        }
                    }
                }
            }
            z = (Real)w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (Real)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30! iterations \n");
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (Real)w[l];
            nm = k - 1;
            y = (Real)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (Real)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (Real)v[jj][j];
                    z = (Real)v[jj][i];
                    v[jj][j] = (Real)(x * c + z * s);
                    v[jj][i] = (Real)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (Real)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (Real)a[j][jj];
                    z = (Real)a[i][jj];
                    a[j][jj] = (Real)(y * c + z * s);
                    a[i][jj] = (Real)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (Real)x;
        }
    }
    //free((void*) rv1);
    free (rv1);
    return(1);
}
