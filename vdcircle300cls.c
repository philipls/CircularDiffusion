/* ==========================================================================
%  Circular diffusion model, indepedendent Gaussian drift rates. 
%
%      [T, Gt, Theta, Ptheta, Mt] = vdcircle300cls(P, tmax, badix);
%      P = [v1, v2, eta1, eta2, sigma, a]
%  [-pi:pi] closed domain so don't lose last bin in interpolation.
%  Building:
%            mex vdcircle300cls.c -lgsl -lgslcblas  -lm 
% ===========================================================================
*/

#include <mex.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#define kmax 50  /* Maximum number of eigenvalues in dhamana */
#define nw 50    /* Number of steps on circle */
#define sz 300   /* Number of time steps */
#define NP 6     /* Number of input parameters */

const double pi = 3.141592653589793;

void dhamana(double *T, double *Gt0, double *P0, double h, int badix) {
    /* 
      ---------------------------------------------------------------
       First-passage-time density for Bessel process.
       Computes roots of J0 using Gnu GSL library.
      ----------------------------------------------------------------
    */
    double J0k[kmax], J0k_squared[kmax], J1k[kmax]; 
    double a, a2, sigma, sigma2, scaler; 
    int i, k;

    a = P0[0];
    sigma = P0[1];
    sigma2 = sigma * sigma;
    a2 = a * a;
    scaler = sigma2 / a2;
    
    /* Roots of J0k */
    for (k = 0; k < kmax; k++) { 
        J0k[k] = gsl_sf_bessel_zero_J0(k+1);
        /* mexPrintf("k = %6d  J0k[k] = %6.4f\n", k, J0k[k]); - OK */
    }
 
    /* Evaluate Bessel function at the roots */
    for (k = 0; k < kmax; k++) {
        J0k_squared[k] = J0k[k] * J0k[k];
        /* J1k[k] = j1(J0k[k]); */  /* besselj in Gnu library */
        J1k[k] = gsl_sf_bessel_J1(J0k[k]); /* GSL library */
    }
    T[0] = 0;
    Gt0[0] = 0;
    for (i = 1; i < sz; i++) {    
        T[i] = i * h;
        Gt0[i] = 0;
        for (k = 0; k < kmax; k++) {
            Gt0[i] += J0k[k] * exp(-J0k_squared[k] * sigma2 * T[i] / (2.0 * a2)) / J1k[k]; 
        }
        Gt0[i] *= scaler;
        if (i <= badix || Gt0[i] < 0) {
            Gt0[i] = 0;
        }
    }
};

void vdcircle300cls(double *T, double *Gt, double *Theta, double *Ptheta, double *Mt, 
              double *P, double tmax, int badix) {
    /* -----------------------------------------------------------------------------------------------
       Calculate first-passage-time density and response probabilities for circular diffusion process
      ------------------------------------------------------------------------------------------------ */
    
    double Gt0[sz], P0[2];
    double w, two_pi, h, v1, v2, eta1, eta2, sigma, a, sigma2, eta1onsigma2, eta2onsigma2, mt, 
           G11, G12, G21, G22, Girs1, Girs2, tscale, 
           mtscale, totalmass;
    double munorm;
    int i, k;

    two_pi = 2.0 * pi;
    w = 2.0 * pi / nw;
  
    /* Parameters */
    h = tmax / sz; 
    v1 = P[0];
    v2 = P[1];
    eta1 = P[2];
    eta2 = P[3];
    if (eta1 <1e-5) {
       eta1 = 0.01;
    }
    if (eta2 <1e-5) {
       eta2 = 0.01;
    }     
    sigma = P[4];
    a = P[5];
    /*mexPrintf("w= %6.4f h = %6.4f \n", w, h);    */
    munorm = sqrt(v1 * v1 + v2 * v2);

    /* Assume same diffusion in both directions */
    sigma2 = sigma * sigma;  
    eta1onsigma2 = (eta1 * eta1) / sigma2;
    eta2onsigma2 = (eta2 * eta2) / sigma2;
    P0[0] = a;
    P0[1] = sigma;
    /* Density of zero-drift process */
    dhamana(T, Gt0, P0, h, badix); 
    
    /* Response circle (1 x nw + 1) */
    Theta[0] = -pi;
    /* Close the domain */
    for (i = 1; i <= nw; i++) {
         Theta[i] = Theta[i-1] + w;
    }

   /* Joint RT distribution (nw * sz) - make Matlab conformant */
   for (k = 0; k < sz; k++) {
        tscale = sqrt(1/(1 + eta1onsigma2 * T[k])) * sqrt(1/(1 + eta2onsigma2 * T[k]));
        G11 = 2 * eta1 * eta1 * (1 + eta1onsigma2 * T[k]);
        G21 = 2 * eta2 * eta2 * (1 + eta2onsigma2 * T[k]);
        for (i = 0; i < nw; i++) {
            G12 = v1 + a * eta1onsigma2 * cos(Theta[i]);
            G22 = v2 + a * eta2onsigma2 * sin(Theta[i]);
            Girs1 = exp((G12 * G12) / G11 - (v1 * v1) / (eta1 * eta1) / 2);
            Girs2 = exp((G22 * G22) / G21 - (v2 * v2) / (eta2 * eta2) / 2);
            Gt[(nw + 1) * k + i] = tscale * Girs1 * Girs2 * Gt0[k] / two_pi; 
        }
        /* Close the domain */
        Gt[(nw + 1) * k + nw] = Gt[(nw + 1) * k];
    } 
    /* Total mass */
    totalmass = 0;
    for (i = 0; i < nw; i++) {
       for (k = 1; k < sz; k++) {
           totalmass += (Gt[nw * k + i] + Gt[nw * (k - 1) + i]) / 2.0;
       } 
    }
    totalmass *= w * h;
    /*mexPrintf("totalmass = %6.4f\n", totalmass);   */
    /* Integrate joint densities to get means hitting probabilities */
    for (i = 0; i < nw; i++) {
       Ptheta[i] = 0;
       Mt[i] = 0;
       for (k = 1; k < sz; k++) {
           Ptheta[i] += (Gt[(nw + 1) * k + i] + Gt[(nw + 1) * (k - 1) + i]) /2.0;
           Mt[i] += (T[k] * Gt[(nw + 1) * k + i] + T[k - 1] * Gt[(nw + 1) * (k - 1)+ i]) / 2.0; 
       }
       Ptheta[i] *= h / totalmass;
       Mt[i] *= h / Ptheta[i] / totalmass; 
   }
   /* Close the domain but don't double-count the mass */
   Ptheta[nw] = Ptheta[0];
   Mt[nw] = Mt[0];

  /* mt = a * gsl_sf_bessel_I1(a * munorm/(sigma * sigma)) 
            / gsl_sf_bessel_I0(a * munorm/(sigma * sigma)) / munorm; 
   mexPrintf("mt = %6.4f\n", mt); */
} /* vdcircle300 */
   
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 /*
     =======================================================================
     Matlab gateway routine.
     =======================================================================
 */
 

int badix; 

double *T, *Gt, *Theta, *Ptheta, *Mt, *P;
 
double tmax, badi;
 
unsigned n, m;

    if (nrhs != 3) {
         mexErrMsgTxt("dcircle300: Requires 3 input args.");
    } else if (nlhs != 5) {
        mexErrMsgTxt("dcircle300: Requires 5 output args."); 
    }

    /*
      -----------------------------------------------------------------------
      Check all input argument dimensions.
      -----------------------------------------------------------------------   
    */

    /* P */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || !(m * n == NP)) {
        mexPrintf("P is %4d x %4d \n", m, n);
        mexErrMsgTxt("dcircle300: Wrong size P");
    } else {
        P = mxGetPr(prhs[0]);
    }
    /* tmax */
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || !(m * n == 1)) {
        mexErrMsgTxt("dcircle300: tmax must be a scalar");
    } else { 
        tmax = mxGetScalar(prhs[1]);
    }
    if (tmax <= 0.0) {
        mexPrintf("tmax =  %6.2f \n", tmax);
        mexErrMsgTxt("tmax must be positive");
    } 

    /* badi */
    m = mxGetM(prhs[2]);
    n = mxGetN(prhs[2]);
    if (!mxIsDouble(prhs[2]) || !(m * n == 1)) {
        mexErrMsgTxt("dcircle300: badi must be a scalar");
    } else {
        badi = mxGetScalar(prhs[2]); 
        badix = (int)(badi+0.5); 
    }  
 
    /*
      -----------------------------------------------------------------------
      Create output arrays. All of the theta-domain structures have an extra
      bin to duplicate the -pi values in +pi for use in interpolation.
      -----------------------------------------------------------------------
    */
 
    /* T */
    plhs[0] = mxCreateDoubleMatrix(1, sz, mxREAL);
    T = mxGetPr(plhs[0]);
    
    /* Gt */
    plhs[1] = mxCreateDoubleMatrix(nw + 1, sz, mxREAL);
    Gt = mxGetPr(plhs[1]);
    
     /* Theta */
    plhs[2] = mxCreateDoubleMatrix(1, nw + 1, mxREAL);
    Theta = mxGetPr(plhs[2]);


    /* Ptheta */
    plhs[3] = mxCreateDoubleMatrix(1, nw + 1, mxREAL);
    Ptheta = mxGetPr(plhs[3]);

    /* Mt */
    plhs[4] = mxCreateDoubleMatrix(1, nw + 1, mxREAL);
    Mt = mxGetPr(plhs[4]);


    /*
      -----------------------------------------------------------------------
      Run the C-function.
      -----------------------------------------------------------------------
    */

    vdcircle300cls(T, Gt, Theta, Ptheta, Mt, P, tmax, badix);
}


