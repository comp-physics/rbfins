// MEX implementation of NACA 4-digit airfoil distance function
// Copyright (C) 2024 Spencer H. Bryngelson, Tianci Chu, Oliver T. Schmidt
// 
// This is a high-performance C++ MEX implementation of the dairfoil.m function
// for computing signed distances to NACA 4-digit airfoils.

#include "mex.h"
#include <algorithm>
#include <cmath>
#include <cstring>

using namespace std;

// Inline helper functions
template<class T> inline T sqr(T x) { return x*x; }
inline double min(double a, double b) { return (a < b) ? a : b; }
inline double max(double a, double b) { return (a > b) ? a : b; }

// Function prototypes
double compute_airfoil_distance(double x, double y, double m, double p_camber, double t);
double dairfoil_single(double x, double y, double x_center, double y_center, 
                      double m, double p_camber, double t, double chord, 
                      double cos_a, double sin_a);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check input arguments
    if (nrhs != 6) {
        mexErrMsgTxt("Usage: d = dairfoil(p, x_center, y_center, naca_digits, chord, angle_deg)");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    // Get input dimensions
    mwSize np = mxGetM(prhs[0]);  // Number of points
    if (mxGetN(prhs[0]) != 2) {
        mexErrMsgTxt("Points array p must be N x 2.");
    }
    
    // Get input pointers
    double *p = mxGetPr(prhs[0]);           // Points [N x 2]
    double x_center = mxGetScalar(prhs[1]); // X center
    double y_center = mxGetScalar(prhs[2]); // Y center
    double *naca_digits = mxGetPr(prhs[3]); // NACA digits [4 x 1]
    double chord = mxGetScalar(prhs[4]);    // Chord length
    double angle_deg = mxGetScalar(prhs[5]); // Angle of attack in degrees
    
    // Check NACA digits array
    if (mxGetNumberOfElements(prhs[3]) != 4) {
        mexErrMsgTxt("naca_digits must be a 4-element array.");
    }
    
    // Extract NACA parameters
    double m = naca_digits[0] / 100.0;        // Maximum camber (as fraction of chord)
    double p_camber = naca_digits[1] / 10.0;  // Position of maximum camber (as fraction of chord)
    double t = (naca_digits[2] * 10.0 + naca_digits[3]) / 100.0; // Maximum thickness (as fraction of chord)
    
    // Convert angle to radians and precompute trigonometric functions
    double angle_rad = angle_deg * M_PI / 180.0;
    double cos_a = cos(-angle_rad);  // Negative angle for coordinate transformation
    double sin_a = sin(-angle_rad);
    
    // Create output array
    plhs[0] = mxCreateDoubleMatrix(np, 1, mxREAL);
    double *d = mxGetPr(plhs[0]);
    
    // Process each point
    for (mwSize i = 0; i < np; i++) {
        double x = p[i];           // x-coordinate of point i
        double y = p[i + np];      // y-coordinate of point i (MATLAB column-major storage)
        
        d[i] = dairfoil_single(x, y, x_center, y_center, m, p_camber, t, chord, cos_a, sin_a);
    }
}

double dairfoil_single(double x, double y, double x_center, double y_center, 
                      double m, double p_camber, double t, double chord, 
                      double cos_a, double sin_a)
{
    // Translate points to airfoil coordinate system (leading edge at origin)
    double x_translated = x - x_center;
    double y_translated = y - y_center;
    
    // Rotate points by negative angle of attack to align with airfoil coordinates
    double x_rotated = cos_a * x_translated - sin_a * y_translated;
    double y_rotated = sin_a * x_translated + cos_a * y_translated;
    
    // Normalize by chord length
    double x_norm = x_rotated / chord;
    double y_norm = y_rotated / chord;
    
    double dist;
    
    // Check if point is within airfoil x-bounds
    if (x_norm < 0.0) {
        // Point is before leading edge - distance to leading edge (0, 0)
        dist = sqrt(sqr(x_norm) + sqr(y_norm));
    } else if (x_norm > 1.0) {
        // Point is after trailing edge - distance to trailing edge (1, 0)
        dist = sqrt(sqr(x_norm - 1.0) + sqr(y_norm));
    } else {
        // Point is within chord bounds - compute distance to airfoil surface
        dist = compute_airfoil_distance(x_norm, y_norm, m, p_camber, t);
    }
    
    // Scale back to physical coordinates
    return dist * chord;
}

double compute_airfoil_distance(double x, double y, double m, double p_camber, double t)
{
    // Compute mean camber line and its derivative
    double yc, dyc_dx;
    
    if (x <= p_camber && p_camber > 0.0) {
        // Forward portion of camber line
        double p_camber_sq = sqr(p_camber);
        yc = (m / p_camber_sq) * (2.0 * p_camber * x - sqr(x));
        dyc_dx = (m / p_camber_sq) * (2.0 * p_camber - 2.0 * x);
    } else if (x > p_camber && p_camber > 0.0) {
        // Aft portion of camber line
        double one_minus_p = 1.0 - p_camber;
        double one_minus_p_sq = sqr(one_minus_p);
        yc = (m / one_minus_p_sq) * ((1.0 - 2.0 * p_camber) + 2.0 * p_camber * x - sqr(x));
        dyc_dx = (m / one_minus_p_sq) * (2.0 * p_camber - 2.0 * x);
    } else {
        // Symmetric airfoil (m = 0 or p_camber = 0)
        yc = 0.0;
        dyc_dx = 0.0;
    }
    
    // Compute thickness distribution using NACA equation
    double sqrt_x = sqrt(x);
    double x2 = sqr(x);
    double x3 = x2 * x;
    double x4 = x3 * x;
    
    double yt = (t / 0.2) * (0.2969 * sqrt_x - 0.1260 * x - 0.3516 * x2 + 0.2843 * x3 - 0.1015 * x4);
    
    // Compute angle of camber line
    double theta = atan(dyc_dx);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    
    // Upper and lower surface coordinates
    double xu = x - yt * sin_theta;
    double yu = yc + yt * cos_theta;
    double xl = x + yt * sin_theta;
    double yl = yc - yt * cos_theta;
    
    // Compute distances to upper and lower surfaces
    double dist_upper = sqrt(sqr(x - xu) + sqr(y - yu));
    double dist_lower = sqrt(sqr(x - xl) + sqr(y - yl));
    
    // Minimum distance to surface
    double dist = min(dist_upper, dist_lower);
    
    // Determine sign (inside/outside) using simple y-coordinate comparison
    // Point is inside if it's between upper and lower surfaces at this x-location
    if (y >= yl && y <= yu) {
        // Point is between upper and lower surfaces - inside airfoil
        dist = -dist;
    }
    // If point is outside (y > yu or y < yl), dist remains positive
    
    return dist;
}
