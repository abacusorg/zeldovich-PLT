#pragma once

#include <cassert>
#include <cstdio>

#include <fmt/base.h>

class SplineFunction {
    // This will be a class to store the points of a spline
    // and return interpolated values
    int n, nmax;
    double *x;
    double *y;
    double *y2;

public:
    SplineFunction() {
        n = nmax = 0;
        x = y = y2 = NULL;
    }
    SplineFunction(int nn) {
        n = nmax = 0;
        x = y = y2 = NULL;
        allocate(nn);
    }
    ~SplineFunction() { deallocate(); }
    void deallocate() {
        delete[] x;
        delete[] y;
        delete[] y2;
        n = nmax = 0;
    }
    int allocate(int nn) {
        if (nmax != 0) deallocate();
        nmax = nn;
        x    = new double[nmax];
        y    = new double[nmax];
        y2   = new double[nmax];
        assert(x != NULL);
        assert(y != NULL);
        assert(y2 != NULL);
        return 0;
    }
    int reallocate(int nn) {
        // Check if we got lucky first
        if (nn == nmax)
            return 0;
        else
            return allocate(nn);
    }
    int size() { return n; };
    // Return N where the tables are 0..N

    void load(double xval, double yval) {
        // Load the next node
        x[n] = xval;
        y[n] = yval;
        n++;
    }
    void load_node(int node, double xval, double yval) {
        // Force a load into a particular node
        assert(node >= 0 && node < nmax);
        x[node] = xval;
        y[node] = yval;
        if (node >= n) n = node + 1;
    }
    void get_node(int node, double *xval, double *yval) {
        assert(node >= 0 && node < n);
        *xval = x[node];
        *yval = y[node];
    }
    void print(FILE *fp) {
        int j;
        for (j = 0; j < n; j++) fmt::print(fp, "{:d}: {:15.10g} {:15.10g}\n", j, x[j], y[j]);
    }

    void sort_arrays() {
        // Sort the (x,y) pairs into increasing x order.
        // This is adapted from Numerical Recipes.
        int i, j, inc;
        double v, w, *a = x - 1, *b = y - 1;
        inc = 1;
        do {
            inc *= 3;
            inc++;
        } while (inc <= n);
        do {
            inc /= 3;
            for (i = inc + 1; i <= n; i++) {
                v = a[i];
                w = b[i];
                j = i;
                while (x[j - inc] > v) {
                    a[j] = a[j - inc];
                    b[j] = b[j - inc];
                    j -= inc;
                    if (j <= inc) break;
                }
                a[j] = v;
                b[j] = w;
            }
        } while (inc > 1);
    }

    void spline() {
        // This will form the interpolation arrays.  Note that we are
        // probably assuming this is monotonic, perhaps monotonically increasing
        // This code is adapted from Numerical Recipes
        // Converted to doubles and to zero-index
        int i, k;
        double yp1 = 1e31, ypn = 1e31;
        double p, qn, sig, un, *u;
        u = new double[n];  // u=vector(1,n-1);
        sort_arrays();
        if (yp1 > 0.99e30)
            y2[0] = u[0] = 0.0;
        else {
            y2[0] = -0.5;
            u[0]  = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
        }
        for (i = 1; i <= n - 2; i++) {
            sig   = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
            p     = sig * y2[i - 1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            u[i]  = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                   - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
        }
        if (ypn > 0.99e30)
            qn = un = 0.0;
        else {
            qn = 0.5;
            un = (3.0 / (x[n - 1] - x[n - 2]))
                 * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
        }
        y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
        for (k = n - 2; k >= 0; k--) y2[k] = y2[k] * y2[k + 1] + u[k];
        delete[] u;
    };

    double val(double v) {
        // This will do the spline interpolation for given v
        // This code is adapted from Numerical Recipes
        // Converted to doubles and zero-index
        int klo, khi, k;
        double h, b, a;
        klo = 0;
        khi = n - 1;
        while (khi - klo > 1) {
            k = (khi + klo) >> 1;
            if (x[k] > v)
                khi = k;
            else
                klo = k;
        }
        h = x[khi] - x[klo];
        assert(h != 0);
        a = (x[khi] - v) / h;
        b = (v - x[klo]) / h;
        return a * y[klo] + b * y[khi]
               + ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * (h * h)
                    / 6.0;
    };
};
