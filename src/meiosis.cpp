#include <Rcpp.h>
using namespace Rcpp;


// single random integer from {low, low+1, ..., high}
int random_int(const int low, const int high)
{
    return (int)R::runif((double)low, double(high+1));
}

// L = length of chr in cM
// m = interference paramater (m=0 is no interference)
// p = prop'n chiasmata from no-interference mechanism
// [[Rcpp::export]]
NumericVector sim_crossovers(const double L, const int m, const double p)
{
    RNGScope scope; // to set/reset random number seed

    // chiasma and intermediate points
    int n_points = R::rpois(L/50.0*(double)(m+1)*(1.0-p));
    NumericVector point_locations = runif(n_points, 0.0, L);
    point_locations.sort();

    // which point is the first chiasma?
    int first;
    if(m==0) first = 0;
    else first = random_int(0, m); // random integer from {0, 1, ..., m}

    int n_nichi; // no. chiasma from no interference process
    if(p > 0) {
        n_nichi = R::rpois(L/100.0*p);
    }
    else n_nichi = 0;
    NumericVector nichi_locations = runif(n_nichi, 0.0, L);
    
    // move every (m+1)st point back to front
    int n_chi=0;
    for(int pt_index=first; pt_index < n_points; pt_index += (m+1), n_chi++)
        point_locations[n_chi] = point_locations[pt_index];
    n_chi++;

    // combine interference and no interference chiasma locations
    NumericVector chi_locations(n_chi + n_nichi);
    for(int i=0; i<n_chi; i++)
        chi_locations[i] = point_locations[i];
    for(int i=0; i<n_nichi; i++)
        chi_locations[i+n_chi] = nichi_locations[i];
    chi_locations.sort();

    // thin by 1/2
    int n_xo=0;
    for(int i=0; i<n_chi+n_nichi; i++) {
        if(R::runif(0.0, 1.0) < 0.5) { // flip coin -> chiasma
            chi_locations[n_xo] = chi_locations[i];
            n_xo++;
        }
    }

    NumericVector xo_locations(n_xo);
    for(int i=0; i<n_xo; i++)
        xo_locations[i] = chi_locations[i];

    return xo_locations;
}
