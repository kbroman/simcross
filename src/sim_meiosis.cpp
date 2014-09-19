#include <Rcpp.h>
using namespace Rcpp;

#include "random.h"

//' Simulate crossover locations using the Stahl model
//'
//'
//' Simulate crossover locations on a single meiotic product using the
//' Stahl model. 
//'
//' @details Chiasma locations are a superposition of two
//' processes: a proportion p exhibiting no interference, and a
//' proportion (1-p) following the chi-square model with interference
//' parameter m.  Crossover locations are derived by thinning the
//' chiasma locations with probability 1/2.
//'
//' @param L length of chr in cM
//' @param m Interference paramater (\code{m=0} is no interference)
//' @param p Proportion of chiasmata from no-interference mechanism
//' (\code{p=0} gives pure chi-square model)
//'
//' @return Numeric vector of crossover locations, in cM
//'
//' @keywords datagen
//'
//' @examples
//' x <- sim_crossovers(200, 10, 0)
//' x <- sim_crossovers(200, 10, 0.04)
//'
//' @references
//' Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
//' interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.
//'
//' Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
//' interference as a function of genetic distance. \emph{Genetics}
//' \bold{133}, 681--691.
//'
//' Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
//' of crossover interference using the chi-square model.  \emph{Genetics}
//' \bold{139}, 1045--1056.
//'
//' @export
//' @useDynLib simcross
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
NumericVector sim_crossovers(const double L, const int m=10, const double p=0)
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
