#include <Rcpp.h>
using namespace Rcpp;

// print a vector
void print_vector(NumericVector x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}
