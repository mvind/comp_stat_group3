#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector epan_den(NumericVector x, double h)
{
	int m = 512; // points in grid
	int n= x.size();

	NumericVector out(n);
	for (int i = 0; i < n; ++i)
	{
		out[i] = sqrt(pow(x[i] + h, 2.0));
	}
	return out;
}