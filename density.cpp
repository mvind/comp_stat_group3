#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector epan_den(NumericVector x, double h, int m = 512)
{
	int n = x.size();
	
	auto max = max_element(begin(x), end(x));
	auto min = min_element(begin(x), end(x));
	
	NumericVector grid(m);
	double l = abs(*max - *min);
		
	for(int i = 1; i < (m - 1); ++i) {
			grid[i] =  *min + (l / m) * i;
	}
	grid[0] = *min;
	grid[m - 1] = *max; // i have no idea
	
	NumericVector y(m);
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			double u = (grid[i] - x[j]) / h;
			if(fabs(u) <= 1.0) {
				y[i] += 3 * (1 - pow(u, 2.0)) / (4);
			}
		}
			y[i] = y[i] / (n * h);
	}
	
	return y;
}