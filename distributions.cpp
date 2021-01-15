using namespace std;
#include "basicmodel.h"

double OffsetGammaCDF (double x, double a, double b, double o) {
	double y = gsl_cdf_gamma_P(x+o+0.5,a,b)-gsl_cdf_gamma_P(x+o-0.5,a,b);
	y=log(y);
	return y;
}

double OffsetGammaCDFFlex (double x, double a, double b, double o) {
	double y=0;
	if (x+o+0.5<0) {
		y=-1e10;
	} else {
		if (x+o-0.5<0) {
			y = gsl_cdf_gamma_P(x+o+0.5,a,b);
		} else {
			y = gsl_cdf_gamma_P(x+o+0.5,a,b)-gsl_cdf_gamma_P(x+o-0.5,a,b);
		}
		y=log(y);
	}
	return y;
}

double LogNormal (double x, double mu, double sigma) {
	double y=0;
	if (x<0) {
		y=-1e10;
	} else {
		if (x>=0.5) {
			y = gsl_cdf_lognormal_P(x+0.5,mu,sigma) - gsl_cdf_lognormal_P(x-0.5,mu,sigma);
		} else {
			y = gsl_cdf_lognormal_P(x+0.5,mu,sigma);
		}
		y=log(y);
	}
	return y;
}

double Poisson (int n, double lambda) {
    double y=1;
    if (lambda!=0||n!=0) {
        y = gsl_ran_poisson_pdf(n,lambda);
        y=log(y);
    }
	return y;
}

double Binomial (int n, double p, int k) {
	double y=1;
	if (n!=0) {
		y = gsl_ran_binomial_pdf(k,p,n);
	}
	//cout << "Binomial " << n << " " << p << " " << k << " " << y << "\n";
	y=log(y);
	return (y);
}
