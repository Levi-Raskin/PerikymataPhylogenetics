#include <cmath>
#include <iomanip>
#include <iostream>
#include "Eigen/Dense"
#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#pragma mark - Beta

double Probability::Beta::rv(RandomVariable* rng, double a, double b) {

    double z0 = Probability::Helper::rndGamma( rng, a );
    double z1 = Probability::Helper::rndGamma( rng, b );
    double sum = z0 + z1;
    double x = z0 / sum;
    return x;
}

double Probability::Beta::pdf(double a, double b, double x) {

    double pdf;
    if ( x < 0.0 || 1.0 < x )
        pdf = 0.0;
    else
        pdf = pow(x, (a - 1.0)) * pow((1.0 - x), (b - 1.0)) / Probability::Helper::beta(a, b);
    return pdf;
}

double Probability::Beta::lnPdf(double a, double b, double x) {

    return ( (Probability::Helper::lnGamma(a + b) - Probability::Helper::lnGamma(a) - Probability::Helper::lnGamma(b)) + (a - 1.0) * log(x) + (b - 1.0) * log(1.0 - x) );
}

double Probability::Beta::cdf(double a, double b, double x) {

    double cdf;
    if ( x <= 0.0 )
        cdf = 0.0;
    else if ( x <= 1.0 )
        cdf = Probability::Helper::incompleteBeta(a, b, x);
    else
        cdf = 1.0;
    return cdf;
}

double Probability::Beta::quantile(double alpha, double beta, double x) {

    int        i, nswitches;
    double    curPos, curFraction, increment;
    
    i = nswitches = 0;
    curPos = 0.5;
    bool stopIter = false;
    bool directionUp;
    increment = 0.25;
    curFraction = Probability::Helper::incompleteBeta(alpha, beta, curPos);
    if (curFraction > x)
        {
        directionUp = false;
        }
    else
        {
        directionUp = true;
        }
    
    while (!stopIter)
        {
        curFraction = Probability::Helper::incompleteBeta(alpha, beta, curPos);
        if (curFraction > x && directionUp == false)
            {
            /* continue going down */
            while (curPos - increment <= 0.0)
                {
                increment /= 2.0;
                }
            curPos -= increment;
            }
        else if (curFraction > x && directionUp == true)
            {
            /* switch directions, and go down */
            nswitches++;
            directionUp = false;
            while (curPos - increment <= 0.0)
                {
                increment /= 2.0;
                }
            increment /= 2.0;
            curPos -= increment;
            }
        else if (curFraction < x && directionUp == true)
            {
            /* continue going up */
            while (curPos + increment >= 1.0)
                {
                increment /= 2.0;
                }
            curPos += increment;
            }
        else if (curFraction < x && directionUp == false)
            {
            /* switch directions, and go up */
            nswitches++;
            directionUp = true;
            while (curPos + increment >= 1.0)
                {
                increment /= 2.0;
                }
            increment /= 2.0;
            curPos += increment;
            }
        else
            {
            stopIter = true;
            }
        if (i > 1000 || nswitches > 20)
            stopIter = true;
        i++;
        }
    
    return (curPos);
}

#pragma mark - Chi Square

double Probability::ChiSquare::rv(RandomVariable* rng, double v) {

	/* Cast the degrees of freedom parameter as an integer. We will see
       if there is a decimal remainder later. */
	int n = (int)(v);
	
	double x2;
	if ( (double)(n) == v && n <= 100 )
		{
		/* If the degrees of freedom is an integer and less than 100, we
		   generate our chi-square random variable by generating v
		   standard normal random variables, squaring each, and taking the
		   sum of the squared random variables. */
		x2 = 0.0;
		for (int i=0; i<n; i++)
			{
			double x = Probability::Normal::rv(rng);
			x2 += x * x;
			}
		}
	else
		{
		/* Otherwise, we use the relationship of the chi-square to a gamma
		   (it is a special case of the gamma) to generate the chi-square
		   random variable. */
		x2 = Probability::Gamma::rv(rng, v/2.0, 0.5);
		}
	return x2;
}

double Probability::ChiSquare::pdf(double v, double x) {

	double pdf;
	if ( x < 0.0 )
		{
		pdf = 0.0;
		}
	else
		{
		double b = v / 2.0;
		pdf = exp ( -0.5 * x ) * pow ( x, ( b - 1.0 ) ) / ( pow ( 2.0, b ) * Probability::Helper::gamma( b ) );
		}
	return pdf;
}

double Probability::ChiSquare::cdf(double v, double x) {

	return Probability::Gamma::cdf( v / 2.0, 0.5, x );
}

double Probability::ChiSquare::lnPdf(double v, double x) {

	double b = v / 2.0;
	return ( -(b * log(2.0) + Probability::Helper::lnGamma(b)) - b + (b - 1.0) * log(x) );
}

double Probability::ChiSquare::quantile(double prob, double v) {

	double 		e = 0.5e-6, aa = 0.6931471805, p = prob,
					a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0,
					x = 0.0, b = 0.0;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0)
		return (-1.0);
	double g = Probability::Helper::lnGamma(v/2.0);
	double xx = v/2.0;
	double c = xx - 1.0;
	double ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (v >= -1.24*log(p))
		goto l1;
	if (ch-e < 0)
		return (ch);
	goto l4;
	l1:
		if (v > 0.32)
			goto l3;
		ch = 0.4;
		a = log(1.0-p);
	l2:
		q = ch;
		p1 = 1.0+ch*(4.67+ch);
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0)
			goto l4;
		else
			goto l2;
	l3:
		x = Probability::Helper::pointNormal(p);
		p1 = 0.222222/v;
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;
		p1 = 0.5*ch;
		if ( (t = Probability::Helper::incompleteGamma(p1, xx, g)) < 0.0 )
			{
			std::cerr << "Error in function \"IncompleteGamma" << std::endl;
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));
		b = t/ch;
		a = 0.5*t-b*c;
		double s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		double s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		double s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		double s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		double s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		double s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e)
			goto l4;
	return (ch);
}

#pragma mark - Dirichlet

double Probability::Dirichlet::pdf(const std::vector<double> &a, const std::vector<double> &z) {
	
	int n = (int)a.size();
	double zSum = 0.0;
	for (int i=0; i<n; i++)
		zSum += z[i];

	double tol = 0.0001;
	if ( tol < fabs( zSum - 1.0 ) )
		{
		std::cout << "Fatal error in dirichletPdf" << std::endl;
		exit(1);
		//ui->error("Fatal error in dirichletPdf");
		//throw(MbException(MbException::ERROR));
		}

	double aSum = 0.0;
	for (int i=0; i<n; i++)
		aSum += a[i];

	double aProd = 1.0;
	for (int i=0; i<n; i++)
		aProd *= Probability::Helper::gamma(a[i]);

	double pdf = Probability::Helper::gamma(aSum) / aProd;

	for (int i=0; i<n; i++)
		pdf = pdf * pow( z[i], a[i] - 1.0 );

	return pdf;
}

double Probability::Dirichlet::lnPdf(const std::vector<double> &a, const std::vector<double> &z) {

	int n = (int)a.size(); //!< we assume that a and z have the same size
	double alpha0 = 0.0;
	for (int i=0; i<n; i++)
		alpha0 += a[i];
	double lnP = Probability::Helper::lnGamma(alpha0);
	for (int i=0; i<n; i++)
		lnP -= Probability::Helper::lnGamma(a[i]);
	for (int i=0; i<n; i++)
		lnP += (a[i] - 1.0) * log(z[i]);
	return lnP;
}

void Probability::Dirichlet::rv(RandomVariable* rng, const std::vector<double> &a, std::vector<double> &z) {

    int n = (int)a.size();
    double sum = 0.0;
    for(int i=0; i<n; i++)
        {
        /* z[i] = rndGamma(a[i]) / 1.0; */
        z[i] = Probability::Helper::rndGamma(rng, a[i]);
        sum += z[i];
        }
    for(int i=0; i<n; i++)
        z[i] /= sum;
}

double Probability::Dirichlet::lnPdf(const Eigen::VectorXd &a, const Eigen::VectorXd &z){
    int n = (int)a.size(); //!< we assume that a and z have the same size
    double alpha0 = 0.0;
    for (int i=0; i<n; i++)
        alpha0 += a(i);
    double lnP = Probability::Helper::lnGamma(alpha0);
    for (int i=0; i<n; i++)
        lnP -= Probability::Helper::lnGamma(a(i));
    for (int i=0; i<n; i++)
        lnP += (a(i) - 1.0) * log(z(i));
    return lnP;
}

void Probability::Dirichlet::rv(RandomVariable* rng, const Eigen::VectorXd &a, Eigen::VectorXd &z){
    int n = (int)a.size();
    double sum = 0.0;
    for(int i=0; i<n; i++)
        {
        /* z[i] = rndGamma(a[i]) / 1.0; */
        z(i) = Probability::Helper::rndGamma(rng, a(i));
        sum += z(i);
        }
    for(int i=0; i<n; i++)
        z(i) /= sum;
}

#pragma mark - Exponential

double  Probability::Exponential::pdf(double lambda, double x) {

    return lambda * exp(-lambda*x);
}

double  Probability::Exponential::lnPdf(double lambda, double x) {

    return log(lambda) - lambda * x;
}

double  Probability::Exponential::rv(RandomVariable* rng, double lambda) {

    return -log( rng->uniformRv() ) / lambda;
}

double  Probability::Exponential::cdf(double lambda, double x) {

    return 1.0 - exp(-lambda*x);
}

#pragma mark - Gamma

double  Probability::Gamma::pdf(double alpha, double beta, double x) {

    return (pow(beta, alpha) / Probability::Helper::gamma(alpha)) * pow(x, alpha - 1.0) * exp(-x * beta);
}

double  Probability::Gamma::lnPdf(double alpha, double beta, double x) {

    return alpha * log(beta) - Probability::Helper::lnGamma(alpha) + (alpha-1.0) * log(x) - x * beta; // TODO: Levi, check this
}

double  Probability::Gamma::rv(RandomVariable* rng, double alpha, double beta) {

    return (Probability::Helper::rndGamma(rng, alpha) / beta);
}

double  Probability::Gamma::cdf(double alpha, double beta, double x) {

    return Probability::Helper::incompleteGamma(beta*x, alpha, Probability::Helper::lnGamma(alpha));
}

void Probability::Gamma::discretization(std::vector<double> &catRate, double a, double b, int nCats, bool median) {

	double factor = a / b * nCats;
 
    if (catRate.size() != nCats)
        catRate.resize(nCats);

	if (median)
		{
		/* the median value for each category is used to represent all of the values
		   in that category */
		double interval = 1.0 / (2.0 * nCats);
		for (int i=0; i<nCats; i++)
			catRate[i] = Probability::ChiSquare::quantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
		double t = 0.0;
		for (int i=0; i<nCats; i++)
			t += catRate[i];
		for (int i=0; i<nCats; i++)
			catRate[i] *= factor / t;
		}
	else
		{
		/* the mean value for each category is used to represent all of the values
		   in that category */
		/* calculate the points in the gamma distribution */
		for (int i=0; i<nCats-1; i++)
			catRate[i] = Probability::ChiSquare::quantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
		/* calculate the cumulative values */
		double lnGammaValue = Probability::Helper::lnGamma(a + 1.0);
		for (int i=0; i<nCats-1; i++)
			catRate[i] = Probability::Helper::incompleteGamma(catRate[i] * b, a + 1.0, lnGammaValue);
		catRate[nCats-1] = 1.0;
		/* calculate the relative values and rescale */
		for (int i=nCats-1; i>0; i--)
			{
			catRate[i] -= catRate[i-1];
			catRate[i] *= factor;
			}
		catRate[0] *= factor;
		}
}

#pragma mark - Geometric

int Probability::Geometric::rv(RandomVariable* rng, double p) {

    int n = 0;
    while (rng->uniformRv() < p)
        n++;
    return n;
}

#pragma mark - Normal

double Probability::Normal::rv(RandomVariable* rng) {

	if ( availableNormalRv == false )
		{
		double v1, v2, rsq;
		do
			{
			v1 = 2.0 * rng->uniformRv() - 1.0;
			v2 = 2.0 * rng->uniformRv() - 1.0;
			rsq = v1 * v1 + v2 * v2;
			} while ( rsq >= 1.0 || rsq == 0.0 );
		double fac = sqrt(-2.0 * log(rsq)/rsq);
		extraNormalRv = v1 * fac;
		availableNormalRv = true;
		return v2 * fac;
		}
	else
		{
		availableNormalRv = false;
		return extraNormalRv;
		}
}

double Probability::Normal::rv(RandomVariable* rng, double mu, double sigma) {

	return ( mu + sigma * Probability::Normal::rv(rng) );
}

double Probability::Normal::pdf(double mu, double sigma, double x) {

	double y = ( x - mu ) / sigma;
	return exp( -0.5 * y * y )  / ( sigma * sqrt ( 2.0 * PI ) );
}

double Probability::Normal::lnPdf(double mu, double sigma2, double x) {
    //mu is mean of distribution
    //sigma2 is variance of distribution
    //x is scope of distribution
    return -0.5 * std::log(2.0 * PI * sigma2) - 0.5 * (x - mu) * (x - mu) / (sigma2);
	/*return -0.5 * std::log(2.0 * PI * sigma) - 0.5 * (x - mu) * (x - mu) / (sigma * sigma);*/ //original
}

double Probability::Normal::cdf(double mu, double sigma, double x) {

	double cdf;
	double q;
	double z = (x - mu) / sigma;

	/* |X| <= 1.28 */
	if ( fabs(z) <= 1.28 )
		{
		double a1 = 0.398942280444;
		double a2 = 0.399903438504;
		double a3 = 5.75885480458;
		double a4 = 29.8213557808;
		double a5 = 2.62433121679;
		double a6 = 48.6959930692;
		double a7 = 5.92885724438;
		double y = 0.5 * z * z;
		q = 0.5 - fabs(z) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) );
		}
	else if ( fabs(z) <= 12.7 )
		{
		double b0 = 0.398942280385;
		double b1 = 3.8052E-08;
		double b2 = 1.00000615302;
		double b3 = 3.98064794E-04;
		double b4 = 1.98615381364;
		double b5 = 0.151679116635;
		double b6 = 5.29330324926;
		double b7 = 4.8385912808;
		double b8 = 15.1508972451;
		double b9 = 0.742380924027;
		double b10 = 30.789933034;
		double b11 = 3.99019417011;
		double y = 0.5 * z * z;
		q = exp(-y) * b0 / (fabs(z) - b1 + b2 / (fabs(z) + b3 + b4 / (fabs(z) - b5 + b6 / (fabs(z) + b7 - b8 / (fabs(z) + b9 + b10 / (fabs(z) + b11))))));
		}
	else
		{
		q = 0.0;
		}
	if ( z < 0.0 )
		{
		/* negative x */
		cdf = q;
		}
	else
		{
		/* positive x */
		cdf = 1.0 - q;
		}
	return cdf;
}

double Probability::Normal::invCdf(double mu, double sigma, double scope){
    // from https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.cpp
    double x;
    double x2;

    x2 = invCdfStandardNormal ( scope );

    x = mu + sigma * x2;

    return x;
}

double Probability::Normal::invCdfStandardNormal(double scope){
    // from https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.cpp
    double a[8] = {
    3.3871328727963666080,     1.3314166789178437745E+2,
    1.9715909503065514427E+3,  1.3731693765509461125E+4,
    4.5921953931549871457E+4,  6.7265770927008700853E+4,
    3.3430575583588128105E+4,  2.5090809287301226727E+3 };
  double b[8] = {
    1.0,                       4.2313330701600911252E+1,
    6.8718700749205790830E+2,  5.3941960214247511077E+3,
    2.1213794301586595867E+4,  3.9307895800092710610E+4,
    2.8729085735721942674E+4,  5.2264952788528545610E+3 };
  double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770E-1,
    2.27238449892691845833E-2,  7.74545014278341407640E-4 };
  double const1 = 0.180625;
  double const2 = 1.6;
  double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550E-1,
    1.48103976427480074590E-1,  1.51986665636164571966E-2,
    5.47593808499534494600E-4,  1.05075007164441684324E-9 };
  double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230E-1,
    2.65321895265761230930E-2,  1.24266094738807843860E-3,
    2.71155556874348757815E-5,  2.01033439929228813265E-7 };
  double f[8] = {
    1.0,                        5.99832206555887937690E-1,
    1.36929880922735805310E-1,  1.48753612908506148525E-2,
    7.86869131145613259100E-4,  1.84631831751005468180E-5,
    1.42151175831644588870E-7,  2.04426310338993978564E-15 };
  double q;
  double r;
  double split1 = 0.425;
  double split2 = 5.0;
  double value;

  if ( scope <= 0.0 )
  {
    value = - HUGE_VAL;
    return value;
  }

  if ( 1.0 <= scope)
  {
    value = HUGE_VAL;
    return value;
  }

  q = scope - 0.5;

  if ( fabs ( q ) <= split1 )
  {
    r = const1 - q * q;
    
    //horners polynomial
    double scratch1;
    scratch1 = a[7];
    for (int i = 7 - 1; 0 <= i; i-- )
    {
        scratch1 = scratch1 * r + a[i];
    }
    
    double scratch2;
    scratch2 = b[7];
    for (int i = 7 - 1; 0 <= i; i-- )
    {
        scratch2 = scratch2 * r + b[i];
    }
    
    value = q * scratch1
              / scratch2;
  }
  else
  {
    if ( q < 0.0 )
    {
      r = scope;
    }
    else
    {
      r = 1.0 - scope;
    }

    if ( r <= 0.0 )
    {
      value = HUGE_VAL;
    }
    else
    {
      r = sqrt ( - log ( r ) );

      if ( r <= split2 )
      {
        r = r - const2;
        
        //horners polynomial
        double scratch1;
        scratch1 = c[7];
        for (int i = 7 - 1; 0 <= i; i-- )
        {
            scratch1 = scratch1 * r + c[i];
        }
        
        double scratch2;
        scratch2 = d[7];
        for (int i = 7 - 1; 0 <= i; i-- )
        {
            scratch2 = scratch2 * r + d[i];
        }
        value = scratch1/scratch2;
       }
       else
       {
         r = r - split2;
         
        //horners polynomial
        double scratch1;
        scratch1 = e[7];
        for (int i = 7 - 1; 0 <= i; i-- )
        {
            scratch1 = scratch1 * r + e[i];
        }
        
        double scratch2;
        scratch2 = f[7];
        for (int i = 7 - 1; 0 <= i; i-- )
        {
            scratch2 = scratch2 * r + f[i];
        }
        value = scratch1/scratch2;
        }
    }

    if ( q < 0.0 )
    {
      value = - value;
    }

  }

  return value;
}

double Probability::Normal::quantile(double mu, double sigma, double p) {
	
	double z = Probability::Helper::pointNormal(p);
	double x = z * sigma + mu;
	return x;
}

#pragma mark - Truncated Normal

double Probability::TruncatedNormal::lnPdf(double x, double mu, double sigma, double lower, double upper){
    //Per: https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.html
    double alpha;
    double alpha_cdf;
    double beta;
    double beta_cdf;
    double pdf;
    double xi;
    double xi_lnPdf;

    if ( x < lower )
    {
        pdf = 0.0;
    }
    else if ( x <= upper )
    {
    alpha = ( lower - mu ) / sigma;
    beta = ( upper - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = Probability::Normal::cdf(0, 1, alpha);
    beta_cdf = Probability::Normal::cdf(0, 1, beta);
    xi_lnPdf = Probability::Normal::lnPdf(0, 1, xi);

    pdf = xi_lnPdf - std::log( beta_cdf - alpha_cdf ) - std::log(sigma);
    }
    else
    {
    pdf = 0.0;
    }

    return pdf;
}

double Probability::TruncatedNormal::rv(RandomVariable* rng, double mu, double sigma, double lower, double upper){
    //Per: https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.html
    double alpha;
    double alpha_cdf;
    double beta;
    double beta_cdf;
    double u;
    double x;
    double xi;
    double xi_cdf;

    alpha = ( lower - mu ) / sigma;
    beta = ( upper - mu ) / sigma;

    alpha_cdf = Probability::Normal::cdf(0, 1, alpha);
    beta_cdf = Probability::Normal::cdf(0, 1, beta);

    u = rng->uniformRv();
    xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
    xi = Probability::Normal::invCdf(0, 1, xi_cdf);

    x = mu + sigma * xi;

    return x;
}

#pragma mark - Uniform

double  Probability::Uniform::pdf(double low, double high, double x) {

    return 1.0 / (high - low);
}

double  Probability::Uniform::lnPdf(double low, double high, double x) {

    return -log(high - low);
}

double  Probability::Uniform::rv(RandomVariable* rng, double low, double high) {

    return low + rng->uniformRv() * (high - low);
}

double  Probability::Uniform::cdf(double low, double high, double x) {

    if (x < low)
        return 0.0;
    else if (x > high)
        return 1.0;
    else
        return 1.0 - ((high - x) / (high - low));
}

#pragma mark - Multivariate T distribution

double Probability::MultivariateT::lnPdf(Eigen::MatrixXd* dat, Eigen::MatrixXd* u, Eigen::MatrixXd* scale, double dof){
    if(dof <= 0)
        Msg::error("dof needs to be greater than 0");
    if(scale->rows() != scale->cols())
        Msg::error("expected square scale matrix");
        
    //we need u and x to be p x 1
    int p = (int)scale->cols();
    
    Eigen::MatrixXd mu;
    Eigen::MatrixXd x;
    if(dat->rows() == 1){
        x = dat->transpose();
    }else if(dat->rows() == p){
        x = *dat;
    }else{
        Msg::error("expected dat to be 1xp or px1 and it was neither");
    }
    if(u->rows() == 1){
        mu = u->transpose();
    }else if(u->rows() == p){
        mu = *dat;
    }else{
        Msg::error("expected u to be 1xp or px1 and it was neither");
    }
    
    double res = std::lgamma((dof + p)/2);
    res -= std::lgamma(dof/2);
    res -= 0.5 * p * std::log(dof);
    res -= 0.5 * p * std::log(PI);
    res -= 0.5 * std::log(scale->determinant());
    res += -0.5 * (dof + p) * std::log(1 + (1/dof * (x-mu).transpose() *  scale->inverse() * (x-mu))(0,0));
    
    return res;
}

#pragma mark - Multivariate Normal
std::vector<double> Probability::MultivariateNormal::rv(RandomVariable* rng, std::vector<double> u, Eigen::MatrixXd* VCV){
    //u + lltVCV * z where z is n x 1 vector of standard normal draws
    Eigen::LDLT<Eigen::MatrixXd> check(*VCV);
    if(check.info() == Eigen::NumericalIssue)
        throw std::invalid_argument("vcv apparently not pos def");
        
    Eigen::LLT<Eigen::MatrixXd> cholDec(*VCV);
    Eigen::MatrixXd lltVCV = cholDec.matrixL();
    
    Eigen::MatrixXd z = Eigen::MatrixXd::Zero(u.size(), 1);
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(u.size(), 1);
    
    for(int i = 0; i < z.rows(); i++){
        z(i, 0) = Probability::Normal::rv(rng);
        U(i, 0) = u[i];
    }
    
    Eigen::MatrixXd res = U + lltVCV * z; // expecting n x 1 matrix here
    
    //converting back
    std::vector<double> resVec;
    resVec.clear();
    for(int i = 0; i < res.rows(); i++)
        resVec.push_back(res(i, 0));
    return resVec;
}

Eigen::VectorXd Probability::MultivariateNormal::rv(RandomVariable* rng, Eigen::VectorXd u, Eigen::MatrixXd* VCV){
    //u + lltVCV * z where z is n x 1 vector of standard normal draws
    Eigen::LLT<Eigen::MatrixXd> cholDec(*VCV);
    Eigen::MatrixXd lltVCV = cholDec.matrixL();
    Eigen::VectorXd z = Eigen::VectorXd::Zero(u.size());
    for(int i = 0; i < z.rows(); i++)
        z(i) = Probability::Normal::rv(rng);
    Eigen::VectorXd res = u + lltVCV * z;
    
    return res;
}

Eigen::VectorXd Probability::MultivariateNormal::rv(RandomVariable* rng, const Eigen::VectorXd& u, const Eigen::MatrixXd& VCV){
    //u + lltVCV * z where z is n x 1 vector of standard normal draws
    Eigen::LLT<Eigen::MatrixXd> cholDec(VCV);
    Eigen::MatrixXd lltVCV = cholDec.matrixL();
    Eigen::VectorXd z = Eigen::VectorXd::Zero(u.size());
    for(int i = 0; i < z.rows(); i++)
        z(i) = Probability::Normal::rv(rng);
    Eigen::VectorXd res = u + lltVCV * z;
    
    return res;
}

void Probability::MultivariateNormal::rv(RandomVariable* rng, Eigen::VectorXd& out, const Eigen::VectorXd& means, const Eigen::MatrixXd& vcvCholLower){
    assert(&out != &means && "out and means must be distinct vectors");
    for(int i = 0; i < out.rows(); i++)
        out(i) = Probability::Normal::rv(rng);
    out = vcvCholLower * out;
    out += means;
}

double Probability::MultivariateNormal::lnPdf(Eigen::MatrixXd* vec, Eigen::MatrixXd* distribVec, Eigen::MatrixXd* VCV){
    //Checked good LYR Aug 28, 2025
    Eigen::VectorXd scope = *vec;
    Eigen::VectorXd mean = *distribVec;
    Eigen::MatrixXd sigma = *VCV;
    double prob = std::log(1 / std::sqrt( (2*PI* sigma).determinant() )) + (-0.5 * (scope-mean).transpose() * sigma.inverse() * (scope - mean) );
    return prob;
    
    /*
    //need: vec to be nx1, distribVec to be nx1, vcv to be nxn; vec is scope, distribVec/VCV are dist. params
    Eigen::MatrixXd x = *vec; //need to be nx1
    Eigen::MatrixXd m = *distribVec; //need to be nx1
    Eigen::MatrixXd sigma = *VCV;
//    Eigen::LLT<Eigen::MatrixXd> lltOfV(sigma);
    Eigen::LDLT<Eigen::MatrixXd> lltOfV(sigma); //more numerically stable
    if(lltOfV.info() == Eigen::NumericalIssue){
        printEigen(sigma);
        Msg::error("non posiative definite VCV");
    }
    if(sigma.cols() != sigma.rows())
        Msg::error("MVN lnPdf expecting sigma to be square matrix");
    if(x.cols() > x.rows())
        x.transposeInPlace(); //if 1xn, transpose to nx1
    if(m.cols() > m.rows())
        m.transposeInPlace(); //if 1xn, transpose to nx1
    if(x.cols() != 1)
        Msg::error("Scope of mvn is expecting nx1");
    if(m.cols() != 1)
        Msg::error("Mean of mvn is expecting nx1");
    Eigen::MatrixXd sigInv = sigma.inverse();
//    Eigen::MatrixXd sigInv = lltOfV.solve(Eigen::MatrixXd::Identity(sigma.rows(), sigma.cols()));
    Eigen::MatrixXd scatter = x-m;
    
    double logP = -0.5 * (scatter.transpose() * sigInv * scatter)(0,0);
//    logP -= (0.5 * (sigma.cols()*std::log(2 * PI) + std::log(sigma.determinant())));
    logP -= (0.5 * (sigma.cols()*std::log(2 * PI) + lltOfV.vectorD().array().log().sum()));
    return logP;
    */
}

double Probability::MultivariateNormal::lnPdf(Eigen::VectorXd* vec, Eigen::VectorXd* u, Eigen::MatrixXd* VCV){
    //Checked good LYR Aug 28, 2025
    Eigen::VectorXd scope = *vec;
    Eigen::VectorXd mean = *u;
    Eigen::MatrixXd sigma = *VCV;
    double prob = std::log(1 / std::sqrt( (2*PI* sigma).determinant() )) + (-0.5 * (scope-mean).transpose() * sigma.inverse() * (scope - mean) );
    return prob;
    
    /*
    //need: vec to be nx1, distribVec to be nx1, vcv to be nxn; vec is scope, distribVec/VCV are dist. params
    Eigen::MatrixXd x = *vec; //need to be nx1
    Eigen::MatrixXd m = *distribVec; //need to be nx1
    Eigen::MatrixXd sigma = *VCV;
//    Eigen::LLT<Eigen::MatrixXd> lltOfV(sigma);
    Eigen::LDLT<Eigen::MatrixXd> lltOfV(sigma); //more numerically stable
    if(lltOfV.info() == Eigen::NumericalIssue){
        printEigen(sigma);
        Msg::error("non posiative definite VCV");
    }
    if(sigma.cols() != sigma.rows())
        Msg::error("MVN lnPdf expecting sigma to be square matrix");
    if(x.cols() > x.rows())
        x.transposeInPlace(); //if 1xn, transpose to nx1
    if(m.cols() > m.rows())
        m.transposeInPlace(); //if 1xn, transpose to nx1
    if(x.cols() != 1)
        Msg::error("Scope of mvn is expecting nx1");
    if(m.cols() != 1)
        Msg::error("Mean of mvn is expecting nx1");
    Eigen::MatrixXd sigInv = sigma.inverse();
//    Eigen::MatrixXd sigInv = lltOfV.solve(Eigen::MatrixXd::Identity(sigma.rows(), sigma.cols()));
    Eigen::MatrixXd scatter = x-m;
    
    double logP = -0.5 * (scatter.transpose() * sigInv * scatter)(0,0);
//    logP -= (0.5 * (sigma.cols()*std::log(2 * PI) + std::log(sigma.determinant())));
    logP -= (0.5 * (sigma.cols()*std::log(2 * PI) + lltOfV.vectorD().array().log().sum()));
    return logP;
    */
}

double Probability::MultivariateNormal::lnPrecisionMatrixPdf(Eigen::VectorXd* vec, Eigen::VectorXd* u, Eigen::MatrixXd* VCVInv){
    Eigen::VectorXd scope = *vec;
    Eigen::VectorXd mean = *u;
    Eigen::MatrixXd sigma = *VCVInv;
    double prob = std::log(sigma.determinant() / std::sqrt( 2*PI )) + (-0.5 * (scope-mean).transpose() * sigma * (scope - mean) );
    return prob;
    
    /*
    //need: vec to be nx1, distribVec to be nx1, vcv to be nxn; vec is scope, distribVec/VCV are dist. params
    Eigen::MatrixXd x = *vec; //need to be nx1
    Eigen::MatrixXd m = *distribVec; //need to be nx1
    Eigen::MatrixXd sigma = *VCV;
//    Eigen::LLT<Eigen::MatrixXd> lltOfV(sigma);
    Eigen::LDLT<Eigen::MatrixXd> lltOfV(sigma); //more numerically stable
    if(lltOfV.info() == Eigen::NumericalIssue){
        printEigen(sigma);
        Msg::error("non posiative definite VCV");
    }
    if(sigma.cols() != sigma.rows())
        Msg::error("MVN lnPdf expecting sigma to be square matrix");
    if(x.cols() > x.rows())
        x.transposeInPlace(); //if 1xn, transpose to nx1
    if(m.cols() > m.rows())
        m.transposeInPlace(); //if 1xn, transpose to nx1
    if(x.cols() != 1)
        Msg::error("Scope of mvn is expecting nx1");
    if(m.cols() != 1)
        Msg::error("Mean of mvn is expecting nx1");
    Eigen::MatrixXd sigInv = sigma.inverse();
//    Eigen::MatrixXd sigInv = lltOfV.solve(Eigen::MatrixXd::Identity(sigma.rows(), sigma.cols()));
    Eigen::MatrixXd scatter = x-m;
    
    double logP = -0.5 * (scatter.transpose() * sigInv * scatter)(0,0);
//    logP -= (0.5 * (sigma.cols()*std::log(2 * PI) + std::log(sigma.determinant())));
    logP -= (0.5 * (sigma.cols()*std::log(2 * PI) + lltOfV.vectorD().array().log().sum()));
    return logP;
    */
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Probability::MultivariateNormal::productMVN(Eigen::MatrixXd& mean0, Eigen::MatrixXd& var0,Eigen::MatrixXd& mean1, Eigen::MatrixXd& var1){
    Eigen::MatrixXd vcvSumInv = (var0+var1).inverse();
    Eigen::MatrixXd meanProd = var1*vcvSumInv*mean0 + var0*vcvSumInv*mean1;
    Eigen::MatrixXd varProd = var0*vcvSumInv*var1;
//    Eigen::LDLT<Eigen::MatrixXd> ldlt(var0+var1);
//    Eigen::MatrixXd meanProd = var1*ldlt.solve(mean0) + var0*ldlt.solve(mean1);
//    Eigen::MatrixXd varProd = var0*ldlt.solve(var1);
    
    return {meanProd, varProd};
}


#pragma mark - Wishart
  
Eigen::MatrixXd Probability::Wishart::rv(RandomVariable* rng, Eigen::MatrixXd* V, double n){
//    Msg::warning("possibly bad???");
    int dim = (int)V->rows();
    if(n < dim)
        Msg::error("DOF must be greater or equal to dimension of matrix");

    int retryCnt = 0;
    retry:
    retryCnt++;
    
    Eigen::LLT<Eigen::MatrixXd> cholDec(*V);
    Eigen::MatrixXd L = cholDec.matrixL();
    Eigen::MatrixXd Lt = cholDec.matrixL().transpose();
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim, dim);
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            if(i==j)
                A(i,j) = std::sqrt(Probability::ChiSquare::rv(rng, (n - i + 1)));
            else if(i > j)
                A(i,j) = Probability::Normal::rv(rng); //standard normal
        }
    }
//    printEigen(A);
    Eigen::MatrixXd LALAT = L * A * A.transpose() * Lt;
    
    if(retryCnt > 100)
        Msg::error("can't manage to draw from wishart matrix with determinant that is greater than 0");
    if (LALAT.determinant() <= 0)
        goto retry;
    
    return LALAT;
}

double Probability::Wishart::lnPdf(Eigen::MatrixXd* support, Eigen::MatrixXd* scale, double dof){
    //lYR checked 5/28/2025
    Eigen::MatrixXd X = *support;
    Eigen::MatrixXd V = *scale;
    double          n = dof;
    int p = (int)V.cols();
    if(n < p)
        Msg::error("DOF must be greater or equal to dimension of matrix");
    double detX = X.determinant();
    double detV = V.determinant();
//    if (detX <= 0 || detV <= 0)
//        Msg::error("Support and scale matrices must be positive definite");
    
    Eigen::LDLT<Eigen::MatrixXd> ldltV(V);
    Eigen::LDLT<Eigen::MatrixXd> ldltX(X);
    if (ldltX.info() == Eigen::NumericalIssue)
        Msg::error("X is not pos def?");
    if (ldltV.info() == Eigen::NumericalIssue)
        Msg::error("V is not pos def?");

    double logP = (0.5*(n-p-1) * std::log(detX)) - 0.5 * (V.inverse() * X).trace();
//    double logP = (0.5*(n-p-1) * std::log(detX)) - 0.5 * (ldltV.solve(X)).trace();
    logP -= (0.5*n*p*std::log(2) + 0.5*n*std::log(detV));
    
    double nOver2 = n*0.5;
    double multiGammaRes = p * (p - 1) / 4.0 * std::log(PI);
    for (int j = 1; j <= p; ++j) {
        multiGammaRes += std::lgamma(nOver2 - (j - 1) / 2.0);
    }
    logP -= (multiGammaRes);

    return logP;
}

#pragma mark - Inverse Wishart
  
Eigen::MatrixXd Probability::InverseWishart::rv(RandomVariable* rng, Eigen::MatrixXd* psi, double n){
    int dim = psi->rows();
    
    // Generate Wishart(psi^-1, n)
    Eigen::MatrixXd psi_inv = psi->inverse();
    Eigen::LLT<Eigen::MatrixXd> llt(psi_inv);
    Eigen::MatrixXd L = llt.matrixL();
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim, dim);
    for(int i = 0; i < dim; ++i) {
        A(i,i) = std::sqrt(Probability::ChiSquare::rv(rng, n - i));
        for(int j = 0; j < i; ++j) {
            A(i,j) = Probability::Normal::rv(rng);
        }
    }
    
    Eigen::MatrixXd LA = L * A;
    Eigen::MatrixXd W = LA * LA.transpose();
    return W.inverse();

}

void Probability::InverseWishart::rv(RandomVariable* rng, Eigen::MatrixXd& out, const Eigen::MatrixXd& psiInvCholLower, double n){
    out.setZero();
    for(int i = 0; i < psiInvCholLower.rows(); ++i) {
        out(i,i) = std::sqrt(Probability::ChiSquare::rv(rng, n - i));
        for(int j = 0; j < i; ++j) {
            out(i,j) = Probability::Normal::rv(rng);
        }
    }
    out = psiInvCholLower * out;
    out = (out * out.transpose()).eval();
    out = out.inverse();
}

double Probability::InverseWishart::lnPdf(Eigen::MatrixXd* support, Eigen::MatrixXd* scale, double dof){
    int p = support->cols();
    double res = 0.0;
    
    if(dof <= p - 1)
        Msg::error("dof needs to be greater than p - 1");
    
    thread_local static double log2 = std::log(2);
    
    double vHalf = dof/2;

    res += vHalf * std::log(scale->determinant());
    
    res -= vHalf * p * log2;
    
    thread_local static double logPi = std::log(PI);
    double multiGammaRes = p * (p - 1) / 4.0 * logPi;
    for (int j = 1; j <= p; ++j) {
        multiGammaRes += std::lgamma(vHalf - (1-j) / 2.0);
    }
    res -= multiGammaRes;
    
    res += -(dof+p+1)/2 * std::log(support->determinant());
    
    res += -0.5 * (*(scale) * support->inverse()).trace();
    
    return res;
}

Eigen::MatrixXd Probability::InverseWishart::rv(RandomVariable* rng, const Eigen::MatrixXd& psi, double n){
    int dim = psi.rows();
    
    // Generate Wishart(psi^-1, n)
    Eigen::MatrixXd psi_inv = psi.inverse();
    Eigen::LLT<Eigen::MatrixXd> llt(psi_inv);
    Eigen::MatrixXd L = llt.matrixL();
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim, dim);
    for(int i = 0; i < dim; ++i) {
        A(i,i) = std::sqrt(Probability::ChiSquare::rv(rng, n - i));
        for(int j = 0; j < i; ++j) {
            A(i,j) = Probability::Normal::rv(rng);
        }
    }
    
    Eigen::MatrixXd LA = L * A;
    Eigen::MatrixXd W = LA * LA.transpose();  // This is Wishart(psi^-1, n)
    
    // Return inverse for Inverse-Wishart(psi, n)
    return W.inverse();
//    int dim = psi.rows();
//    if(psi.cols() != dim)
//        throw std::runtime_error("Scale matrix psi must be square");
//
//    if(n <= dim - 1)
//        throw std::runtime_error("Degrees of freedom n must be > dim - 1");
//
//    // Cholesky decomposition: psi = L * L^T
//    Eigen::LLT<Eigen::MatrixXd> llt(psi);
//    Eigen::MatrixXd L = llt.matrixL();
//
//    // Build generalized Bartlett matrix A
//    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim, dim);
//    for(int i = 0; i < dim; ++i) {
//        // Diagonal: sqrt of chi-square with real df = n - i
//        double df = n - i;
//        if(df <= 0.0) throw std::runtime_error("Invalid degrees of freedom on diagonal");
//        A(i,i) = std::sqrt(Probability::ChiSquare::rv(rng, df));
//
//        // Lower triangle: standard normal
//        for(int j = 0; j < i; ++j) {
//            A(i,j) = Probability::Normal::rv(rng);
//        }
//    }
//
//    // QR decomposition of A to get upper-triangular R
//    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
//    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
//
//    // Solve R^T * T^T = L^T  => T = ...
//    Eigen::MatrixXd Tt = R.transpose().triangularView<Eigen::Upper>().solve(L.transpose());
//    Eigen::MatrixXd T = Tt.transpose();
//
//    // Return Inverse Wishart sample
//    return T * T.transpose();
}

double Probability::InverseWishart::lnPdf(const Eigen::MatrixXd& x, const Eigen::MatrixXd& psi, double v){
    int p = x.cols();
    double res = 0.0;
    
    if(v <= p - 1)
        Msg::error("dof needs to be greater than p - 1");
    
    thread_local static double log2 = std::log(2);
    
    double vHalf = v/2;

    res += vHalf * std::log(psi.determinant());
    
    res -= vHalf * p * log2;
    
    thread_local static double logPi = std::log(PI);
    double multiGammaRes = p * (p - 1) / 4.0 * logPi;
    for (int j = 1; j <= p; ++j) {
        multiGammaRes += std::lgamma(vHalf - (1-j) / 2.0);
    }
    res -= multiGammaRes;
    
    res += -(v+p+1)/2 * std::log(x.determinant());
    
    res += -0.5 * (psi * x.inverse()).trace();
    
    return res;
}

#pragma mark - LKJ

double Probability::LKJ::lnPdf(Eigen::MatrixXd* support, double shape){
    Eigen::MatrixXd R = *support;
    int dim = R.rows();
    
    double lnPro = 0.0;
    lnPro += (shape - 1) * std::log(R.determinant());
    
    //calcualting normalizing constant C
    double scratch = 0.0;
    for(int k = 0; k < (dim - 1); k++)
        scratch += (2*shape - 2 + dim - k)*(dim - k);
    lnPro += scratch * std::log(2);
    
    for(int k = 0; k < (dim - 1); k++)
        lnPro += (dim - k) * std::log(Probability::Helper::beta(shape + (dim - k - 1)/2, shape + (dim-k-1)/2));
    return lnPro;
}

#pragma mark - Helper Functions

double Probability::Helper::beta(double a, double b) {

    return ( exp(lnGamma(a) + lnGamma(b) - lnGamma(a + b)) );
}

double Probability::Helper::chebyshevEval(double x, const double *a, const int n) {

    double b0, b1, b2, twox;
    int i;
    
    if (n < 1 || n > 1000)
    {
        Msg::error("Cannot compute chebyshev function");
    }
    
    if (x < -1.1 || x > 1.1)
    {
        Msg::error("Cannot compute chebyshev function");
    }
    
    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++)
    {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

/*
 * This function returns the round off unit for floating point arithmetic.
 * The returned value is a number, r, which is a power of 2 with the property
 * that, to the precision of the computer's arithmetic, 1 < 1 + r, but
 * 1 = 1 + r / 2. This function comes from John Burkardt.
 */
double Probability::Helper::epsilon(void) {

    double r = 1.0;
    while ( 1.0 < (double)(1.0 + r)  )
        r = r / 2.0;
    return 2.0 * r;
}

/*!
 * This function calculates the gamma function for real x.
 */
double Probability::Helper::gamma(double x) {

    double c[7] = { -1.910444077728E-03,
                    8.4171387781295E-04,
                    -5.952379913043012E-04,
                    7.93650793500350248E-04,
                    -2.777777777777681622553E-03,
                    8.333333333333333331554247E-02,
                    5.7083835261E-03 };
    double fact;
    int i;
    int n;
    double p[8] = { -1.71618513886549492533811,
                    2.47656508055759199108314E+01,
                    -3.79804256470945635097577E+02,
                    6.29331155312818442661052E+02,
                    8.66966202790413211295064E+02,
                    -3.14512729688483675254357E+04,
                    -3.61444134186911729807069E+04,
                    6.64561438202405440627855E+04 };
    bool parity;
    double q[8] = { -3.08402300119738975254353E+01,
                    3.15350626979604161529144E+02,
                    -1.01515636749021914166146E+03,
                    -3.10777167157231109440444E+03,
                    2.25381184209801510330112E+04,
                    4.75584627752788110767815E+03,
                    -1.34659959864969306392456E+05,
                    -1.15132259675553483497211E+05 };
    double sqrtpi = 0.9189385332046727417803297;
    double sum2;
    double value;
    double xbig = 35.040;
    double xden;
    double xminin = 1.18E-38;
    double xnum;
    double y;
    double y1;
    double ysq;
    double z;

    parity = false;
    fact = 1.0;
    n = 0;
    y = x;

    if ( y <= 0.0 )
        {
        /* argument negative */
        y = -x;
        y1 = ( double ) ( ( int ) ( y ) );
        value = y - y1;

        if ( value != 0.0 )
            {
            if ( y1 != ( double ) ( ( int ) ( y1 * 0.5 ) ) * 2.0 )
                parity = true;
            fact = -PI / sin(PI * value);
            y = y + 1.0;
            }
        else
            {
            //value = d_huge ( );
            value = HUGE_VAL;
            return value;
            }
        }
    if ( y < Probability::Helper::epsilon() )
        {
        /* argument < EPS */
        if ( xminin <= y )
            {
            value = 1.0 / y;
            }
        else
            {
            //value = d_huge ( );
            value = HUGE_VAL;
            return value;
            }
        }
    else if ( y < 12.0 )
        {
        y1 = y;
        /* 0.0 < argument < 1.0 */
        if ( y < 1.0 )
            {
            z = y;
            y = y + 1.0;
            }
        /* 1.0 < argument < 12.0, reduce argument if necessary */
        else
            {
            n = int ( y ) - 1;
            y = y - ( double ) ( n );
            z = y - 1.0;
            }
        /* evaluate approximation for 1.0 < argument < 2.0 */
        xnum = 0.0;
        xden = 1.0;
        for ( i = 0; i < 8; i++ )
            {
            xnum = ( xnum + p[i] ) * z;
            xden = xden * z + q[i];
            }

        value = xnum / xden + 1.0;
        /* adjust result for case  0.0 < argument < 1.0 */
        if ( y1 < y )
            {
            value = value / y1;
            }
        /* adjust result for case  2.0 < argument < 12.0 */
        else if ( y < y1 )
            {
            for ( i = 1; i <= n; i++ )
                {
                value = value * y;
                y = y + 1.0;
                }
            }
        }
    else
        {
        /* evaluate for 12 <= argument */
        if ( y <= xbig )
            {
            ysq = y * y;
            sum2 = c[6];
            for ( i = 0; i < 6; i++ )
                {
                sum2 = sum2 / ysq + c[i];
                }
            sum2 = sum2 / y - y + sqrtpi;
            sum2 = sum2 + ( y - 0.5 ) * log ( y );
            value = exp ( sum2 );
            }
        else
            {
            //value = d_huge ( );
            value = HUGE_VAL;
            return value;
            }

        }
    /* final adjustments and return */
    if ( parity )
        {
        value = -value;
        }
    if ( fact != 1.0 )
        {
        value = fact / value;
        }

    return value;
}

bool Probability::Helper::isFinite(double x) {
  
    return x > -std::numeric_limits<double>::infinity() && x < std::numeric_limits<double>::infinity();
}

double Probability::Helper::lnBeta(double a, double b)
{
    double corr, p, q;
    
    p = q = a;
    if (b < p) p = b;/* := min(a,b) */
    if (b > q) q = b;/* := max(a,b) */
    
    /* both arguments must be >= 0 */
    if (p < 0)
    {

        std::cout << "Cannot compute log-beta function for a = " << a << " and b = " << b;
    }
    else if (p == 0)
    {
        return -1;
        //return RbConstants::Double::inf;
    }
    else if (!isFinite(q))
    { /* q == +Inf */
        return -1;
        //return RbConstants::Double::neginf;
    }
    
    if (p >= 10) {
        /* p and q are big. */
        corr = lnGammacor(p) + lnGammacor(q) - lnGammacor(p + q);
        return log(q) * -0.5 + 0.918938533204672741780329736406 + corr
            + (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    }
    else if (q >= 10) {
        /* p is small, but q is big. */
        corr = lnGammacor(q) - lnGammacor(p + q);
        return lnGamma(p) + corr + p - p * log(p + q) + (q - 0.5) * log1p(-p / (p + q));
    }
    else
        /* p and q are small: p <= q < 10. */
        return log(gamma(p) * (gamma(q) / gamma(p + q)));
    
}

/*
 * This function calculates the log of the gamma function, which is equal to:
 * Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
 * The result is accurate to 10 decimal places. Stirling's formula is used
 * for the central polynomial part of the procedure.
 *
 * Pike, M. C. and I. D. Hill. 1966. Algorithm 291: Logarithm of the gamma
 *      function. Communications of the Association for Computing Machinery, 9:684.
 */
double Probability::Helper::lnGamma(double a) {

    double x = a;
    double f = 0.0;
    double z;
    if (x < 7)
        {
        f = 1.0;
        z = x - 1.0;
        while (++z < 7.0)
            f *= z;
        x = z;
        f = -log(f);
        }
    z = 1.0 / (x*x);
    return  (f + (x-0.5)*log(x) - x + 0.918938533204673 +
            (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
            0.083333333333333)/x);
}

double Probability::Helper::lnGammacor(double x) {

    const thread_local static double algmcs[15] = {
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
        -.3399615005417721944303330599666e-15,
        +.2683181998482698748957538846666e-17,
        -.2868042435334643284144622399999e-19,
        +.3962837061046434803679306666666e-21,
        -.6831888753985766870111999999999e-23,
        +.1429227355942498147573333333333e-24,
        -.3547598158101070547199999999999e-26,
        +.1025680058010470912000000000000e-27,
        -.3401102254316748799999999999999e-29,
        +.1276642195630062933333333333333e-30
    };
    
    double tmp;
    
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     *   xbig = 2 ^ 26.5
     *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#undef  xmax
#define xmax  3.745194030963158e306
    
    if (x < 10)
    {
        Msg::error("Cannot compute log-gammacor function");
    }
    else if (x >= xmax) {
        Msg::error("Cannot compute log-gammacor function");
        /* allow to underflow below */
    }
    else if (x < xbig) {
        tmp = 10 / x;
        return chebyshevEval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}

double Probability::Helper::incompleteBeta(double a, double b, double x) {
    
    double tol = 1.0E-07;
    
    double value;
    if ( x <= 0.0 )
    {
        value = 0.0;
        return value;
    }
    else if ( 1.0 <= x )
    {
        value = 1.0;
        return value;
    }
    
    /* change tail if necessary and determine S */
    double psq = a + b;
    
    double xx, cx, pp, qq;
    bool indx;
    if ( a < (a + b) * x )
    {
        xx = 1.0 - x;
        cx = x;
        pp = b;
        qq = a;
        indx = true;
    }
    else
    {
        xx = x;
        cx = 1.0 - x;
        pp = a;
        qq = b;
        indx = false;
    }
    
    double term = 1.0;
    int i = 1;
    value = 1.0;
    int ns = (int)(qq + cx * (a + b));
    
    /* use Soper's reduction formulas */
    double rx = xx / cx;
    
    double temp = qq - (double)i;
    if ( ns == 0 )
        rx = xx;
    
    int it = 0;
    int it_max = 1000;
    for (;;)
    {
        it++;
        if ( it_max < it )
        {
            //std::cerr << "Error in incompleteBeta: Maximum number of iterations exceeded!" << std::endl;
            return -1;
        }
        term = term * temp * rx / ( pp + ( double ) ( i ) );
        value = value + term;
        temp = fabs(term);
        if ( temp <= tol && temp <= tol * value )
            break;
        i++;
        ns--;
        if ( 0 <= ns )
        {
            temp = qq - (double)i;
            if ( ns == 0 )
                rx = xx;
        }
        else
        {
            temp = psq;
            psq = psq + 1.0;
        }
    }
    
    /* finish calculation */
    value = value * exp(pp * log(xx) + (qq - 1.0) * log(cx) - lnBeta(a, b) - log(pp));
    if ( indx )
        value = 1.0 - value;
    return value;
}

/*
 * This function returns the incomplete gamma ratio I(x,alpha) where x is
 * the upper limit of the integration and alpha is the shape parameter.
 *
 * Bhattacharjee, G. P. 1970. The incomplete gamma integral. Applied
 *      Statistics, 19:285-287.
 */
double Probability::Helper::incompleteGamma (double x, double alpha, double LnGamma_alpha) {

    double            p = alpha, g = LnGamma_alpha,
                    accurate = 1e-8, overflow = 1e30,
                    rn = 0.0, a = 0.0, b = 0.0, an = 0.0,
                    gin, dif = 0.0, term = 0.0, pn[6];

    if (x == 0.0)
        return (0.0);
    if (x < 0 || p <= 0)
        return (-1.0);

    double factor = exp(p*log(x)-x-g);
    if (x > 1 && x >= p)
        goto l30;
    gin = 1.0;
    term = 1.0;
    rn = p;
    l20:
        rn++;
        term *= x/rn;
        gin += term;
        if (term > accurate)
            goto l20;
        gin *= factor/p;
        goto l50;
    l30:
        a = 1.0-p;
        b = a+x+1.0;
        term = 0.0;
        pn[0] = 1.0;
        pn[1] = x;
        pn[2] = x+1;
        pn[3] = x*b;
        gin = pn[2]/pn[3];
    l32:
        a++;
        b += 2.0;
        term++;
        an = a*term;
        for (int i=0; i<2; i++)
            pn[i+4] = b*pn[i+2]-an*pn[i];
        if (pn[5] == 0)
            goto l35;
        rn = pn[4]/pn[5];
        dif = fabs(gin-rn);
        if (dif>accurate)
            goto l34;
        if (dif<=accurate*rn)
            goto l42;
    l34:
        gin = rn;
    l35:
        for (int i=0; i<4; i++)
            pn[i] = pn[i+2];
        if (fabs(pn[4]) < overflow)
            goto l32;
        for (int i=0; i<4; i++)
            pn[i] /= overflow;
        goto l32;
    l42:
        gin = 1.0-factor*gin;
    l50:
        return (gin);
}

void Probability::Helper::normalize(std::vector<double>& vec) {

    double sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        {
        sum += vec[i];
        if (vec[i] < 0.0)
            Msg::error("Cannot normalize a vector with negative elements");
        }
    for (int i=0; i<vec.size(); i++)
        vec[i] /= sum;
}

double Probability::Helper::pointNormal(double prob) {

	double a0 = -0.322232431088;
	double a1 = -1.0;
	double a2 = -0.342242088547;
	double a3 = -0.0204231210245;
 	double a4 = -0.453642210148e-4;
 	double b0 = 0.0993484626060;
 	double b1 = 0.588581570495;
 	double b2 = 0.531103462366;
 	double b3 = 0.103537752850;
 	double b4 = 0.0038560700634;
 	double p = prob;
	double p1 = ( p < 0.5 ? p : 1.0 - p);
	if (p1 < 1e-20)
	   return (-9999.0);
	double y = sqrt( log(1.0/(p1*p1)) );
	double z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return ( p < 0.5 ? -z : z );
}

double Probability::Helper::rndGamma(RandomVariable* rng, double s) {

    double r = 0.0;
    if (s <= 0.0)
        std::cout << "Gamma parameter less than zero" << std::endl;
    else if (s < 1.0)
        r = Probability::Helper::rndGamma1(rng, s);
    else if (s > 1.0)
        r = Probability::Helper::rndGamma2(rng, s);
    else
        r = -log(rng->uniformRv());
    return (r);
}

double Probability::Helper::rndGamma1(RandomVariable* rng, double s) {

    double            r, x = 0.0, small = 1e-37, w;
    thread_local static double   a, p, uf, ss = 10.0, d;
    
    if (s != ss)
        {
        a  = 1.0 - s;
        p  = a / (a + s * exp(-a));
        uf = p * pow(small / a, s);
        d  = a * log(a);
        ss = s;
        }
    for (;;)
        {
        r = rng->uniformRv();
        if (r > p)
            {
            x = a - log((1.0 - r) / (1.0 - p));
            w = a * log(x) - d;
            }
        else if (r>uf)
            {
            x = a * pow(r / p, 1.0 / s);
            w = x;
            }
        else
            return (0.0);
        r = rng->uniformRv();
        if (1.0 - r <= w && r > 0.0)
        if (r*(w + 1.0) >= 1.0 || -log(r) <= w)
            continue;
        break;
        }
    
    return (x);
}

double Probability::Helper::rndGamma2(RandomVariable* rng, double s) {

    double            r, d, f, g, x;
    thread_local static double    b, h, ss = 0.0;
    
    if (s != ss)
        {
        b  = s - 1.0;
        h  = sqrt(3.0 * s - 0.75);
        ss = s;
        }
    for (;;)
        {
        r = rng->uniformRv();
        g = r - r * r;
        f = (r - 0.5) * h / sqrt(g);
        x = b + f;
        if (x <= 0.0)
            continue;
        r = rng->uniformRv();
        d = 64.0 * r * r * g * g * g;
        if (d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))
            break;
        }
    
    return (x);
}
