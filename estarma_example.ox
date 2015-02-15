#include <oxstd.h>
#include "estarma.ox"

/*
to demonstarate the use of our algorithm to estimate an arma model with a small example
*/

main()
{

decl vY = funcsimarma(<0.5;0.3;1;2>,250,1,1,1); //simulate a 250 x 1 ARMA(1, 1) vector with phi1 = 0.5, alpha1 = 0.3, mu = 1 and sigma2 = 2
decl vTheta, mW;
decl p = 1, q = 1, r = 0;	//but I estimate an ARMA(1, 1) with r = 0

vTheta = <0.6;0.2;1.2;2.25>;	//initial value to start the algorithm
decl bi;
//bi = 1; //if 1 use simulated annealing, else use BFGS
bi = 1;

funcestarma(vY, &vTheta, &mW, p, q, r, bi);	 //the function that does the estimation

println("vTheta = ", vTheta);
println("mW = ", mW);
	

}
