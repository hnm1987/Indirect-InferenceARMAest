#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#import <maximize>
#import <lib/hacest>
#include "funcprepar.ox"
#include "armasim.ox" //version 5 is a lot faster
#include "estar.ox"  //version 2 is faster by a second or two for each loop, makes a lot of difference for large number of loops

//global constants
static decl n=250;
static decl iS =100;
static decl seed=1;

decl g_para, g_ar, g_ma;

decl g_r, g_i, g_beta, g_sigma;

funcmeanbeta(const adfval, const vp)
{
/*
purpose- 
to approximate the binding function b(theta) by computing the sample mean of beta(s)

input-
adfval: address to store the function value evaluated at vp_temp
vp_temp: column vector of inputs in the form AR paras|MA paras|mu|sigma2

output-
return 1, function value stored in the address of adfval
*/
    ranseed(seed+g_i^2);
    decl vX,i,vp_temp,mSigma_temp,ilags,mStore,varti_para,irowvp;
   
    irowvp=sizer(vp);
    ilags=g_ar+g_ma+g_r;
   
    //take absolute value of sigma2 to prevent negative value problem
	varti_para=vp;
    varti_para[irowvp-1]=fabs(varti_para[irowvp-1]);
	  
	mStore=zeros(iS,ilags+2);
	
	//generate iS series of ARMA(p, q)
	vX = funcsimarma(varti_para, n, g_ar, g_ma, iS);
   
	//estimate an AR(ilags=g_ar+g_ma+g_r) for each series of ARMA(p, q), store the estimated parameters
	for (i=0;i<iS;i++)
	{

		funcestar(vX[][i],ilags,&vp_temp,&mSigma_temp);
		mStore[i][]=vp_temp';
   
   }
   
   //calculate the sample mean of the estimated parameters and return through an address
   adfval[0]=meanc(mStore)';

   return 1;
   
}	 

funcmin(const vp_temp, const adfunc, const avscore, const amhessian)
{
/*
purpose- 
the quadratic function, Q(theta) to minimize

input-
vp_temp: column vector of inputs in the form AR paras|MA paras|mu|sigma2
adfunc: address to store the function value evaluated at vp_temp
avscore: address to store the score vector
amhessian: address to store the hessian

output-
return 1, function value, score and hessian stored in their respective addresses
*/

	decl dfval;

	funcmeanbeta(&dfval,vp_temp)  ;

	adfunc[0]=-(g_beta-dfval)'*invertgen(g_sigma)*(g_beta-dfval);

	//safeguarding
	if(ismissing(adfunc[0]))
	{
	 	adfunc[0] = -10;
	}

	return 1;
	
}

callmin(const avp_opt, const amW, const amJacobian)
{
/*
purpose- 
to return the optimized theta, W matrix and Jacobian in addresses

input-
vp_temp: an address to store the estimated theta
amW: an address to store the estimated W matrix
amJacobian: an address to store the estimated Jacobian matrix

output-
return 1, matrix W and Jacobian stored in their respective addresses
*/

	oxwarning(0); // all warnings off

	decl dfval,mjacob,mW;
	decl mh, vp_temp;
			
	vp_temp=g_para; //the initial theta
	mh = unit(sizer(vp_temp));	//the initial hessian
	
	MaxControl(50, 0, 0);
	MaxBFGS(funcmin, &vp_temp, &dfval, &mh, TRUE);
	
	avp_opt[0]=vp_temp;
	
	NumJacobian(funcmeanbeta, vp_temp, &mjacob);	//compute the Jacobian at the estimated theta
	amJacobian[0]=mjacob;
	mW=(1+1/iS)*invertgen(mjacob'*invertgen(g_sigma)*mjacob);
	amW[0]=mW;

	return 1;
}

saveMC(savetheta, saveW, iR)
{

/*
purpose- 
to save the Monte Carlo experiment's parameters vector and variance matrix

input-
savetheta: filename to save the parameters vector
saveW: filename to save the diagonal of the variance matrix
iR: number of replications

output-
saved files of parameters vector and variance matrix

*/

decl vTheta, mW, mjacob, mStoreTheta, mStoreW, vY, ilags;
decl time0, timeloop1, timeloop2;

ilags=g_ar+g_ma+g_r;
time0=timer();	//start timer for the whole Monte Carlo
mStoreTheta=zeros(g_ar+g_ma+2,iR);
mStoreW=zeros(g_ar+g_ma+2,iR);

	 for (g_i=1;g_i<=iR;g_i++)
	 {
	 
		vY=funcsimarma(g_para, n,g_ar,g_ma, 1) ;
		funcestar(vY,ilags,&g_beta, &g_sigma); //estimate beta hat and sigma hat to be used by callmin later
		timeloop1=timer();	//start timer for each Monte Carlo replication
		callmin(&vTheta,&mW,&mjacob);
		 
		mStoreTheta[][g_i-1]=vTheta;
		mStoreW[][g_i-1]=diagonal(mW)';
		 
		//printings
		println("i= ",g_i);
		println("g_r= ", g_r);
		timeloop2=timer();	//stop timer for each Monte Carlo replication
		println("time elapsed each loop=",timespan(timeloop1,timeloop2));
	
	 }
	
	//savings
	savemat(savetheta,mStoreTheta);
	savemat(saveW,mStoreW);
	
	//print the total time needed for the whole Monte Carlo
	println("total time elapsed for MC=",timespan(time0,timer())); 

}



