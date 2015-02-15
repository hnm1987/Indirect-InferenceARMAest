#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#import <maximize>
#import <lib/hacest>
#include "funcprepar.ox"
#include "armasim.ox" //version 5 is a lot faster
#include "estar.ox"  //version 2 is faster by a second or two for each loop, makes a lot of difference for large number of loops
//#import <packages/maxsa/maxsa>
#import <maxsa/maxsa>

//global constants
static decl g_S =100;
static decl seed = 1;

decl g_n=250;

decl g_para, g_ar, g_ma;

decl g_r, g_i, g_beta, g_sigma;

decl g_bi;

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
	  
	mStore=zeros(g_S,ilags+2);
	
	//generate g_S series of ARMA(p, q)
	vX = funcsimarma(varti_para, g_n, g_ar, g_ma, g_S);
   
	//estimate an AR(ilags=g_ar+g_ma+g_r) for each series of ARMA(p, q), store the estimated parameters
	for (i=0;i<g_S;i++)
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

	funcmeanbeta(&dfval,vp_temp);

	adfunc[0]=-(g_beta-dfval)'*invertgen(g_sigma)*(g_beta-dfval);

	//safeguarding
	if(ismissing(adfunc[0]))
	{
	 	adfunc[0] = -1e5;
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
	decl ir, iQ;
			
	vp_temp=g_para; //the initial theta
	mh = unit(sizer(vp_temp));	//the initial hessian

	if(g_bi==1)
	{
	decl dT = 3;
	
	MaxSAControl(1e5, 2);
	MaxSAControlStep(10, 2, 0.5, 0.1, 2);
	ir = MaxSA(funcmin, &vp_temp, &dfval, &dT);
	
	println("\n",MaxSAConvergenceMsg(ir), " at parameters ", vp_temp[:], "\nwith function value ", double(dfval));
	}
	else
	{
	MaxControl(1000, 1, 0);
	ir = MaxBFGS(funcmin, &vp_temp, &dfval, &mh, TRUE);
	println("\n",MaxConvergenceMsg(ir), " at parameters ", vp_temp[:], "\nwith function value ", double(dfval));
	}


	
	avp_opt[0]=vp_temp;
	
	NumJacobian(funcmeanbeta, vp_temp, &mjacob);	//compute the Jacobian at the estimated theta
	amJacobian[0]=mjacob;
	mW=(1+1/g_S)*invertgen(mjacob'*invertgen(g_sigma)*mjacob);
	amW[0]=mW;

	return 1;
}

funcestarma(const vY, const avp, const amW, const ip, const iq, const ir, const bi)
{
/*
purpose- 
to return the optimized theta and W matrix in addresses

input-
vY: the column vector of interest
avp: an adress to store the estimated theta
amW: an address to store the estimated W matrix
ip: number of AR terms
iq: number of MA terms
ir: number of overidentifying terms
bi: if 1 use simulated annealing, else use BFGS

output-
return 1, vector theta and matrix W stored in their respective addresses
*/

decl ilags, vTheta, mW, mjacob;

g_ar = ip;
g_ma = iq;
g_r = ir;

g_n = sizer(vY);
ilags=g_ar+g_ma+g_r;
g_para = avp[0];

g_bi = bi; //if 1 use simulated annealing, else use BFGS

	funcestar(vY,ilags,&g_beta, &g_sigma); //estimate beta hat and sigma hat to be used by callmin later
	callmin(&vTheta,&mW,&mjacob);

	//assign theta and W to their addresses
	avp[0] = vTheta;
	amW[0] = mW;

	return 1;

}