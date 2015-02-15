#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#import <maximize>
#import <lib/hacest> //important

funcestar(const vX, const ilags, const avp, const amSigma)
{

/*

purpose- 
to estimate an ar model

input-
vX: the column vector of interest to estimate
ilags: number of lagged terms to estimate
avp: address of vector of parameters
amSigma: address of the heteroskedastic consistent sigma matrix

output-
return 1 and the vector of parameters and sigma matrix in their respective addresses

*/

decl mXX,mX,mxx_temp,mhac,mHACSE;
decl vbeta,vY,vres,vZt;
decl iN,iK;
decl ds_ols,ds_wht;

//initialization
mxx_temp= funcprepar(vX,ilags);
vY=mxx_temp[][0];
mX=mxx_temp[][1:];
iN=sizer(mxx_temp);
iK=sizec(mxx_temp);
mXX  = (mX'mX);

	//estimation of betas
	vbeta    = invertgen(mXX)*mX'vY;
	
	//estimate heteroskedastic consistent sigma
	vres   = vY - mX*vbeta;
	ds_ols = vres'vres/(iN-iK);  //adjusts degrees of freedom
	vZt    = mX.*vres;
	ds_wht = vZt'vZt;
					
	mhac=HACest_tr(vZt,vZt,-1)*iN;
	mHACSE = zeros(iK, iK);
	mHACSE[0:iK-2][0:iK-2]  = iN/(iN-ilags-iK)*invertsym(mXX)*mhac*invertsym(mXX)';

	//final outputs
	mHACSE[iK-1][iK-1]=2*ds_ols^4/(iN-iK);	//variance of sigma2 = 2*sigma4/(n-k)
	avp[0]=vbeta|ds_ols;
	amSigma[0]=mHACSE;

	return 1;
	
}

