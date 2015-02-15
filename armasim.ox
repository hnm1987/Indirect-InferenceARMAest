#include <oxstd.h>

funcsimarma(const vp, const iN, const iAR, const iMA, const iS)
{

/*
purpose-
to generate a matrix of arma(p, q) random variables

input-
vp: the column vector of ARMA(p, q) parameters in the form iAR(1)|iAR(2)|...|iAR(p)|iMA(1)|iMA(2)|...|iMA(q)|dMu|dSigma2 
iN: number of observations per ARMA vector
iAR: number of AR parameters
iMA: number of MA parameters
iS: number of  ARMA series
		
output-	
an ARMA matrix of size iN x iS
*/

decl vEps, mY; //vector of errors and matrix of final output
decl i; //loop variable 
decl dMu, dSigma2, vAR_para, vMA_para;	//initializing arma parameters
decl m, p, q, iT, unc_mean, iDiscard; //other initializing constants
	
	//initialize constants
	dMu = vp[sizer(vp)-2];
	dSigma2 = vp[sizer(vp)-1];
	m = max(iAR, iMA);
	iDiscard = 2000;	//will be used to avoid biases due to initial conditions	
	iT = iN+iDiscard;
	
	//input checking for AR part
	if(iAR>=1)
	{
	vAR_para = vp[0:iAR-1];
	p = iAR;
	}
	else
	{
	vAR_para = 0;
	p = 1;
	}

	//input checking for MA part
	if(iMA>=1)
	{
	vMA_para = 1|vp[iAR:iAR+iMA-1];
	q = iMA;
	}
	else
	{
	vMA_para = 1;
	q = iMA;
	}

	unc_mean = dMu/(1-sumc(vAR_para));	//unconditional mean of ARMA, used to initialize the first few lags
	vEps = sqrt(dSigma2)*rann(iN+iDiscard, iS);
	mY = zeros(iN+iDiscard, iS);	
	mY[0:m-1][:]= unc_mean*ones(m, 1);	//initialize the first few lags
	 
		
	for(i=m;i<iT;i++)
	{
	  	mY[i][:] = vAR_para'*reversec(mY[i-p:i-1][:]) + vMA_para'*reversec(vEps[i-q:i][:]);					   
	}

	mY = dMu+mY;	//add mu		

	return mY[iDiscard:iT-1][];	//to avoid some biases due to initial conditions, start taking observations after the discard-th observation			
	
}


