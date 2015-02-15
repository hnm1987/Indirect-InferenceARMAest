#include <oxstd.h>

funcprepar(const vx,const lags) 
{

/*
purpose-
to append a vector with its own lagged values and a vector of ones for regression

input-
vx: the column vector of interest
ilags: number of lags to append
		
output-	
the original column vector appended with lagged vectors	and a vector of ones
*/

decl mxx,i,iT;
iT=rows(vx);
mxx=zeros(iT-lags,1+lags);

	for (i=0;i<=lags;i++)

	{
	mxx[][i]=vx[lags-i:iT-i-1];
	 
	}

return mxx~1;	 

}


