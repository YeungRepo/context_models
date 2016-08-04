function[mOut] = m_gyrtopfunc(sigmaX,pDNA,s0,T0,G0,tau,gamma,TLX,pLen,kMgyr)
s0 = s0*pLen/TLX;

s0m = abs(s0);
pX = (sigmaX>s0);%
nX = (sigmaX<s0);
Hill = abs(sigmaX-s0)/kMgyr/(s0m+((sigmaX-s0)/kMgyr)^2);
mOut = nX*Hill*(T0*tau/pDNA)-pX*Hill*(G0*gamma/pDNA);

