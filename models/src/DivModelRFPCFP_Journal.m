
function[dx] = DivModelJournal(t,x,x0,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
    PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer)
%State definitions 
mS = x(1);
MG = x(2);
ECS = x(3);
ECG = x(4);
stS = x(5);
stG = x(6);
spS = x(7);
spG = x(8);
CCS = x(9);
CCG = x(10);
LacI  = x(11);
TetR = x(12);
IPTG = x(13);
aTc = x(14); 
PmS = x(15);
PMG = x(16);
fmS = x(17);
fMG = x(18);
ffMG = x(19);


%auxiliary functions/definitions
pLacC = LacItot-LacI-IPTGtot+IPTG;
pTetC = TetRtot-TetR-aTctot+aTc;
pLac = pLactot-ECS-CCS-pLacC;
pTet = pTettot-ECG-CCG-pTetC;
R = Rtot-ECS-ECG-CCS-CCG;

DeltaKink = (spS+spG)*h0;

nfS = PLS+ NS/2 -DeltaKink;
nfG = PLG + NS/2 +DeltaKink;

                if nfS <0 % indicates complete override of MG on mSpinach, no more mSpinach transcription. 
                    nfS = -Inf;
                    nfG = PLS+TLS + TLG + NS;
                end
                if nfG <0 % indicates complete override of mSpinach coils on MG, no more MG transcription. 
                    nfG = -Inf; 
                    nfS = PLG+TLS + TLG+ NS;
                end



kfS = kfmaxfunc(PLS,spS,s0,kfmaxl,std_coil,pLen);
kfG = kfmaxfunc(PLG,spG,s0,kfmaxt,std_coil,pLen);



kcatS = kcatmaxfunc(t,TLS,stS,s0,kcatmax,std_coil,pLen);
kcatG = kcatmaxfunc(t,TLG,stG,s0,kcatmax,std_coil,pLen);



%state update
%dmS 
dx(1) = kcatS*ECS-dms*mS;%(Dtot*mS/(kMd+mS));
%dmG
dx(2) = kcatG*ECG-dmg*MG;%(Dtot*MG/(kMd+MG));
%dECs
dx(3) = kopen*CCS-(kcatS)*ECS;
%dECG
dx(4) = kopen*CCG-(kcatG)*ECG;



%CCS 
dx(9) = kfS*(Rtot-ECS-ECG-CCS-CCG)*(pLactot-CCS-ECS-pLacC) - (kr+kopen)*CCS; 

%CCG
dx(10) = kfG*(Rtot-ECS-ECG-CCS-CCG)*(pTettot-CCG-ECG-pTetC) - (kr+kopen)*CCG;


%LacI
dx(11) = rhol -kaL*LacI*IPTG+ kuaL*(LacItot-LacI-pLacC)+kuL*pLacC-kseqL*pLac*LacI-dp*LacI;
%TetR 
dx(12) = rhot -kaL*TetR*aTc + kuaT*(TetRtot-TetR-pTetC)+kuT*pTetC-kseqT*pTet*TetR-dp*TetR;
%IPTG
dx(13) = -kaL*(LacI+pLacC)*IPTG + kuaL*(LacItot-LacI-pLacC);
%aTc 
dx(14) = -kaT*(TetR + pTetC)*aTc + kuaT*(TetRtot-TetR - pTetC);


%PmS 
dx(15) = BCD16*kTLfunc(t,kTL)/(TLS/3)*Ribotot*mS -kfoldcfp*PmS;
%PMG
dx(16) = BCD_Bgl*kTLfunc(t,kTL)/(TLG/3)*Ribotot*MG-kfoldrfp*PMG;

%fmS
dx(17) = kfoldcfp*PmS;
%fMG
dx(18) = kfoldrfp*PMG-kfolddimer*fMG;
dx(19) = kfolddimer*fMG;

dECS = dx(3);
dECG = dx(4);
dCCS = dx(9);
dCCG = dx(10);
dmScr = dx(1)+dms*mS;
dMGcr = dx(2) + dmg*MG;



%sigmatS
dx(5) = Omega*(1/(2*h0)*((dECS-dCCS)*PLS/TLS) ... 
     + m_gyrtopfunc(stS,pLactot,s0,T0,G0,tau,gamma,TLS,pLen, kMgyr)); 

%sigmatG
dx(6) = Omega*(-1/(2*h0)*((dECG-dCCG)*PLG/TLG)...
    + m_gyrtopfunc(stG,pTettot,s0,T0,G0,tau,gamma,TLG,pLen, kMgyr));

%sigmapS
dx(7) =  Omega*(-1/(2*h0)*((dECS-dCCS)*PLS/nfS)  ...
    +m_gyrtopfunc(spS,pLactot,s0,T0,G0,tau,gamma,PLS,pLen, kMgyr));
%sigmapG
dx(8) =Omega*(1/(2*h0)*((dECG-dCCG)*PLG/nfG )...
    + m_gyrtopfunc(spG,pTettot,s0,T0,G0,tau,gamma,PLG,pLen, kMgyr));

dx = dx';

