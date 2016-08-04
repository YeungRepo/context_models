
clc
digits 32

maint_on = 1;
%Physiological scaling factors --- concentration conversion, time unit
%conversion, and volume conversion (dilation factors) 
Omega = 1; %dilfac and Omega are functionally equivalent. 
%1.4e1;%1e9; %10.5 works.%uLiters %(conversion to nL)
csf = 1; % 1e-3 means everything is in µM, 1 is nM, all orig. units in nM. 
tdf = 1; % time-dilation factor,1 is in seconds, 60 is in minutes, 60*60 is in hours. all original units in /s 
dilfac = 1/30; %Noireaux paper ranges 1/20 to 1/30 
cgain = 1;%1.4e1; %fold-change factor, measure of parametric uncertainty; cgain = 1+x defines a x*100% parametric uncertainty

%Length parameters (of DNA) and turn parameters. 
h0 = 10.5; %Stryer, Wang Liu 1984 

% - - - - - - - - - - - - -
% - - - - - - - - - - - - - 
% TLS and TLG are the transcript length parameters

TLS = 818; %CFP gene length 
TLG = 758; %RFP gene length 

% - - - - - - - - - - - - 
% - - - - - - - - - - - - 



NS = 150; %as per our synthesis
PLS = 40; %pLac promoter length BBR reference 
PLG = 44; %pTet promoter length BBR reference 
pLen = 2500; %plasmid length, used to calculate global supercoiling state. 
std_coil = 50; %125 is a good estimate for mSpinach MG/ 

%concentrations of enzymes, DNA  etc. 
Rtot = dilfac*18931*csf; % or 2490-18931 nM Bremer Dennis 1996 polymerases in log phase e coli
Ribotot = dilfac*11291*csf; % nM ribosomes in log phase e coli (Bremer Dennis 1996)
pLactot = 11*csf; %nM from XPT107.2 TXTL template using the OLD-TRI stock (seq verified with some intact pLac)
pTettot = pLactot; %identical since on the same plasmid.
k_catm = 1e0;
%% 
G0 = k_catm*maint_on*dilfac*12*csf; %nM Maier, Schmidt  .. Serrano est. Mcoplasmisa pneumoniaefrom Mol Syst Biol. 2011 
T0 = 2; %nM handfitted 

LacItot = 10*dilfac*csf ; % nM  in vitro TXTL cond.est. from Kalisky, Dekel, cAlon Phys. Biol. 2007  
TetRtot = 10*csf; % nM  in vitro TXTL cond. 
IPTGtot = 0*csf;
aTctot = 100*csf;
Dtot= 2441*dilfac; %nM Mackie, Nature. Review. MicroBiology 2013 

%transcription/prod. rates,   550 = kr/kfmax = k 
kfmaxl = cgain*tdf*7*1e-2; %nM/s taken from Tusa, Gaskins, Kim, and Sdeveyanski  % 
kr = cgain*tdf*550*kfmaxl;%   per s   Calculated from Bintu, Buchler,  ... Phillips,Curr. Op. Genetics 2005 
kcatmax = cgain*tdf*85/((TLS+TLG)/2);%>70 is okay %40-90 nt per second see  Condon c, French s, Squires c. etclibrary.bib Bremer1996
kopen = cgain*tdf*.04; % Buc H, McClure, WR Biochemistry 1985 
kl = cgain*tdf*.02; % Larson2008 applied 
kw = cgain*tdf*1;
rhol = cgain*csf*dilfac*tdf*0; %in vitro case, no production of LacI 
rhot = cgain*csf*dilfac*tdf*0;  %in vitro case, no production of TetR  
kTL = cgain*tdf*21/(mean([TLS TLG]))*3.0e-1;
kfoldcfp = 1/(20*60);
kfoldrfp = 1/(70*60); %should be 100 min maturation time, 
kfolddimer = 1/(70*60);
BCD16 = .007;
BCD_Bgl = .018;
kfmaxt = kfmaxl*.22; %estimate 22 percent strength of Lac promoter at full on. 
%degradation rates. 
dms =log(2)/(30*60);%;assumes an exponential decay rate , 10 min from Taniguchi !!including cell division, 
dmg = log(2)/(60*60);
% the corrected rate satisfies : log(2)/log(10*60) = log(2)/log(dm*60) + log(2)/log(30*60)  
%dm  = log(2)/60*1/(log(2)/(10*60) - log(2)/(30*60)); 
kMd = 10;

%gyrase/topoisomerase rates. 

tau = cgain*tdf*.5; %Wang and Liu 1984
gamma =cgain*tdf*.5; %Wang and LIu 1984
s0 = -0.065; %est. Rhee & Hatfield 
kMgyr = 200; 
%LacI binding, induction, sequestration rates 
kaL = cgain/csf*tdf*6*1e3; %/nM /s  %Uses parameter for allolactose substrate, ...
%% assuming IPTG and allolactose have similar affinity/association constant. Kalisky, Dekel, Alon Phys. Biol. 2007  
kuaL = cgain*tdf*1; % /s Xie, Choi, Li, ... Lia, 2008  Annua. Rev. Biophy 
kseqL = cgain*tdf*1e1; % nM/s  Calculated using assoc. constant from Kalisky, Dekel, Alon Phys. Biol. 2007  
kuL = cgain*tdf*.022; % Nelson Sauer, 1985 

%TetR binding, induction, sequestration rates 
kaT = kaL; %no lit. values
kuaT = kuaL; %no lit. values 
kseqT  = kseqL; % no lit. values 
kuT = kuL; % no lit values. 
dp = cgain*dilfac*tdf*30; % in vitro case, no degradation of LacI 



options=odeset('AbsTol',1e-5); 
%reset time back to its original state. 
protein_hours =7 ;
t =  [0 protein_hours*60*60/tdf];%
gyrpet= -4.5 ;%-.225;%-.225*s0;
mS0 = 0;% nM
mG0 = 0; % nM
ECS0=0; % nM
ECG0 = 0; %nM
CCS0 = 0; % nM
CCG0 = 0; % nM;
stS0 =-4; %unitless, supercoiling density is unitless
stG0 =5;
spS0 =-4;
spG0 =5;
ECGECS0= 0;
LacI0  = LacItot; % nM; in vitro 
TetR0 = TetRtot; % nM;  in vitro 
IPTG0 = 0; % 0 nM, cons. with exp. conditions. 
aTc0 = 0; % 0 nM, cons. with exp. conditions.
PmS0 = 0;
PMG0 = 0;


stC0 =  -4.5;
stR0 = -4.5;
spC0 =  -4.5;
spR0 = -1; 
cd src;

x0 = define_newx0(stC0,stR0,spC0,spR0,LacItot,TetRtot);
[t,xcp]=ode23s('ConvModelRFPCFP_Journal',t,x0,options,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
   PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer);
clear x0

stC0 = -30;
stR0 = 30;
spC0 = 20;
spR0 = 12; 
x0 = define_newx0(stC0,stR0,spC0,spR0,LacItot,TetRtot);


[t,xdp] = ode23s('DivModelRFPCFP_Journal',t,x0,options,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
    PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer);

clear x0

stC0 = 25;
stR0 = -4;
spC0 = -1;
spR0 = 3; 
x0 = define_newx0(stC0,stR0,spC0,spR0,LacItot,TetRtot);

[t,xtp] = ode23s('TandModelRFPCFP_Journal',t,x0,options,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
    PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer);

 figure(3);
%axes3 = axes('Parent',fig3,'FontSize',fs)
%Convergent expression dynamics 
hold on


t = t*tdf/(60*60);
lws = 3;
ms = 3;
cmat1 = {[0 0 .9]; [.7 .7 .7]; [.45 .65 .65]; [0 0 0]};
cmat2 = {[.9 0 0]; [.7 .7 .7]; [155/255 90/255 90/255]};
plot(t,xcp(:,17)/csf,'s-','Color',cmat1{1},'LineWidth',lws,'MarkerSize',ms);
plot(t,xdp(:,17)/csf,'s-','Color',cmat1{2},'LineWidth',lws,'MarkerSize',ms);
plot(t,xtp(:,17)/csf,'s-','Color',cmat1{3},'LineWidth',lws,'MarkerSize',ms);
%title('mSpinach gene expression versus time','FontSize',24)
h_leg = legend('Con','Div','Tan');
set(h_leg,'Location','Best')
set(h_leg,'EdgeColor','White','FontSize',24)
%xlabel('Time','FontSize',24)
%ylabel('Expression (nM)','FontSize',24)
set(gca,'FontSize',24)
hold off

 figure(4);
hold on

plot(t,xcp(:,19)/csf,'s-','Color',cmat2{1},'LineWidth',lws,'MarkerSize',ms);
plot(t,xdp(:,19)/csf,'s-','Color',cmat2{2},'LineWidth',lws,'MarkerSize',ms);
plot(t,xtp(:,19)/csf,'s-','Color',cmat2{3},'LineWidth',lws,'MarkerSize',ms);
h_leg = legend('Con','Div','Tan');
set(h_leg,'Location','Best')
set(h_leg,'EdgeColor','White','FontSize',24)
%xlabel('Time','FontSize',24)
%ylabel('Expression (nM)','FontSize',24)
set(gca,'FontSize',24)
hold off




fig5 = figure(5);
%axes3 = axes('Parent',fig3,'FontSize',fs)
%Convergent Supercoiling
hold all
t = t*tdf/(60*60);

lis = 5;
ms = 50;
fs= 36;


plot(t,xcp(:,5),'LineWidth',lis)
plot(t,xcp(:,6),'LineWidth',lis)
plot(t,xcp(:,7),'LineWidth',lis)
plot(t,xcp(:,8),'LineWidth',lis)
title('Convergent supercoiling versus time','FontSize',24)
% xlabel('Time','FontSize',24)
% ylabel('Supercoiling Density','FontSize',24)
 legend('stS','stG','spS','spG');
% set(gca,'FontSize',24)
% hold off



%%Divergent Supercoiling
fig6 = figure(6);
%axes4 = axes('Parent',fig4,'FontSize',fs)
hold all 
plot(t,xdp(:,5)/Omega,'LineWidth',ms/4)
plot(t,xdp(:,6)/Omega,'LineWidth',ms/4)
plot(t,xdp(:,7)/Omega,'LineWidth',ms/4)
plot(t,xdp(:,8)/Omega,'LineWidth',ms/4)
hold off
title('Divergent supercoiling versus time','FontSize',24)
xlabel('Time','FontSize',14)
ylabel('Supercoiling Density','FontSize',24)
legend('stS','stG','spS','spG');
set(gca,'fontsize',24)
% 
% %Tandem Supercoiling
fig7 = figure(7);
%axes5 = axes('Parent',fig5,'FontSize',fs)
hold all
plot(t,xtp(:,5)/Omega,'LineWidth',ms/4)
plot(t,xtp(:,6)/Omega,'LineWidth',ms/4)
plot(t,xtp(:,7)/Omega,'LineWidth',ms/4)
plot(t,xtp(:,8)/Omega,'LineWidth',ms/4)
hold off
title('Tandem supercoiling versus time','FontSize',24)
xlabel('Time','FontSize',24)
ylabel('Supercoiling Density','FontSize',24)
legend('stS','stG','spS','spG')
set(gca,'fontsize',24)

cd .. 
