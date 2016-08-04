
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
% TLS and TLG are the RNA transcript length parameters

TLS = 203;% mSpinach aptamer, Paige Science 2011
TLG = 68; % MG aptamer, Babendure 2003 

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


hours =3.1;
t =  [0 hours*60*60/tdf];%linspace(0,hours*60*60/tdf,25)';
options=odeset('AbsTol',1e-5); 

% 
% 
% %%mSpinach Simulation ---- 


 %localized supercoiling density is much higher than global supercoiling density, 
% see Supplementary for comments 

stS0 = -6;
stG0 = -3;
spS0 = -6;
spG0 = -3; 

cd src;

x0 = define_newx0(stS0,stG0,spS0,spG0,LacItot,TetRtot);
[t,xc]=ode23s('ConvModelRFPCFP_Journal',t,x0,options,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
   PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer);
clear x0


stS0 = 25;
stG0 = -5;
spS0 = 11;
spG0 = 5; 

x0 = define_newx0(stS0,stG0,spS0,spG0,LacItot,TetRtot);
[t,xd] = ode23s('DivModelRFPCFP_Journal',t,x0,options,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
    PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer);

clear x0

stS0 = -2;
stG0 = 4;
spS0 = 4;
spG0 = 2; 

x0 = define_newx0(stS0,stG0,spS0,spG0,LacItot,TetRtot);
[t,xt] = ode23s('TandModelRFPCFP_Journal',t,x0,options,dms,dmg,Rtot,pLactot,kr,pTettot,Omega,h0,TLS,TLG,NS,PLS,...
    PLG,tau,gamma,kfmaxl,kfmaxt,kcatmax,s0,kl,kw,T0,G0,rhol,kaL,kuaL,kseqL,kuL,rhot,kaT,kuaT,kseqT,kuT,dp,...
    LacItot,TetRtot,IPTGtot,aTctot,kopen,Dtot,kMd,Ribotot,kTL,BCD16,BCD_Bgl,pLen,kMgyr,std_coil,kfoldcfp,kfoldrfp,kfolddimer);



set(0,'DefaultFigureWindowStyle','docked')


fig1 = figure(1)
fig2 = figure(2)

%axes1 = axes('Parent',fig1,'FontSize',24)

AFUscale = 1;
AFUscaleMG = 1; 

lis = 5;
ms = 50;
fs= 36;




fig5 = figure(5);

%Convergent Supercoiling
hold all
t = t*tdf/(60*60);

plot(t,xc(:,5),'LineWidth',lis)
plot(t,xc(:,6),'LineWidth',lis)
plot(t,xc(:,7),'LineWidth',lis)
plot(t,xc(:,8),'LineWidth',lis)
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
plot(t,xd(:,5)/Omega,'LineWidth',ms/4)
plot(t,xd(:,6)/Omega,'LineWidth',ms/4)
plot(t,xd(:,7)/Omega,'LineWidth',ms/4)
plot(t,xd(:,8)/Omega,'LineWidth',ms/4)
hold off
title('Divergent supercoiling versus time','FontSize',24)
xlabel('Time','FontSize',14)
ylabel('Supercoiling Density','FontSize',24)
legend('stS','stG','spS','spG');
set(gca,'fontsize',24)
% 
% %Tandem Supercoiling
fig7 = figure(7)
%axes5 = axes('Parent',fig5,'FontSize',fs)
hold all
plot(t,xt(:,5)/Omega,'LineWidth',ms/4)
plot(t,xt(:,6)/Omega,'LineWidth',ms/4)
plot(t,xt(:,7)/Omega,'LineWidth',ms/4)
plot(t,xt(:,8)/Omega,'LineWidth',ms/4)
hold off
title('Tandem supercoiling versus time','FontSize',24)
xlabel('Time','FontSize',24)
ylabel('Supercoiling Density','FontSize',24)
legend('stS','stG','spS','spG')
set(gca,'fontsize',24)

cmat1 = {[0 .9 0 ],[.7 .7 .7],[0 .4 .1]};
figure(fig1)
hold on

lws = 4;
plot(t,xc(:,1)/csf,'s-','Color',cmat1{1},'LineWidth',lws,'MarkerSize',10);
plot(t,xd(:,1)/csf,'s-','Color',cmat1{2},'LineWidth',lws,'MarkerSize',10);
plot(t,xt(:,1)/csf,'s-','Color',cmat1{3},'LineWidth',lws,'MarkerSize',10);
%title('mSpinach gene expression versus time','FontSize',24)
h_leg = legend('Conv','Div','Tand');
set(h_leg,'Location','Best')
set(h_leg,'EdgeColor','White','FontSize',24)
xlabel('Time','FontSize',24)
ylabel('Expression (nM)','FontSize',24)
set(gca,'FontSize',24)
hold off
cmat2 = {[.9 0 0 ],[.7 .7 .7],[.5 .65 .65]};

figure(fig2)
hold on
plot(t,xc(:,2)/csf,'s-','Color',cmat2{1},'LineWidth',lws);

plot(t,xd(:,2)/csf,'s-','Color',cmat2{2},'LineWidth',lws);
plot(t,xt(:,2)/csf,'s-','Color',cmat2{3},'LineWidth',lws);
h_leg = legend('Conv','Div','Tand');
set(h_leg,'Location','Best')
set(h_leg,'EdgeColor','White','FontSize',24)
xlabel('Time','FontSize',24)
ylabel('Expression (nM)','FontSize',24)
set(gca,'FontSize',24)
hold off
% 

cd .. 
