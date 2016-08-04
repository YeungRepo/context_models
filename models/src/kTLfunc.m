function[kTLeff] = kTLfunc(t,kTL)

nhours = 2;
if (t>=60*60*nhours)
   half_life = 80*60;
   alpha_decay = log(2)/half_life;
   ksc = exp(-alpha_decay*(t-nhours*60*60)); 
   kTLeff = ksc*kTL;
else
    kTLeff = kTL;
end
