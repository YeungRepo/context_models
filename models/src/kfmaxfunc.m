function[kfX] = kfmaxfunc(XLength,sigmaX,s0,kfmax,std_coil,pLen)



kfX = std_coil*kfmax*1/(std_coil+(sigmaX - s0*pLen/XLength)^2);

