function[kcatX] = kcatmaxfunc(t,XLength,sigmaX,s0,kcatmax,std_coil,pLen)

kcatX = std_coil*kcatmax*1/(std_coil+(sigmaX - s0*pLen/XLength)^2);
