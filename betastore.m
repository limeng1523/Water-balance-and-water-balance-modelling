% function to compute direct runoff  Qdir
% Input Prec Precipitation
% d_t time step
% S_M soil moisture
% S_max maximum soil moisture
% b_eta beta parameter 
function Qdir=betastore(Prec,d_t,S_M,S_max,b_eta);

Qdir=Prec*d_t*(S_M/S_max)^b_eta;
