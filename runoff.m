% function to compute runoff from upper and lower reservoir
% S_1 storage upper reservoir
% S_2 storage lower reservoir
% G_s percolation from soil store
% d_t time step
% k0 upper reservoir constant upper outlet
% k1 upper reservoir constant lower outlet
% k2 lower reservoir constant 
% kperc upper reservoir constant for percolation to groundwater store 
% [Q0(i) Q1(i) Qperc(i) Q2(i)]
% (G(i),S1(i),S2(i),dt,k0,k1,kperc,k2,L)


function [Q_0, Q_1, Q_perc, Q_2]=runoff(G_s,S_1,S_2,d_t,k0,k1,kperc,k2,L);
% Runoff from upper reservoir
if S_1 > L
      Q_0=k0*(S_1-L)*d_t;
else 
      Q_0=0.;
end
% Runoff from upper reservoir
Q_1=k1*S_1*d_t;

% Percolation into lower reservoir (from upper reservoir and Groundwater
% recharge G
Q_perc=kperc*S_1*d_t+G_s;

% Runoff from lower reservoir
Q_2=k2*S_2*d_t;

