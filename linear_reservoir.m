% function to compute runoff from upper and lower reservoir
% S_1 storage  reservoir
% d_t time step
% k1 upper reservoir constant lower outlet
% 

function [Q_1]=linear_reservoir(S_1,d_t,k1);
% Runoff from linear reservoir

Q_1=k1*S_1*d_t;

