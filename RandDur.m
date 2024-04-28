function [t1,t2,Ts,tf]=RandDur(Mw,Rjb,FN,FR,Vs30)

% Generate random temporal characteristics based on a significant duration GMM
% t1 = the beginning of the strong motion phase
% Ts = significant duration
% tf = total duration

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

% start of strong motion
t1=3+2*rand;

% Strong motion duration
[Ts] = StrongMotionDur(Mw,Rjb,FN,FR,Vs30);

% End of strong motion
t2 = t1+Ts;

% Total duration
tf = t2+30;
end