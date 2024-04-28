% Code that computes the response spectrum of accelerograms

% Copyright (cc) to Ferreira F, Moutinho C, Cunha Á, Caetano E. 
% An artificial accelerogram generator code written in Matlab. 
% Engineering Reports. 2020;2:e12129. https://doi.org/10.1002/eng2.12129
% Open Access

% DETERMINE THE ACCELEROGRAM RESPONSE SPECTRA (ARS)
function [S]=ARS(accelerogram,z,T)
% Function to determine the linear elastic spectra of a SDOF oscillator 
% INPUT
% T is the oscillator periods vector (s)
% z is the oscillator damping
% accelerogram is the input excitation: [time; acceleration]
% OUTPUT
% S is the response spectra: in the form: 
% S= [period; (displacement, velocity and acceleration response spectra)]
S=zeros(4,length(T));
S(1,:)=T;
w=1./T*2*pi;
%Total_time=accelerogram(1,end);
dt=accelerogram(1,2)-accelerogram(1,1);
time=accelerogram(1,:);
Forces=interp1(accelerogram(1,:),accelerogram(2,:),time);
dForces_dt=Forces(2:end)-Forces(1:(end-1));
for jj=1:length(w); %parfor using parallel computing toolbox
% State Space Differential Equations Matrixes
A=[ 0         1;
    -w(jj)^2, -2*z*w(jj)];
B=[ 0; 1];
% Analytical Time Step Solver of State Space 
Kb=real(exp(1)^(A*dt));
Kf=A\(Kb-eye(2));
KfB=Kf*B;
KdfB=A\(Kf/dt-eye(2))*B;
% Diffenrential Equation Solver
x=zeros(2,length(time)); %1st line -displacement; 2nd line - velocity
acel=zeros(1,length(time)); % acceleration 
aux_acel=-[w(jj)^2,2*z*w(jj)]; 
for j=1:(length(time)-1)   
   x(:,j+1)=Kb*x(:,j)+KfB*Forces(1,j)+KdfB*dForces_dt(1,j);
   acel(j+1)=aux_acel*x(:,j+1); %acceleration
end
%Displacement, velocity and acceleration response spectra
S(2:4,jj)=[max(abs(x(1,:)));max(abs(x(2,:)));max(abs(acel))]; 
end
end