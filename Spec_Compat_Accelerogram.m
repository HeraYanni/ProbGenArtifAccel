% Evolutionary power spectral density method
%
% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

function [G_EPSD,A_wt] = Spec_Compat_Accelerogram(Sa,zeta,N,omega,dom,t,Env_t,t1,Ts)
g = 9.81; %m/sec2

%% Stationary One sided Power Spectral Density (PSD) calculation
Sa = flip(Sa); %m/sec2
[Gstat] = PSDstationary(zeta,N,omega,Sa,Ts,dom);

%% Evolutionary PSD: Non-separable process 
% Modulating function constants
a_mod = [2/(t1+Ts/2), 0.01 , 0.0, 2]; % [p0, p1, p2, gamma]
[G_EPSD,A_wt] = EPSD(a_mod,omega,t,N,Gstat,Env_t);

end
