function [H] = HeavisideStep(x)
% Heaviside step functtion

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

H=zeros(length(x),1);
for i=1:length(x)
if x(i)>0
    H(i)=1;
else
    H(i)=0;
end
end
end