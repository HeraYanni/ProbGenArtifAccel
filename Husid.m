% Husid function computation 
% time point t1 is the beginning of the strong motion phase at 5%  
% time point t2 is the end of the strong motion phase at 95%

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

function [Hus,t1,t2] = Husid(t,acc)

% Husid plot
Hus = zeros(length(t),1); 
Ht = cumsum(acc.^2);

for i=1:length(t)
   Hus(i) = Ht(i)/Ht(end);
end    

% Time values were H(t1) = 5% , H(t2)=95%

h1 = abs(Hus-0.05);
h2 = abs(Hus-0.95);

[~,i1] = min(h1);
[~,i2] = min(h2);

t1=t(i1);
t2=t(i2);

end