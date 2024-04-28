function [a_G,t,SpecGenAcc] = Evol_Spec(nrecs,Sa_var,N,omega,dom,T,dt,zeta,N_iter,ag,Mw,Rjb,FN,FR,Vs30)

% Generate ground motions with the evolutionary spectrum based variant

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

%% Target Spectrum Variability
SpecGenAcc = zeros(N,nrecs);

a_G = cell(nrecs,1);
t = cell(nrecs,1);

t1=zeros(nrecs,1);
t2=zeros(nrecs,1);
Ts=zeros(nrecs,1);
tf=zeros(nrecs,1);

for jj=1:nrecs
     %% Envelope function creation
    
    [t1(jj),t2(jj),Ts(jj),tf(jj)]=RandDur(Mw,Rjb,FN,FR,Vs30);
       
    t{jj} = (0:dt:tf(jj)); % Acceleration Time step
    
    temp_t = t{jj,1};
    temp_t=temp_t';
     
    % Time history envelope function (Jennings et al.)
    tenv = [0, t1(jj), t2(jj), tf(jj)];
    beta_env = 3/(tenv(end)-tenv(end-1));
    [Env_t] = THmodulating(temp_t,tenv,beta_env);
   
    %% Target spectrum
    [G_EPSD,A_wt] = Spec_Compat_Accelerogram(Sa_var(jj,:),zeta,N,omega,dom,temp_t,Env_t,t1(jj),Ts(jj));
   
    %%  Accelerogram generation           
    [temp_aG] = AccGen(A_wt,G_EPSD,dom,omega,temp_t,N); 

    % Baseline correction
    [temp_aG] = baseline_correction(temp_aG','Q');
    
    %% Spectrum matching iterations
    [temp_aG,PSa_aG]=FT_iter([temp_t'; temp_aG],N_iter,[T';Sa_var(jj,:)],T',zeta,ag(jj));
    a_G{jj} = temp_aG;
    SpecGenAcc(:,jj)= PSa_aG(N_iter,:)';
end

