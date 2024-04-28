function [a_G,t,SpecGenAcc] = Site_Spec(nrecs,conv_factor,Sa_var,N,omega,dom,T,zeta,N_iter,ag)
% Generate ground motions with the site and spectrum based variant

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com

% Load the seed records
recData = dir('*.dat') ; 
NrecData = length(recData) ; 

%% Accelerogram generation
SpecGenAcc = zeros(N,nrecs);
a_G = cell(nrecs,1);
t = cell(nrecs,1);


for jj=1:nrecs
        %% Seed record selection
    load_rec_acc = recData(randi(NrecData)).name; 
        
    % Recorded accelerogram file:
    % first column = time values, second column = acceleration values
    % time in sec, acceleration in g
    rec_aG = load(load_rec_acc);

    % Acceleration Time
    temp_t = rec_aG(:,1);
    t{jj} = temp_t';
    rec_acc = rec_aG(:,2).*conv_factor;

    % Husid function and plot
    [Hus,t1,t2]=Husid(temp_t,rec_acc);

    %% Time history envelope function (Jennings et al.)
    tenv = [0, t1, t2, temp_t(end)];
    beta_env = 3/(tenv(end)-tenv(end-1));
    [Env_t] = THmodulating(temp_t,tenv,beta_env);
  
    %% Recorded Accelerogram response spectrum
    [Sa_rec]=ARS([temp_t'; rec_acc'],zeta,T');
    Sa_rec = Sa_rec(4,:)';
  
 
    %% Coefficient a
    a_coeff = min(Sa_var(jj,:)'./Sa_rec);
    if a_coeff>1
        a_coeff=1;
    end
    
    % PSD
    Sa1 = flip(Sa_var(jj,:)'); % m/sec2
    Sa2 = flip(Sa_rec);  % m/sec2
    Ts = tenv(3)-tenv(2); % strong motion phase duration
    [Gstat] = PSD_rec_artif(zeta,N,omega,Sa1,Sa2,Ts,dom,a_coeff);
    
    %%  Accelerogram generation
    [temp_aG] = AccGen_rec_artif(a_coeff,rec_acc,Env_t,Gstat,dom,omega,temp_t,N);
    
    % Baseline correction
    [temp_aG] = baseline_correction(temp_aG','Q');
  
    [temp_aG,PSa_aG]=FT_iter([temp_t'; temp_aG],N_iter,[T';Sa_var(jj,:)],T',zeta,ag);
    
    a_G{jj} = temp_aG;
    SpecGenAcc(:,jj)= PSa_aG(N_iter,:)';   
end

end
