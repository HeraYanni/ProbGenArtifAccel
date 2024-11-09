clc; clearvars; close all;
set(groot,'defaultFigureColor','w')

% This code generates target spectrum-compatible fully-nonstationary 
% artificial seismic ground motions 
% The produced suites of ground motions can match a given target mean 
% spectrum and target variability for the whole period range of interest.

% works also for zero variability, i.e. only mean matching

% Further details are provided in the following document:
%
% Yanni H., Fragiadakis M., and Mitseas I.P. 
% "Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records". 
% Earthquake Engineering and Structural Dynamics.
% doi: https://doi.org/10.1002/eqe.4153

% Version 1.0 created by Hera Yanni, first release: 28th of April, 2024 

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nrecs = number of ground motions to be generated
% f0 = lowest frequency (Hz) of the frequency range of the generated signals 
% fmax = highest frequency (Hz) of the frequency range of the generated signals
% N_iter = number of corrective iterations
% dt = time-step of the generated signals
% N = number of harmonics to be superposed
% targSpec = type of target spectrum to be used
% T_corr = desired correlation between pairs of periods
% variant = variant used for the generation of the ground motion suites

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of User inouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the number of ground motions to be generated
nrecs = 7; % should be >1

% Define the desired frequency range of the ground motions
f0 = 1/(2*pi); % lowest frequency (Hz), can't be less than 0.36/(2*pi)
fmax = 100/(2*pi); % highest frequency (Hz)

% Define the desired number of corrective iterations
% range from 1 to 12 is ok 
% for perfect correlation high values produce good results, because the
% produced spectra are "smooth"
% for a specified correlation structure a value around 5 is recommended
% because the produced spectra are not "smooth"
N_iter = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the frequency/period range for target spectrum generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81; %m/sec2
dt = 0.01; %sec, signal time-step

Fs = 1/dt;   % sampling frequency (Hz)
Fn = Fs/2;   % nyquist frequency  (Hz)
omega_max = 2*pi*Fn; % maximum frequency
omega_m = 2*pi*fmax; % rad/sec

omega_u = min(2*pi*fmax, omega_max); % rad/sec
omega_0 = 2*pi*f0; % rad/sec

% Cut-off frequency upper limit check
if omega_m > omega_max
    msgbox('Maximum frequency must be less than the Nyquist frequency')
      return
end

% Cut-off frequency lower limit check
if omega_0 < 0.36 
    msgbox('Minimum frequency must be larger than 0.36 rad/s')
      return
end

N = 1000; % number of harmonics to be superposed
omega = linspace(omega_0,omega_u,N)'; % frequencies for the harmonics
dom = mean(diff(omega)); % frequency domain integration step
T = flip(2*pi./omega); % periods for the target spectrum spectral values 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the target spectrum, and target variability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GMM
% The median and logarithmic standard deviation values are obtained from 
% the BSSA14 NGA-West2 model, based on the following:
%
% Boore DM, Stewart JP, Seyhan E, Atkinson GM. 
% NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for 
% Shallow Crustal Earthquakes. Earthquake Spectra. 2014;30:1057-1085.
% https://doi.org/10.1193/070113EQS184M.

% The significant duration GMM is the D5%-95% model of the following:
%
% Sandıkkaya MA, Akkar S. Cumulative absolute velocity, Arias intensity 
% and significant duration predictive models from a pan-European
% strong-motion dataset. Bulletin of Earthquake Engineering. 
% 2017;15:1881-1898. https://doi.org/10.1007/s10518-016-0066-6.

% The correlation between pairs of periods are obtained from the following:
%
% Baker JW, Jayaram N. Correlation of Spectral Acceleration Values from 
% NGA Ground Motion Models. Earthquake Spectra. 2008;24(1):299-317.
% https://doi.org/10.1193/1.2857544

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EC8
% EC8 (EN1998-1-2004) code spectrum 
% works only with user defined coefficient of variation (CoV),
% the site and spectrum based variant, and perfect positive correlation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For targSpec = 'GMM', the model uses a GMM to generate the spectra and
% the accelerograms
% For targSpec = 'EC8', the model uses the EC8 spectrum to generate the 
% spectra and the accelerograms
targSpec = 'GMM';
%targSpec = 'EC8';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(targSpec,'GMM')

% GMM BSSA 14 model for the median and standard deviation values

Mw = 6.5; % Moment Magnitude
Rjb = 10; % Joyner-Boore distance (km)
% Fault_Type    = 0 for unspecified fault
%               = 1 for strike-slip fault
%               = 2 for normal fault
%               = 3 for reverse fault
Fault_Type = 2; 

Vs30 =450; % shear wave velocity averaged over top 30 m in m/s

% region        = 0 for global (incl. Taiwan)
%               = 1 for California
%               = 2 for Japan
%               = 3 for China or Turkey
%               = 4 for Italy
region = 4;

% z1            = Basin depth (km); depth from the groundsurface to the
%                 1km/s shear-wave horizon.
%               = 999 if unknown
z1 = 999;

% Peak groung acceleration
[ag sigma_ag] = gmpe_bssa_2014(Mw, 0, Rjb, Fault_Type, region, z1, Vs30);

ag = ag*g;
zeta = 0.05;

% Ground Motion Model median and standard deviation values (g)
for i=1:N
    [Sa_mean(i) sigma_ln(i)] = gmpe_bssa_2014(Mw, T(i), Rjb, Fault_Type, region, z1, Vs30);
end

Sa_mean = Sa_mean.*g; % values in m/s2

% Fault type parameters for the significant duration GMM
% Akkar and Sandikaya
FN=1; % =1 for normal fault, 0 otherwise
FR=0; % =1 for reverse fault, 0 otherwise

% Define the desired correlation between pairs of periods
%T_corr       = 1 for perfect positive correlation, i.e. rho(Ti,Tj)=1
%T_corr       = 2 for correlation rho(Ti,Tj) obtained from the Baker and
                % Jayaram 2008 model
T_corr = 1;
 
% weights=[ w1 w2 w3 ]        
weights =  [1 2.5 0.1]; 

elseif strcmp(targSpec,'EC8')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Eurocode 8 elastic spectrum (EN1998-1-2004)
    PGA = 0.24;         % peak ground acceleration in rock bed, in g
    gamma_I=1.00;       % importance factor
    ag = PGA*gamma_I*g; % peak ground acceleration, in m/s^2
    soil = 'B';         % EC8 soil type
    zeta= 0.05;         % damping ratio

    % Target mean spectrum
    [Sa_mean,~,~] = EC8spectrumElastic(ag,soil,T,zeta); % m/s^2

    % Target variability 
    CoV_targ = 0.3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the variant used for the generation of the ground motion suites
% Only for the GMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variant = 'site' for the site and spectrum based variant
% variant = 'evol' for the evolutionary spectrum based variant
 variant = 'evol';
% variant = 'site';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user inputs. 
% Calculations
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(targSpec,'GMM')
    %% Generate the target spectra
    if T_corr == 1      % perfect positive correlation
        % GMM 
        [Sa_var,a,pga] = Targ_Spec_GMM_rho1(nrecs,Sa_mean,sigma_ln,ag,sigma_ag,N);
        
    elseif T_corr == 2  % correlations as Baker and Jayaram 2008
        % GMM
        [Sa_var,a,pga] = Targ_Spec_GMM_rho(nrecs,Sa_mean,sigma_ln,ag,sigma_ag,N,T,weights);
    
    end

%% Generate the ground motions

if strcmp(variant,'evol')

    [a_G,t,SpecGenAcc] = Evol_Spec(nrecs,Sa_var,N,omega,dom,T,dt,zeta,N_iter,pga,Mw,Rjb,FN,FR,Vs30);

elseif strcmp(variant,'site')

    % Put the seed records files in .dat format and
    % in the same file as the matlab scripts
    % the file's first column should be the time and the second column the
    % acceleration values

    % the seed records time-step must be the same as the dt signal 
    % time-step in line 44
    
    % Put the seed record convertion factor in order to convert the record 
    % values to m/s2
    % for example if records are in g, put conv_factor=9.81
    conv_factor = 9.81;

    [a_G,t,SpecGenAcc] = Site_Spec(nrecs,conv_factor,Sa_var,N,omega,dom,T,zeta,N_iter,ag);    
end

%% Compute the spectral ordinates

SpecGenAcc_G = log(SpecGenAcc./g);

% Mean Spectrum of the generated accelerograms
MeanAccSpec = mean(SpecGenAcc_G,2);

% Standard deviation
sigma = std(SpecGenAcc_G,0,2);

medianSa = median(SpecGenAcc,2);
medianSa_var = median(Sa_var,1);

%% Final results error values
% Mean square error

MSE_mean=(sum(abs(medianSa'- medianSa_var).^2)/length(medianSa))
MSE_sigma=(sum(abs(sigma'-sigma_ln).^2)/length(sigma))

%% Outputs

% Save the target spectrum and variability
% save('Sa_mean.mat','Sa_mean');
% save('sigma_ln.mat','sigma_ln');

% Save the produced target spectra
% save('Sa_var.mat','Sa_var'); % Save the produced target spectra
% save('a.mat','a');           % Save the random variables a

% Save the generated accelerograms
% save('a_G.mat','a_G');       % Save the produced accelerograms
% save('t.mat','t');           % Save the time history

% Plot the produced target spectra
figure()
hold on;  box on; grid on;
for ii=1:nrecs
    % plot(T,Sa_var(ii,:)./g,'Color',[0.5 0.5 0.5]); for grey color
    plot(T,Sa_var(ii,:)./g); 
    plot(T,SpecGenAcc(:,ii)./g,'--'); 
end
p1 = plot(T,medianSa_var./g,'k','LineWidth',2);
plot(T,Sa_mean./g,'r--','LineWidth',3);
plot(T,medianSa./g,'b--','LineWidth',2);
set(gca,'FontSize',18,'FontName','Times New Roman')
xlim([0 4])
ylim([0 2])
xlabel('Period [s]','FontSize',20,'interpreter','latex');
ylabel('Spectral accelerations [g]','FontSize',20,'interpreter','latex');
legend(p1,{'${S^*_a}(T_i,\zeta)$'},'FontSize',18,'interpreter','latex','Location','southwest')
set(gca, 'XScale','log')
set(gca, 'YScale','log')
%      print('MedianTarget.eps','-depsc')
%      print('MedianTarget.jpg','-djpeg')

% Plot the mean/median comparison
figure()
hold on; grid on; box on;
plot(T,Sa_mean./g,'k','LineWidth',2)
plot(T,medianSa./g,'-.','Color',[0.4 0.4 0.4],'LineWidth',4)
set(gca,'FontSize',18,'FontName','Times New Roman')
xlabel('Period [s]','FontSize',20,'interpreter','latex');
ylabel('Spectral accelerations [g]','FontSize',20,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','southwest')
set(gca, 'XScale','log')
set(gca, 'YScale','log')
xlim([0 4])
%ylim([0 2])
%     print('MedianGMMEPSD.eps','-depsc')
%     print('MedianGMMEPSD.jpg','-djpeg')

% Plot the variability comparison
figure()
hold on; grid on; box on;
set(gca,'FontSize',18,'FontName','Times New Roman')
plot(T,sigma_ln,'k','LineWidth',2)
plot(T,sigma,'Color',[0.4 0.4 0.4],'LineWidth',3)
xlabel('Period [s]','FontSize',20,'interpreter','latex');
ylabel('Variability $\beta$','FontSize',20,'interpreter','latex');
legend('$\beta^*(T_i,\zeta)$','$\beta(T_i,\zeta)$','FontSize',17,'interpreter','latex','Location','southeast')   %'Inherent $\sigma_{ln}$'
set(gca, 'XScale','log')
set(gca, 'YScale','log')
%ylim([0.55 0.8])
xlim([0 4])
%     print('CoVGMMEPSD.eps','-depsc')
%     print('CoVGMMEPSD.jpg','-djpeg')

% plot the accelerograms
for i=1:nrecs
    figure()
hold on;  box on; grid on;
plot(t{i},a_G{i}./g,'k')
xlim([0,t{i}(end)]);
ylim([-0.41 0.41]);    
set(gca,'FontSize',18,'FontName','Times New Roman')
xlabel('Time [s]','FontSize',20,'interpreter','latex');
ylabel('$\alpha(t)$ [g]','FontSize',20,'interpreter','latex');
yline(0,'k')
end

elseif strcmp(targSpec,'EC8')
    [Sa_var,a,pga] = Targ_Spec_EC8(nrecs,Sa_mean,CoV_targ,ag,N);

    % Put the seed records files in .dat format and
    % in the same file as the matlab scripts
    % the file's first column should be the time and the second column the
    % acceleration values

    % the seed records time-step must be the same as the dt signal 
    % time-step in line 44
    
    % Put the seed record convertion factor in order to convert the record 
    % values to m/s2
    % for example if records are in g, put conv_factor=9.81
    conv_factor = 9.81;

    [a_G,t,SpecGenAcc] = Site_Spec(nrecs,conv_factor,Sa_var,N,omega,dom,T,zeta,N_iter,ag); 

 %% Compute the spectral ordinates

 % Mean Spectrum of the generated accelerograms
MeanAccSpec2 = mean(SpecGenAcc,2);

% Standart deviation
sigma = std(SpecGenAcc,0,2);

% Coefficient of Variation
CoV = sigma./MeanAccSpec2;

meanSa_var = mean(Sa_var,1);

%% Final results error values
% Mean square error

MSE_mean=(sum(abs(MeanAccSpec2'- meanSa_var).^2)/length(MeanAccSpec2))
MSE_sigma=(sum(abs(CoV'-CoV_targ).^2)/length(CoV_targ))

%% Outputs

% Save the target spectrum and variability
% save('Sa_mean.mat','Sa_mean');
% save('CoV_targ.mat','CoV_targ');

% Save the produced target spectra
% save('Sa_var.mat','Sa_var'); % Save the produced target spectra
% save('a.mat','a');           % Save the random variables a

% Save the generated accelerograms
% save('a_G.mat','a_G');       % Save the produced accelerograms
% save('t.mat','t');           % Save the time history

% Plot the produced target spectra
figure()
hold on;  box on; grid on;
for ii=1:nrecs
    % plot(T,Sa_var(ii,:)./g,'Color',[0.5 0.5 0.5]); for grey color
    plot(T,Sa_var(ii,:)./g); 
    plot(T,SpecGenAcc(:,ii)./g,'--'); 
end
p1 = plot(T,meanSa_var./g,'k','LineWidth',2);
plot(T,Sa_mean./g,'r--','LineWidth',3);
plot(T,MeanAccSpec2./g,'b--','LineWidth',2);
set(gca,'FontSize',18,'FontName','Times New Roman')
xlim([0 4])
% ylim([0 2])
xlabel('Period [s]','FontSize',20,'interpreter','latex');
ylabel('Spectral accelerations [g]','FontSize',20,'interpreter','latex');
legend(p1,{'${S^*_a}(T_i,\zeta)$'},'FontSize',18,'interpreter','latex','Location','southwest')
%      print('MeanTarget.eps','-depsc')
%      print('MeanTarget.jpg','-djpeg')

% Plot the mean/median comparison
figure()
hold on; grid on; box on;
plot(T,Sa_mean./g,'k','LineWidth',2)
plot(T,MeanAccSpec2./g,'-.','Color',[0.4 0.4 0.4],'LineWidth',4)
set(gca,'FontSize',18,'FontName','Times New Roman')
xlabel('Period [s]','FontSize',20,'interpreter','latex');
ylabel('Spectral accelerations [g]','FontSize',20,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','southwest')
xlim([0 4])
%ylim([0 2])
%     print('MeanEC8.eps','-depsc')
%     print('MeanEC8.jpg','-djpeg')

% Plot the variability comparison
figure()
hold on; grid on; box on;
set(gca,'FontSize',18,'FontName','Times New Roman')
yline(CoV_targ,'k','LineWidth',2)
plot(T,CoV,'Color',[0.4 0.4 0.4],'LineWidth',3)
xlabel('Period [s]','FontSize',20,'interpreter','latex');
ylabel('Variability $\beta$','FontSize',20,'interpreter','latex');
legend('$\beta^*(T_i,\zeta)$','$\beta(T_i,\zeta)$','FontSize',17,'interpreter','latex','Location','southeast')   %'Inherent $\sigma_{ln}$'
xlim([0 4])
%     print('CoVEC8.eps','-depsc')
%     print('CoVEC8.jpg','-djpeg')

% plot the accelerograms
for i=1:nrecs
    figure()
hold on;  box on; grid on;
plot(t{i},a_G{i}./g,'k')
xlim([0,t{i}(end)]);
ylim([-0.41 0.41]);    
set(gca,'FontSize',18,'FontName','Times New Roman')
xlabel('Time [s]','FontSize',20,'interpreter','latex');
ylabel('$\alpha(t)$ [g]','FontSize',20,'interpreter','latex');
yline(0,'k')
end
end
