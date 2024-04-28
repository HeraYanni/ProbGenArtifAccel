% Time modulating function of Jennings PC, Housner GW, Tsai C. 
% Simulated earthquake motions for design purpose. 
% In: Proceedings of the 4th world conference earth Eng Santiago, 
% vol. A-1; 1969. pp. 145–160.

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

function [Env_t] = THmodulating(t,tmod,b_env)
%% Modulating function


Env_t = zeros(1,length(t));

for ti = 1:length(t)
    
    if t(ti)< tmod(2)
        Env_t(ti) = (t(ti)/tmod(2))^2;

    elseif tmod(2)<=t(ti) && t(ti)<=tmod(3)
        Env_t(ti) = 1;
        
    else
        Env_t(ti) = exp(-b_env*(t(ti)-tmod(3)));
    end

end

% if iplot == 1 
%     figure()
%     hold on; grid on; box on;
%     plot(t,Env_t,'b','Linewidth',2)
%     title('Time history envelope function (Jennings et al.)');
%     xlabel('Time [sec]');
%     ylabel('\phi [t]');
%     xlim([0 t(end)])
% end

end