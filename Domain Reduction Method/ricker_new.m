function [Acce,Velo,Disp] = ricker_new(fricker,Aricker,t,t0,dt)


for i = 1:length(t) 
    tau = t(i)-t0;
    Disp(i,1) = Aricker*(1.0-2.0*tau.*tau*fricker*fricker*pi*pi).*exp(-tau.*tau*pi*pi*fricker*fricker);
end
    
% Velo = cumtrapz(t,Acce);
% Disp = cumtrapz(t,Velo);

Velo = [0; diff(Disp)/dt];
Acce = [0; diff(Velo)/dt];
end







