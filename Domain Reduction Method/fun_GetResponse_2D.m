function [U,V] = fun_GetResponse_2D(t,coords,angle,Vs,poi,Node_0,Surface_z)

global Disp Velo Acce

Vp = Vs*sqrt(2*(1-poi)/(1-2*poi));
theta_s = angle/180.0*pi;
theta_p = asin(Vp/Vs*sin(theta_s));
      
x0 = Node_0(1) + (Surface_z-Node_0(2)) * tan(theta_s);
y0 = Surface_z;

x_rela = coords(1) - x0;
y_rela = y0 - coords(2);  

t_ini = (Surface_z-Node_0(2))/Vs/cos(theta_s);

k = (Vp/Vs);
      
U_si = 1.0;
U_sr = (sin(2*theta_s)*sin(2*theta_p)-k^2*cos(2*theta_s)^2)/...
       (sin(2*theta_s)*sin(2*theta_p)+k^2*cos(2*theta_s)^2)*U_si;
U_pr = -(-2*k^2*sin(2*theta_s)*cos(2*theta_s))/...
       (sin(2*theta_s)*sin(2*theta_p)+k^2*cos(2*theta_s)^2)*U_si*(Vs/Vp);
  
t1 = -x_rela/Vs*sin(theta_s)+y_rela/Vs*cos(theta_s)+t-t_ini;
t2 = -x_rela/Vs*sin(theta_s)-y_rela/Vs*cos(theta_s)+t-t_ini;
t3 = -x_rela/Vp*sin(theta_p)-y_rela/Vp*cos(theta_p)+t-t_ini;


U(1,:)  = U_si*cos(theta_s)*interp1(t,Disp,t1,'pchip',0) + ...
        -U_sr*cos(theta_s)*interp1(t,Disp,t2,'pchip',0) + ...
        U_pr*sin(theta_p)*interp1(t,Disp,t3,'pchip',0);
    
U(2,:)  = U_si*cos(theta_s)*interp1(t,Velo,t1,'pchip',0) + ...
        -U_sr*cos(theta_s)*interp1(t,Velo,t2,'pchip',0) + ...
        U_pr*sin(theta_p)*interp1(t,Velo,t3,'pchip',0);
    
U(3,:)  = U_si*cos(theta_s)*interp1(t,Acce,t1,'pchip',0) + ...
        -U_sr*cos(theta_s)*interp1(t,Acce,t2,'pchip',0) + ...
        U_pr*sin(theta_p)*interp1(t,Acce,t3,'pchip',0);    
    
V(1,:)  = U_si*sin(theta_s)*interp1(t,Disp,t1,'pchip',0) + ...
          U_sr*sin(theta_s)*interp1(t,Disp,t2,'pchip',0) + ...
          U_pr*cos(theta_p)*interp1(t,Disp,t3,'pchip',0);
    
V(2,:)  = U_si*sin(theta_s)*interp1(t,Velo,t1,'pchip',0) + ...
          U_sr*sin(theta_s)*interp1(t,Velo,t2,'pchip',0) + ...
          U_pr*cos(theta_p)*interp1(t,Velo,t3,'pchip',0);
    
V(3,:)  = U_si*sin(theta_s)*interp1(t,Acce,t1,'pchip',0) + ...
          U_sr*sin(theta_s)*interp1(t,Acce,t2,'pchip',0) + ...
          U_pr*cos(theta_p)*interp1(t,Acce,t3,'pchip',0);       

V = -V;

end

