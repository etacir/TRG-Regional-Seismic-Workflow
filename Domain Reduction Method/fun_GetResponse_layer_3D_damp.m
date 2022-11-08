function [U,V,W] = fun_GetResponse_layer_3D_damp(t,t_int,coords,angle_h,angle_v,Node_0,Surface_z,Layer_info,u,w,Ralpha,Rbeta)

T = [ cosd(angle_h)  sind(angle_h) 0;
     -sind(angle_h)  cosd(angle_h) 0;
                  0              0 1];
Node_0 = T*Node_0';
Node_0 = [Node_0(1) Node_0(3)];

coords = T*coords';
coords = [coords(1) coords(3)];


nu = Layer_info(:,4);
G = Layer_info(:,2).*Layer_info(:,3).^2;
E = G*2.*(1+nu);
L = E.*nu./(1+nu)./(1-2*nu);
Vs = Layer_info(:,3);
Vp = sqrt(2*(1-nu)./(1-2*nu)).*Vs;
h = Layer_info(:,1);
n_Layer = size(Layer_info,1);

t0 = (Surface_z - sum(h(1:end-1)) - Node_0(2))*cosd(angle_v)/Vs(end);

for i = 1:n_Layer
    if i == 1
        Layer_vertical_coord(i) = Surface_z;
    else
        Layer_vertical_coord(i) = Layer_vertical_coord(i-1) - (h(i-1));
    end 
end

for i = n_Layer:-1:1
    
    if coords(2) <= Layer_vertical_coord(i)
        Layer_ID = i;
        break
    end
end

Beta = real(asind(Vs(Layer_ID)/Vs(end)*sind(angle_v)));
Alpha = real(asind(Vp(Layer_ID)/Vs(Layer_ID)*sind(Beta)));

x0 = coords(1) - Node_0(1);

if Layer_ID == n_Layer
    y0 = Layer_vertical_coord(Layer_ID) - coords(2);
else
    y0 = Layer_vertical_coord(Layer_ID+1) - coords(2);
end

n = length(t);
dt = t(2) - t(1);

f = linspace(0,1,n)/dt;

for i = 1:n
    n_Acce = i;
    omega = 2*pi*f(n_Acce);
        
    xi = 1/2*(Ralpha/omega + Rbeta*omega);
    if xi >= 0.2
        xi = 0.2;
    end
    
    if omega == 0
        xi = 0;
    end
            
    
    kx = omega*sind(Beta)/(Vs(Layer_ID)*sqrt(1+2*1i*xi));
    kpz = omega*cosd(Alpha)/(Vp(Layer_ID)*sqrt(1+2*1i*xi));
    ksz = omega*cosd(Beta)/(Vs(Layer_ID)*sqrt(1+2*1i*xi));
    
    amp1(i) = exp(1j*(-imag(kx)*1j)*x0)*exp(1j*(imag(kpz)*1j)*y0);
    amp2(i) = exp(1j*(-imag(kx)*1j)*x0)*exp(1j*(-imag(kpz)*1j)*y0);
    amp3(i) = exp(1j*(-imag(kx)*1j)*x0)*exp(1j*(imag(ksz)*1j)*y0);
    amp4(i) = exp(1j*(-imag(kx)*1j)*x0)*exp(1j*(-imag(ksz)*1j)*y0);

    Kx(i) = kx;
    Kpz(i) = kpz;
    Ksz(i) = ksz;
    

end


Beta = real(asind(Vs(Layer_ID)/Vs(end)*sind(angle_v)));
Alpha = real(asind(Vp(Layer_ID)/Vs(Layer_ID)*sind(Beta)));

u{Layer_ID}(1,:) = u{Layer_ID}(1,:).*amp1;
u{Layer_ID}(2,:) = u{Layer_ID}(2,:).*amp2;
u{Layer_ID}(3,:) = u{Layer_ID}(3,:).*amp3;
u{Layer_ID}(4,:) = u{Layer_ID}(4,:).*amp4;

w{Layer_ID}(1,:) = w{Layer_ID}(1,:).*amp1;
w{Layer_ID}(2,:) = w{Layer_ID}(2,:).*amp2;
w{Layer_ID}(3,:) = w{Layer_ID}(3,:).*amp3;
w{Layer_ID}(4,:) = w{Layer_ID}(4,:).*amp4;


for i = size(Layer_info,1):-1:1
    ut{i} = ifft(u{i},n,2,'symmetric');
    wt{i} = ifft(w{i},n,2,'symmetric');
end

Velo_info_X = ut;
Velo_info_Y = wt;

Velo_info_X = Velo_info_X{Layer_ID};
Velo_info_Y = Velo_info_Y{Layer_ID};

t1 = -x0/Vp(Layer_ID)*sind(Alpha) + y0/Vp(Layer_ID)*cosd(Alpha) + t_int - t0;
t2 = -x0/Vp(Layer_ID)*sind(Alpha) - y0/Vp(Layer_ID)*cosd(Alpha) + t_int - t0;
t3 = -x0/Vs(Layer_ID)*sind(Beta) + y0/Vs(Layer_ID)*cosd(Beta) + t_int - t0;
t4 = -x0/Vs(Layer_ID)*sind(Beta) - y0/Vs(Layer_ID)*cosd(Beta) + t_int - t0;

VeloX = interp1(t,Velo_info_X(1,:),t1,'pchip',0) + ...
        interp1(t,Velo_info_X(2,:),t2,'pchip',0) + ...
        interp1(t,Velo_info_X(3,:),t3,'pchip',0) + ...
        interp1(t,Velo_info_X(4,:),t4,'pchip',0);
VeloY = interp1(t,Velo_info_Y(1,:),t1,'pchip',0) + ...
        interp1(t,Velo_info_Y(2,:),t2,'pchip',0) + ...
        interp1(t,Velo_info_Y(3,:),t3,'pchip',0) + ...
        interp1(t,Velo_info_Y(4,:),t4,'pchip',0);

U(1,:) = cumtrapz(t_int,VeloX);
U(2,:) = VeloX;
U(3,:) = [0 diff(VeloX)/dt];

V(1,:) = -cumtrapz(t_int,VeloY);
V(2,:) = -VeloY;
V(3,:) = -[0 diff(VeloY)/dt];

W = V;
V = U*sind(angle_h);
U = U*cosd(angle_h);


end



