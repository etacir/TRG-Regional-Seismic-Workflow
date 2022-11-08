function [ut,wt] = fun_GenerateAnalyticalSolu_SV_2D_Layer(t,dt,Velo,Layer_info,Angle)

n = length(t);
f = linspace(0,1,n)/dt;
Velo_FFT = fft(Velo,n);

nu = Layer_info(:,4);
G = Layer_info(:,2).*Layer_info(:,3).^2;
E = G*2.*(1+nu);
L = E.*nu./(1+nu)./(1-2*nu);
Vs = Layer_info(:,3);
Vs = Vs.*(1+2*1i*0);

Vp = sqrt(2*(1-nu)./(1-2*nu)).*Vs;
h = Layer_info(:,1);
n_Layer = size(Layer_info,1);

K = eye(4);

%% Computer Q matrix

for ii = 1:n
    n_Acce = ii;
    omega = 2*pi*f(n_Acce);
    
    K = eye(4);
    
    for i = size(Layer_info,1):-1:1

        if i == size(Layer_info,1)
            beta(i) = Angle;        
        else
            beta(i) = asind(Vs(i)/Vs(i+1)*sind(beta(i+1)));
        end
        alpha(i) = asind(Vp(i)/Vs(i)*sind(beta(i)));

        kx = omega*sind(beta(i))/Vs(i);

        kpz = omega*cosd(alpha(i))/Vp(i);
        ksz = omega*cosd(beta(i))/Vs(i);

        P = -L(i)*kx^2 - (L(i)+2*G(i))*kpz^2;
        eta = G(i)*(ksz^2-kx^2);

        Q{i} = [       -1j*kx*1j*omega         -1j*kx*1j*omega        -1j*ksz*1j*omega         1j*ksz*1j*omega;
                       1j*kpz*1j*omega        -1j*kpz*1j*omega         -1j*kx*1j*omega         -1j*kx*1j*omega;
                            P              P  2*G(i)*ksz*kx -2*G(i)*ksz*kx;
                2*G(i)*kpz*kx -2*G(i)*kpz*kx            eta            eta];  

        if i ~= size(Layer_info,1)
            d{i} = diag([exp(1j*kpz*h(i)), exp(-1j*kpz*h(i)) exp(1j*ksz*h(i)) exp(-1j*ksz*h(i))]);
            if omega == 0
                T{i} = eye(4);
            else
                T{i} = inv(Q{i+1})*Q{i}*d{i};
            end
            K = K*T{i};
        else
            if omega == 0
                Es_n1 = 0;
            else
                Es_n1 = -(Velo_FFT(n_Acce))*cosd(Angle)/ksz/1j/omega/1j;
            end

            Ep_n1 = 0;

        end

        if i == 1
            P1 = P;
            kx1 = kx;
            kpz1 = kpz;
            ksz1 = ksz;
        end
        
        Kx(i) = kx;
        Kpz(i) = kpz;
        Ksz(i) = ksz;

    end

    K1 = [             P1               P1     2*G(1)*ksz1*kx1    -2*G(1)*ksz1*kx1;
          2*G(1)*kpz1*kx1 -2*G(1)*kpz1*kx1 G(1)*(ksz1^2-kx1^2) G(1)*(ksz1^2-kx1^2)];

    K2 = [1 0 0 0;
          0 0 1 0];

    if omega == 0  
        E_1_n1 = zeros(8,1);
    else
        E_1_n1 = ([K diag([-1 -1 -1 -1]);K1 zeros(2,4); zeros(2,4) K2]\[zeros(4,1); zeros(2,1); 0; Es_n1]);
    end


    E1 = E_1_n1(1:4);
    En1 = E_1_n1(5:8);
    E_vec{1,ii} = E1;
    E_vec{n_Layer,ii} = En1;
    
    u{n_Layer}(:,ii) = [Q{n_Layer}(1,1)*En1(1); Q{n_Layer}(1,2)*En1(2); Q{n_Layer}(1,3)*En1(3); Q{n_Layer}(1,4)*En1(4)];
    w{n_Layer}(:,ii) = [Q{n_Layer}(2,1)*En1(1); Q{n_Layer}(2,2)*En1(2); Q{n_Layer}(2,3)*En1(3); Q{n_Layer}(2,4)*En1(4)];

    for i = 1:n_Layer-1
        if i >= 2
            E_vec{i,ii} = (T{i-1}*E_vec{i-1,ii});
        end
        u{i}(:,ii) = [Q{i}(1,1)*d{i}(1,1)*E_vec{i,ii}(1); Q{i}(1,2)*d{i}(2,2)*E_vec{i,ii}(2); Q{i}(1,3)*d{i}(3,3)*E_vec{i,ii}(3); Q{i}(1,4)*d{i}(4,4)*E_vec{i,ii}(4)];
        w{i}(:,ii) = [Q{i}(2,1)*d{i}(1,1)*E_vec{i,ii}(1); Q{i}(2,2)*d{i}(2,2)*E_vec{i,ii}(2); Q{i}(2,3)*d{i}(3,3)*E_vec{i,ii}(3); Q{i}(2,4)*d{i}(4,4)*E_vec{i,ii}(4)];
    end
    
    
end


for i = size(Layer_info,1):-1:1
    ut{i} = ifft(u{i},n,2,'symmetric');
    wt{i} = ifft(w{i},n,2,'symmetric');
end

end

