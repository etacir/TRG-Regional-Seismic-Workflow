function [u_history,udot_history,uddot_history] = Site_response_analysis_seismic_SV_wave

%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% This code generates required DRM information for doing scattered type
% analysis.
% The output is equivalent nodal forces along the DRM interface
%
% before running this code:
% 1) update the L_bottom in line 32
% 2) update shear wave velocity and density of the soil below the DRM 
%    interface in 35
% 3) update element size in line 38
% 4) update rayleigh damping coefficients in line 41-45
% 5) modify the incident wave function in line 70-86
% 6) update density of the soil under certain elevation in line 182, this
%    can be defined as a piecewise function by using the coordinates of
%    that integration points (x1, x2, x3)
% 7) update poisson's ratio of the soil under certain elevation in line 
%    184, this can be defined as a piecewise function by using the 
%    coordinates of that integration points (x1, x2, x3)
% 8) update shear wave velocity of the soil under certain elevation in line 
%    186, this can be defined as a piecewise function by using the 
%    coordinates of that integration points (x1, x2, x3)
%
%--------------------------------------------------------------------------
% this part should be given by the user
%--------------------------------------------------------------------------
% Elevation of the bottom of the soil column 
L_bottom = -20; 

% damshpot coefficient for the soil below the DRM interface
Vs_ER = 400; rho_ER = 1800; nu = 0.3;
damper_c_ER = rho_ER * Vs_ER;
V0 = Vs_ER; VH = Vs_ER; nvs = 0.5; b = (V0/VH)^(1/nvs);

% element size
El_size = 2;
%--------------------------------------------------------------------------
% defining rayleigh damping coefficients
xi1 = 0.05; xi2 = 0.05;
f1  = 3 ; f2  = 5*f1;
AC  =1/4/pi*[1/f1 4*pi^2*f1;1/f2 4*pi^2*f2];
aC = AC\[xi1;xi2];
aC = [0;0]; %undamped system

%--------------------------------------------------------------------------
%Gauss Quadrature order
QuadOrder = 3;
% number of nodes per elemens
NEN = 4;

%--------------------------------------------------------------------------
% reading mesh info for the soil column
%--------------------------------------------------------------------------
[EssentialBCX,EssentialBCY,NodeXYCoordinate, ElementConnectivity] = fun_MeshInfo_Soil_Column_small;

% node numbers
nnodes = size(NodeXYCoordinate,1); 
% element numbers
nelems = size(ElementConnectivity,1);

% total dof (ux, uy)
nGlDof = 2*nnodes;
%--------------------------------------------------------------------------
% defining the incident velocities along the lysmer boundary
% the following is for ricker pulse. it can be modified to read an
% arbitrary signal
%--------------------------------------------------------------------------
% Defining the global variables
global dt Disp Velo Acce

% some time info for the input motion
t1 = (length(Acce)-1)*dt;t = (0:dt:t1)';n = length(t);

%--------------------------------------------------------------------------%
% prescribing the input motions
Input_ug = Disp; Input_ugdot  = Velo; Input_ugddot = Acce;

%--------------------------------------------------------------------------
% reading gauss quadrature info
[GaussQuad_coef, GaussQuad_wgt] = fun_Order3GaussQuadDataOrder;
%--------------------------------------------------------------------------
% initialization of the mass, damping, and stiffness matrices
M = zeros(nGlDof); C = zeros(nGlDof); K = zeros(nGlDof); CKL = zeros(nGlDof);
%--------------------------------------------------------------------------
% computing elemental matrices and assembling
jC = 0;

for e = 1:nelems
    
    %disp(e)
    
    ElementDof = ElementConnectivity(e,2:end);

    ElementDof2GlDofu1  =            ElementDof;
    ElementDof2GlDofu2  =   nnodes + ElementDof;
    
    ElementXY = NodeXYCoordinate(ElementDof,2:end);
        
    %Mu1u1
    M_local_u1u1 = zeros(NEN);
    
    %Mu2u2
    M_local_u2u2 = zeros(NEN);    
    %--------------------------------------------------------------------------
    % C matrix
    %--------------------------------------------------------------------------
    % Cu1u1
    CKL_local_u1u1 = zeros(NEN); 
    
    C_local_u1u1   = zeros(NEN); 
    
    %Cu1u2
    C_local_u1u2 = zeros(NEN);
    
    %Cu2u1
    C_local_u2u1 = zeros(NEN); 
    
    %Cu2u2
    C_local_u2u2 = zeros(NEN);        
    %--------------------------------------------------------------------------
    % K matrix
    %--------------------------------------------------------------------------
    % Ku1u1
    K_local_u1u1 = zeros(NEN);
    
    % Ku1u1
    K_local_u1u2 = zeros(NEN);
   
    % Ku2u2
    K_local_u2u1 = zeros(NEN);
    
    % Ku2u2
    K_local_u2u2 = zeros(NEN);
    %--------------------------------------------------------------
    % C matrix
    %--------------------------------------------------------------
    % defining absorbing lysmer boudary at the bottom of the soil column
    % this will be used to compute the force coming from the underlying
    % elastic bedrcok
    
    for i = 1:NEN
        if NodeXYCoordinate(ElementConnectivity(e,i+1),3) == L_bottom
            jC = jC+1;
            indexC(jC) = NodeXYCoordinate(ElementConnectivity(e,i+1),1);
            CKL_local_u1u1(i,i) = CKL_local_u1u1(i,i) + damper_c_ER*El_size/2; % mid points
        end
    end
    %--------------------------------------------------------------  
    % start integration for mass and stiffness matrices
    for i_GQ = 1:QuadOrder
        
        xi_GQ = GaussQuad_coef(i_GQ); wt_xi = GaussQuad_wgt(i_GQ);
        
        for j_GQ = 1:QuadOrder
            
            eta_GQ = GaussQuad_coef(j_GQ); wt_eta = GaussQuad_wgt(j_GQ);
            
            [SF, DSFDxieta] = fun_ShapeFunction(xi_GQ, eta_GQ);
            
            Jacobian = DSFDxieta'*ElementXY; J = det(Jacobian);
            
            DSFDX = Jacobian\(DSFDxieta'); DSFDX = DSFDX';
            
            x1 = SF'*ElementXY(:,1);
            x2 = SF'*ElementXY(:,2);
            
            % the density of the soil
            El_rho = rho_ER;
            % the poisson's ratio of the soil
            El_nu = nu;
            % the shear wave velocity of the soil
            if x2 > -10
               vs = 200;
            else
               vs = 400;
            end            


            El_mu = El_rho*vs^2;
            El_E = El_mu*2*(1+El_nu);
            El_lambda = El_E *El_nu/(1+El_nu)/(1-2*El_nu);            
           
            %--------------------------------------------------------------
            % M matrix
            %--------------------------------------------------------------
            %Mu1u1
            M_local_u1u1 = M_local_u1u1...
                + El_rho*(SF)*SF'*J*wt_xi*wt_eta;
            
            %Mu2u2
            M_local_u2u2 = M_local_u2u2...
                + El_rho*(SF)*SF'*J*wt_xi*wt_eta;  
            %--------------------------------------------------------------
            % K matrix
            %--------------------------------------------------------------
            % Ku1u1
            K_local_u1u1 = K_local_u1u1...
                + ((El_lambda + 2*El_mu)*DSFDX(:,1)*DSFDX(:,1)'...
                + El_mu*DSFDX(:,2)*DSFDX(:,2)')*J*wt_xi*wt_eta;
                
            % Ku1u2
            K_local_u1u2 = K_local_u1u2...
                + (El_lambda*DSFDX(:,1)*DSFDX(:,2)'...
                + El_mu*DSFDX(:,2)*DSFDX(:,1)')*J*wt_xi*wt_eta;
            
            % Ku2u1
            K_local_u2u1 = K_local_u2u1...
                + (El_lambda*DSFDX(:,2)*DSFDX(:,1)'...
                + El_mu*DSFDX(:,1)*DSFDX(:,2)')*J*wt_xi*wt_eta;
            
            % Ku2u2
            K_local_u2u2 = K_local_u2u2...
                + (El_mu*DSFDX(:,1)*DSFDX(:,1)'...
                + (El_lambda + 2*El_mu)*DSFDX(:,2)*DSFDX(:,2)')*J*wt_xi*wt_eta;
            
        end
    end
    %--------------------------------------------------------------
    %C matrix (Rayleigh damping)
    %--------------------------------------------------------------
    C_local_u1u1 = aC(1)*M_local_u1u1 + aC(2)*K_local_u1u1;
    C_local_u1u2 = aC(2)*K_local_u1u2;
    C_local_u2u1 = aC(2)*K_local_u2u1;
    C_local_u2u2 = aC(1)*M_local_u2u2 + aC(2)*K_local_u2u2;
    % M matrix
    %----------------------------------------------------------------------
    %Mu1u1
    M(ElementDof2GlDofu1,ElementDof2GlDofu1) = ...
        M(ElementDof2GlDofu1,ElementDof2GlDofu1) + M_local_u1u1;
    
    %Mu2u2
    M(ElementDof2GlDofu2,ElementDof2GlDofu2) = ...
        M(ElementDof2GlDofu2,ElementDof2GlDofu2) + M_local_u2u2;
    %----------------------------------------------------------------------
    % C matrix
    %----------------------------------------------------------------------
    % Cu1u1
    C(ElementDof2GlDofu1,ElementDof2GlDofu1) = ...
        C(ElementDof2GlDofu1,ElementDof2GlDofu1) + C_local_u1u1;
    
    CKL(ElementDof2GlDofu1,ElementDof2GlDofu1) = ...
        CKL(ElementDof2GlDofu1,ElementDof2GlDofu1) + CKL_local_u1u1;
    
    % Cu1u2
    C(ElementDof2GlDofu1,ElementDof2GlDofu2) = ...
        C(ElementDof2GlDofu1,ElementDof2GlDofu2) + C_local_u1u2;
    
    % Cu2u1
    C(ElementDof2GlDofu2,ElementDof2GlDofu1) = ...
        C(ElementDof2GlDofu2,ElementDof2GlDofu1) + C_local_u2u1;
    
    % Cu2u2
    C(ElementDof2GlDofu2,ElementDof2GlDofu2) = ...
        C(ElementDof2GlDofu2,ElementDof2GlDofu2) + C_local_u2u2;   
    %----------------------------------------------------------------------
    % K matrix
    %----------------------------------------------------------------------
    % Ku1u1
    K(ElementDof2GlDofu1,ElementDof2GlDofu1) = ...
        K(ElementDof2GlDofu1,ElementDof2GlDofu1) + K_local_u1u1;
    
    % Ku1u2
    K(ElementDof2GlDofu1,ElementDof2GlDofu2) = ...
        K(ElementDof2GlDofu1,ElementDof2GlDofu2) + K_local_u1u2;
    
    % Ku2u1
    K(ElementDof2GlDofu2,ElementDof2GlDofu1) = ...
        K(ElementDof2GlDofu2,ElementDof2GlDofu1) + K_local_u2u1;
    
    % Ku2u2
    K(ElementDof2GlDofu2,ElementDof2GlDofu2) = ...
        K(ElementDof2GlDofu2,ElementDof2GlDofu2) + K_local_u2u2;   
end

% damping matrix due to lysmer boundary. it will be the same as C if we
% don't define any other source of damping in structure. we use this matrix
% to define the effective force at the bottom edge of the soil column.
C = C+CKL;

%disp('done with assembly')
%--------------------------------------------------------------------------
% initialization of the force vector
F = zeros(nGlDof,1); 
%--------------------------------------------------------------------------
% partitioning of matrices and vectors
% whatch out: here we only have constraints along the vertical direction!!!
KnownDof = [nnodes+EssentialBCY(:,1)];
KnownST  = [EssentialBCY(:,2)];

UnKnownDof = (1:nGlDof)'; UnKnownDof(KnownDof) = [];
%--------------------------------------------------------------------------
%disp('start partitioning')

Kuu = K(UnKnownDof,UnKnownDof);
Kun = K(UnKnownDof,KnownDof);
Knu = K(KnownDof,UnKnownDof);
Knn = K(KnownDof,KnownDof);

Cuu = C(UnKnownDof,UnKnownDof);
Cun = C(UnKnownDof,KnownDof);
Cnu = C(KnownDof,UnKnownDof);
Cnn = C(KnownDof,KnownDof);

Muu = M(UnKnownDof,UnKnownDof);
Mun = M(UnKnownDof,KnownDof);
Mnu = M(KnownDof,UnKnownDof);
Mnn = M(KnownDof,KnownDof);

%disp('end partitioning')    
%--------------------------------------------------------------------------
% Newmark time integration. 
%--------------------------------------------------------------------------
K_eff = Muu + 0.5 * Cuu * dt + 0.25 * Kuu * dt^2; 

%initialization
u_curt     = zeros(length(UnKnownDof), 1);
udot_curt  = zeros(length(UnKnownDof), 1);
uddot_curt = zeros(length(UnKnownDof), 1);

u_prev     = zeros(length(UnKnownDof), 1);
udot_prev  = zeros(length(UnKnownDof), 1);
uddot_prev = zeros(length(UnKnownDof), 1);

u_history     = zeros(nGlDof,n);
udot_history  = zeros(nGlDof,n);
uddot_history = zeros(nGlDof,n);

%--------------------------------------------------------------------------
%disp('time integration started')
FCArray = zeros(nGlDof,1);
FCArray(indexC,1) = 1;

for i = 1:length(t)
    
    cur_t = t(i) ; %disp(cur_t)
    
    % define the force vector at each time
    F = 2*CKL*Input_ugdot(i,1)*FCArray;
    
    Fn = F(UnKnownDof,1);
  
    %copying u_curt & udot_curt & uddot_curt
    % to     u_prev & udot_prev & uddot_prev
    
    u_prev = u_curt; 
    udot_prev = udot_curt; 
    uddot_prev = uddot_curt; 
    
    if (i == 1) %solving the initial acceleration
        RHS = Fn - Cuu * udot_curt - Kuu * u_curt; 
        uddot_curt = Muu\RHS ; 
        
    else %time integration.  
        RHS = Fn - Cuu * (udot_curt + uddot_curt * 0.5 * dt) ...
            - Kuu * (u_curt + udot_curt*dt+ uddot_curt * 0.25 * dt^2);
        
        uddot_curt = K_eff\RHS ; 
    end 
    
    %updating the u_curt & udot_curt
    u_curt = u_prev + udot_prev * dt + ...
            0.5 * (uddot_prev * 0.5 + uddot_curt* 0.5)*dt^2;  
    
    udot_curt = udot_prev + ...
            (uddot_prev * 0.5 + uddot_curt* 0.5)*dt; 
    
    u_history(UnKnownDof,i)     = u_curt; 
    udot_history(UnKnownDof,i)  = udot_curt; 
    uddot_history(UnKnownDof,i) = uddot_curt;
    
end 


end

