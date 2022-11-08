clc; clear; close all;
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
% 1) update element size in line 35
% 2) update rayleigh damping coefficients in line 39-43
% 3) update density of the soil under certain elevation in line 210, this
%    can be defined as a piecewise function by using the coordinates of
%    that integration points (x1, x2, x3)
% 4) update poisson's ratio of the soil under certain elevation in line 
%    212, this can be defined as a piecewise function by using the 
%    coordinates of that integration points (x1, x2, x3)
% 5) update shear wave velocity of the soil under certain elevation in line 
%    214, this can be defined as a piecewise function by using the 
%    coordinates of that integration points (x1, x2, x3)
% 6) update the free-field motion computed from the 1D site response
%    analysis in line 417
% 7) update the mesh info for the soil column that is used for computing 1D
%    site response analysis in line 420
% 8) update the timestep info for the incident wave you applied in line 423
% 9) update material properties within gaussian integration
% 10) update Layer info
% 11) update angle, surface, Node_0
% 12) Run GenerateAnalyticalSolu_SV_3D_Layer.m first to get analytical
% solution
% 13) update the results name (.mat file)
%
%--------------------------------------------------------------------------
% this part should be given by the user
%--------------------------------------------------------------------------
% damshpot coefficient for the soil below the DRM interface
Vs_ER = 400; rho_ER = 2000; damper_c_ER = rho_ER * Vs_ER;
V0 = 200; VH = Vs_ER; nvs = 0.5; b = (V0/VH)^(1/nvs);
L_bottom = -25; 

% element size
El_size = 20;

%--------------------------------------------------------------------------
% defining rayleigh damping coefficients
xi1 = 0.1; xi2 = 0.1;
f1  = 0.35 ; f2  = 20;
AC  =1/4/pi*[1/f1 4*pi^2*f1;1/f2 4*pi^2*f2];
aC = AC\[xi1;xi2];
aC = [0;0]; %undamped system

%--------------------------------------------------------------------------
% Gauss Quadrature order
QuadOrder = 3;
% number of nodes per elemens
NEN = 8;

%--------------------------------------------------------------------------
% reading mesh info from the input file
%--------------------------------------------------------------------------
fid = fopen('Job-3D-Layer-canyon-input.inp');
cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true ); 
cac = cac{1};

% reading node coordinates info from the input file    
NodeXYCoordinate = str2num(char(cac(10:163008)));
NodeXYCoordinate(:,4) = NodeXYCoordinate(:,4)-400;


% reading element connectivity info from the input file    
ElementConnectivity = str2num(char(cac(163010:213409)));

% reading inside DRM interface nodes' number info from the input file    
% NodeData_DRM_inner = str2num(char(cac(238749:239186)));
% NodeData_DRM_inner = reshape(NodeData_DRM_inner',[size(NodeData_DRM_inner,1)*size(NodeData_DRM_inner,2),1]);
% NodeData_DRM_inner = [NodeData_DRM_inner; str2num(char(cac(239187)))'];
NodeData_DRM_inner = str2num(char(cac(300133)));
NodeData_DRM_inner = (NodeData_DRM_inner(1):NodeData_DRM_inner(2))';

% reading outside DRM interface nodes' number info from the input file    
NodeData_DRM_outer = str2num(char(cac(300135:300459)));
NodeData_DRM_outer = reshape(NodeData_DRM_outer',[size(NodeData_DRM_outer,1)*size(NodeData_DRM_outer,2),1]);
NodeData_DRM_outer = [NodeData_DRM_outer; str2num(char(cac(300460)))'];

% reading DRM interface elements' number info from the input file    
% ElementData_DRM = str2num(char(cac(237820:238269)));
% ElementData_DRM = reshape(ElementData_DRM',[size(ElementData_DRM,1)*size(ElementData_DRM,2),1]);
% ElementData_DRM = [ElementData_DRM; str2num(char(cac(238270)))'];
ElementData_DRM = str2num(char(cac(300131)));
ElementData_DRM = (ElementData_DRM(1):ElementData_DRM(2))';

% total number of nodes along the DRM interface
nnodesDRM = length((NodeData_DRM_inner))+length((NodeData_DRM_outer));

% total number of DOFs for the DRM nodes, i.e., DOFs = 3 for each node 
nGlDRMDof = 3*nnodesDRM;

% sorting the DRM nodes' numbers and create a vector for those numbers
DRM_Node_Data = sort([NodeData_DRM_inner;NodeData_DRM_outer]);

% getting the nodes' coordinates for all DRM nodes 
DRM_Node_Coords = find(ismember(NodeXYCoordinate(:,1), DRM_Node_Data));
DRM_Node_Coords = NodeXYCoordinate(DRM_Node_Coords,:);
%--------------------------------------------------------------------------
% reading gauss quadrature info
[GaussQuad_coef, GaussQuad_wgt] = fun_Order3GaussQuadDataOrder;
%--------------------------------------------------------------------------
% initialization of the mass, damping, and stiffness matrices
M = zeros(nGlDRMDof); C = zeros(nGlDRMDof); K = zeros(nGlDRMDof);

%--------------------------------------------------------------------------
% computing elemental matrices and assembling
for e = 1:length(ElementData_DRM)
    
    disp(e)
    
    DRM_Ele_Index = find(ElementConnectivity(:,1) == ElementData_DRM(e));
    ElementDof = ElementConnectivity(DRM_Ele_Index,2:end);

    [~,rank1] = sort(ElementDof);
    [~,rank2] = sort(rank1);
    DRM_Node_Index = find(ismember(DRM_Node_Coords(:,1), ElementDof));
    DRM_Node_Index = DRM_Node_Index(rank2);
    
    ElementXY = DRM_Node_Coords(DRM_Node_Index,2:4);
    
    ElementDof3GlDofu1  =                 DRM_Node_Index;
    ElementDof3GlDofu2  =     nnodesDRM + DRM_Node_Index;
    ElementDof3GlDofu3  =   nnodesDRM*2 + DRM_Node_Index;
        
    %Mu1u1
    M_local_u1u1 = zeros(NEN);
    
    %Mu2u2
    M_local_u2u2 = zeros(NEN);    
    
    %Mu3u3
    M_local_u3u3 = zeros(NEN);    
    %--------------------------------------------------------------------------
    % C matrix
    %--------------------------------------------------------------------------
    %Cu1u1
    C_local_u1u1 = zeros(NEN); 
    
    %Cu1u2
    C_local_u1u2 = zeros(NEN);
    
    %Cu1u3
    C_local_u1u3 = zeros(NEN);
    
    %Cu2u1
    C_local_u2u1 = zeros(NEN); 
    
    %Cu2u2
    C_local_u2u2 = zeros(NEN);  
    
    %Cu2u3
    C_local_u2u3 = zeros(NEN);
    
    %Cu3u1
    C_local_u3u1 = zeros(NEN); 
    
    %Cu3u2
    C_local_u3u2 = zeros(NEN);  
    
    %Cu3u3
    C_local_u3u3 = zeros(NEN);
    %--------------------------------------------------------------------------
    % K matrix
    %--------------------------------------------------------------------------
    %Ku1u1
    K_local_u1u1 = zeros(NEN); 
    
    %Ku1u2
    K_local_u1u2 = zeros(NEN);
    
    %Ku1u3
    K_local_u1u3 = zeros(NEN);
    
    %Ku2u1
    K_local_u2u1 = zeros(NEN); 
    
    %Ku2u2
    K_local_u2u2 = zeros(NEN);  
    
    %Ku2u3
    K_local_u2u3 = zeros(NEN);
    
    %Ku3u1
    K_local_u3u1 = zeros(NEN); 
    
    %Ku3u2
    K_local_u3u2 = zeros(NEN);  
    
    %Ku3u3
    K_local_u3u3 = zeros(NEN);

    %--------------------------------------------------------------  
    % start integration for mass and stiffness matrices
    for i_GQ = 1:QuadOrder
        
        xi_GQ = GaussQuad_coef(i_GQ); wt_xi = GaussQuad_wgt(i_GQ);
        
        for j_GQ = 1:QuadOrder
            
            eta_GQ = GaussQuad_coef(j_GQ); wt_eta = GaussQuad_wgt(j_GQ);
                  
            for k_GQ = 1:QuadOrder
            
                mu_GQ = GaussQuad_coef(k_GQ); wt_mu = GaussQuad_wgt(k_GQ);
            
                [SF, DSFDxieta] = fun_ShapeFunction_3D(xi_GQ, eta_GQ, mu_GQ);
            
            	Jacobian = DSFDxieta'*ElementXY; J = det(Jacobian);
            
                DSFDX = Jacobian\(DSFDxieta'); DSFDX = DSFDX';
            
                x1 = SF'*ElementXY(:,1);
                x2 = SF'*ElementXY(:,2);
                x3 = SF'*ElementXY(:,3);
                           
%                 % the density of the soil
%                 El_rho = 2000;
%                 % the poisson's ratio of the soil
%                 El_nu = 0.3;
%                 % the shear wave velocity of the soil
%                 Vs = VH*(b+(1-b)*abs(x3)/abs(L_bottom))^nvs;
                
                if x3 >= -200
                    Vs = 400;
                    El_rho = 2000;
                    El_nu = 1/3;
                    aC = [0; 2*0.05/2/pi/1.5];
                else
                    Vs = 800;
                    El_rho = 2000;
                    El_nu = 1/3;
                    aC = [0; 2*0.02/2/pi/1.5];
                end

                El_mu = El_rho*Vs^2;
                El_E = El_mu*2*(1+El_nu);
                El_lambda = El_E *El_nu/(1+El_nu)/(1-2*El_nu);

                %--------------------------------------------------------------
                % M matrix
                %--------------------------------------------------------------
                %Mu1u1
                M_local_u1u1 = M_local_u1u1...
                    + El_rho*(SF)*SF'*J*wt_xi*wt_eta*wt_mu;

                %Mu2u2
                M_local_u2u2 = M_local_u2u2...
                    + El_rho*(SF)*SF'*J*wt_xi*wt_eta*wt_mu; 
                
                %Mu3u3
                M_local_u3u3 = M_local_u3u3...
                    + El_rho*(SF)*SF'*J*wt_xi*wt_eta*wt_mu; 
                
                %--------------------------------------------------------------
                % K matrix
                %--------------------------------------------------------------
                % Ku1u1
                K_local_u1u1 = K_local_u1u1...
                    + ((El_lambda + 2*El_mu)*DSFDX(:,1)*DSFDX(:,1)'...
                    + El_mu*DSFDX(:,2)*DSFDX(:,2)'....
                    + El_mu*DSFDX(:,3)*DSFDX(:,3)')*J*wt_xi*wt_eta*wt_mu;

                % Ku1u2
                K_local_u1u2 = K_local_u1u2...
                    + (El_lambda*DSFDX(:,1)*DSFDX(:,2)'...
                    + El_mu*DSFDX(:,2)*DSFDX(:,1)')*J*wt_xi*wt_eta*wt_mu;
                
                % Ku1u3
                K_local_u1u3 = K_local_u1u3...
                    + (El_lambda*DSFDX(:,1)*DSFDX(:,3)'...
                    + El_mu*DSFDX(:,3)*DSFDX(:,1)')*J*wt_xi*wt_eta*wt_mu;

                % Ku2u1
                K_local_u2u1 = K_local_u2u1...
                    + (El_lambda*DSFDX(:,2)*DSFDX(:,1)'...
                    + El_mu*DSFDX(:,1)*DSFDX(:,2)')*J*wt_xi*wt_eta*wt_mu;

                % Ku2u2
                K_local_u2u2 = K_local_u2u2...
                    + (El_mu*DSFDX(:,1)*DSFDX(:,1)'...
                    + (El_lambda + 2*El_mu)*DSFDX(:,2)*DSFDX(:,2)'...
                    +  El_mu*DSFDX(:,3)*DSFDX(:,3)')*J*wt_xi*wt_eta*wt_mu;
                
                % Ku2u3
                K_local_u2u3 = K_local_u2u3...
                    + (El_lambda*DSFDX(:,2)*DSFDX(:,3)'...
                    + El_mu*DSFDX(:,3)*DSFDX(:,2)')*J*wt_xi*wt_eta*wt_mu;

                % Ku3u1
                K_local_u3u1 = K_local_u3u1...
                    + (El_lambda*DSFDX(:,3)*DSFDX(:,1)'...
                    + El_mu*DSFDX(:,1)*DSFDX(:,3)')*J*wt_xi*wt_eta*wt_mu;

                % Ku3u2
                K_local_u3u2 = K_local_u3u2...
                    + (El_lambda*DSFDX(:,3)*DSFDX(:,2)'...
                    + El_mu*DSFDX(:,2)*DSFDX(:,3)')*J*wt_xi*wt_eta*wt_mu;      
                
                % Ku3u3
                K_local_u3u3 = K_local_u3u3...
                    + (El_mu*DSFDX(:,1)*DSFDX(:,1)'...
                    + El_mu*DSFDX(:,2)*DSFDX(:,2)'...
                    + (El_lambda + 2*El_mu)*DSFDX(:,3)*DSFDX(:,3)')*J*wt_xi*wt_eta*wt_mu;     
            end
        end
    end
    %--------------------------------------------------------------
    % C matrix (Rayleigh damping)
    %--------------------------------------------------------------
    C_local_u1u1 = aC(1)*M_local_u1u1 + aC(2)*K_local_u1u1;
    C_local_u1u2 = aC(2)*K_local_u1u2;
    C_local_u1u3 = aC(2)*K_local_u1u3;
    C_local_u2u1 = aC(2)*K_local_u2u1;
    C_local_u2u2 = aC(1)*M_local_u2u2 + aC(2)*K_local_u2u2;
    C_local_u2u3 = aC(2)*K_local_u2u3;
    C_local_u3u1 = aC(2)*K_local_u3u1;
    C_local_u3u2 = aC(2)*K_local_u3u2;
    C_local_u3u3 = aC(1)*M_local_u3u3 + aC(2)*K_local_u3u3;
    % M matrix
    %----------------------------------------------------------------------
    %Mu1u1
    M(ElementDof3GlDofu1,ElementDof3GlDofu1) = ...
        M(ElementDof3GlDofu1,ElementDof3GlDofu1) + M_local_u1u1;
    
    %Mu2u2
    M(ElementDof3GlDofu2,ElementDof3GlDofu2) = ...
        M(ElementDof3GlDofu2,ElementDof3GlDofu2) + M_local_u2u2;
    
    %Mu3u3
    M(ElementDof3GlDofu3,ElementDof3GlDofu3) = ...
        M(ElementDof3GlDofu3,ElementDof3GlDofu3) + M_local_u3u3;
    
    %----------------------------------------------------------------------
    % C matrix
    %----------------------------------------------------------------------
    % Cu1u1
    C(ElementDof3GlDofu1,ElementDof3GlDofu1) = ...
        C(ElementDof3GlDofu1,ElementDof3GlDofu1) + C_local_u1u1;
    
    % Cu1u2
    C(ElementDof3GlDofu1,ElementDof3GlDofu2) = ...
        C(ElementDof3GlDofu1,ElementDof3GlDofu2) + C_local_u1u2;
    
    % Cu1u3
    C(ElementDof3GlDofu1,ElementDof3GlDofu3) = ...
        C(ElementDof3GlDofu1,ElementDof3GlDofu3) + C_local_u1u3;
    
    % Cu2u1
    C(ElementDof3GlDofu2,ElementDof3GlDofu1) = ...
        C(ElementDof3GlDofu2,ElementDof3GlDofu1) + C_local_u2u1;
    
    % Cu2u2
    C(ElementDof3GlDofu2,ElementDof3GlDofu2) = ...
        C(ElementDof3GlDofu2,ElementDof3GlDofu2) + C_local_u2u2;   
    
    % Cu2u3
    C(ElementDof3GlDofu2,ElementDof3GlDofu3) = ...
        C(ElementDof3GlDofu2,ElementDof3GlDofu3) + C_local_u2u3;  
    
    % Cu3u1
    C(ElementDof3GlDofu3,ElementDof3GlDofu1) = ...
        C(ElementDof3GlDofu3,ElementDof3GlDofu1) + C_local_u3u1;
    
    % Cu3u2
    C(ElementDof3GlDofu3,ElementDof3GlDofu2) = ...
        C(ElementDof3GlDofu3,ElementDof3GlDofu2) + C_local_u3u2;   
    
    % Cu3u3
    C(ElementDof3GlDofu3,ElementDof3GlDofu3) = ...
        C(ElementDof3GlDofu3,ElementDof3GlDofu3) + C_local_u3u3; 
    
    %----------------------------------------------------------------------
    % K matrix
    %----------------------------------------------------------------------
    % Ku1u1
    K(ElementDof3GlDofu1,ElementDof3GlDofu1) = ...
        K(ElementDof3GlDofu1,ElementDof3GlDofu1) + K_local_u1u1;
    
    % Ku1u2
    K(ElementDof3GlDofu1,ElementDof3GlDofu2) = ...
        K(ElementDof3GlDofu1,ElementDof3GlDofu2) + K_local_u1u2;
    
    % Ku1u3
    K(ElementDof3GlDofu1,ElementDof3GlDofu3) = ...
        K(ElementDof3GlDofu1,ElementDof3GlDofu3) + K_local_u1u3;
    
    % Ku2u1
    K(ElementDof3GlDofu2,ElementDof3GlDofu1) = ...
        K(ElementDof3GlDofu2,ElementDof3GlDofu1) + K_local_u2u1;
    
    % Ku2u2
    K(ElementDof3GlDofu2,ElementDof3GlDofu2) = ...
        K(ElementDof3GlDofu2,ElementDof3GlDofu2) + K_local_u2u2;   
    
    % Ku2u3
    K(ElementDof3GlDofu2,ElementDof3GlDofu3) = ...
        K(ElementDof3GlDofu2,ElementDof3GlDofu3) + K_local_u2u3;  
    
    % Ku3u1
    K(ElementDof3GlDofu3,ElementDof3GlDofu1) = ...
        K(ElementDof3GlDofu3,ElementDof3GlDofu1) + K_local_u3u1;
    
    % Ku3u2
    K(ElementDof3GlDofu3,ElementDof3GlDofu2) = ...
        K(ElementDof3GlDofu3,ElementDof3GlDofu2) + K_local_u3u2;   
    
    % Ku3u3
    K(ElementDof3GlDofu3,ElementDof3GlDofu3) = ...
        K(ElementDof3GlDofu3,ElementDof3GlDofu3) + K_local_u3u3;    
end

disp('done with assembly')

% %-------------------------------------------------------------------------- 
% finding and sorting the inside and outside DRM nodes  
DRM_Inner_Node_Seq = find(ismember(DRM_Node_Coords(:,1), NodeData_DRM_inner));
DRM_Outer_Node_Seq = find(ismember(DRM_Node_Coords(:,1), NodeData_DRM_outer));

% creating the DOFs for inside and outside DRM nodes 
InnerDof = [DRM_Inner_Node_Seq; DRM_Inner_Node_Seq+nnodesDRM; DRM_Inner_Node_Seq+nnodesDRM*2];
OuterDof = [DRM_Outer_Node_Seq; DRM_Outer_Node_Seq+nnodesDRM; DRM_Outer_Node_Seq+nnodesDRM*2];

%--------------------------------------------------------------------------
% initialization of the force vector
Mbe = M(InnerDof,OuterDof);
Meb = M(OuterDof,InnerDof);
Cbe = C(InnerDof,OuterDof);
Ceb = C(OuterDof,InnerDof);
Kbe = K(InnerDof,OuterDof);
Keb = K(OuterDof,InnerDof);

% %-------------------------------------------------------------------------- 
Layer_info = [200 2000 400 1/3;
              200 2000 800 1/3];

           
Angle_h = 45;
Angle_v = 29.9999999999999;
Node_0 = [-700 -300 -300];
Surface_z = 0;

dt = 0.01;
t0 = 0:dt:100;
t = 0:dt:15;
n = length(t);

[Acce,Velo,Disp] = ricker_new(1,1e-1,t0,1,dt);

[ut,wt] = fun_GenerateAnalyticalSolu_SV_2D_Layer(t0,dt,Velo,Layer_info,Angle_v);

Velo_info_X = ut;
Velo_info_Y = wt;

% initialization of acceleration, velocity and displacement response    
Ueddot = zeros(nGlDRMDof,n);
Uedot = zeros(nGlDRMDof,n);
Ue = zeros(nGlDRMDof,n);

% computing the response for all DRM nodes
for j = 1:nnodesDRM 
        disp(j)

        coords = DRM_Node_Coords(j,2:4);        
        [U,V,W] = fun_GetResponse_layer_3D(t0,t,coords,Angle_h,Angle_v,Node_0,Surface_z,Layer_info,Velo_info_X,Velo_info_Y);

        Ueddot(j,:) = U(3,:);
        Uedot(j,:) = U(2,:);
        Ue(j,:) = U(1,:);
        Ueddot(j+nnodesDRM,:) = V(3,:);
        Uedot(j+nnodesDRM,:) = V(2,:);
        Ue(j+nnodesDRM,:) = V(1,:);
        Ueddot(j+nnodesDRM*2,:) = W(3,:);
        Uedot(j+nnodesDRM*2,:) = W(2,:);
        Ue(j+nnodesDRM*2,:) = W(1,:);
end  

% computing the equivalent nodal forces for all DRM nodes
Feff = zeros(nGlDRMDof,n);
Feff(InnerDof,:) = -Mbe*Ueddot(OuterDof,:) - Cbe*Uedot(OuterDof,:) - Kbe*Ue(OuterDof,:);
Feff(OuterDof,:) = Meb*Ueddot(InnerDof,:) + Ceb*Uedot(InnerDof,:) + Keb*Ue(InnerDof,:);
%StaticLoad=interp1(Static_t,StaticLoad,t);
% Feff(InnerDof,:) = Feff(InnerDof,:) + StaticLoad'*1;
    
%--------------------------------------------------------------------------
% writing the files that include the equivalent nodal forces for all DRM
% nodes in 3 directions (Fx, Fy, Fz)
FileFolderName = sprintf('Files need to be imported into ABAQUS 3D Layer');
mkdir(FileFolderName);

scale_factor = 1;

OutputName = sprintf('%s/EqNodalForceInfo.txt',FileFolderName);
fid = fopen(OutputName,'wt');
for i = 1:nnodesDRM
        fprintf(fid,'*Amplitude, name=AMP-Node%d_x, TIME=TOTAL TIME\n',i);
        fprintf(fid,'%.2f, %.1f\n',[[t(1:scale_factor:end)]; [Feff(i,1:scale_factor:end)]]);
        fprintf(fid,'*Amplitude, name=AMP-Node%d_y, TIME=TOTAL TIME\n',i);
        fprintf(fid,'%.2f, %.1f\n',[[t(1:scale_factor:end)]; [Feff(i+nnodesDRM,1:scale_factor:end)]]);
        fprintf(fid,'*Amplitude, name=AMP-Node%d_z, TIME=TOTAL TIME\n',i);
        fprintf(fid,'%.2f, %.1f\n',[[t(1:scale_factor:end)]; [Feff(i+nnodesDRM*2,1:scale_factor:end)]]);
end
fclose(fid);

% writing the files that teach ABAQUS to understand which file belongs to
% which node
OutputName = sprintf('%s/ForcesInfo.txt',FileFolderName);
fid = fopen(OutputName,'wt');
for i = 1:nnodesDRM
        fprintf(fid,'*Cload, amplitude=AMP-Node%d_x\n',i);
        fprintf(fid,'PART-1-1.%d, 1, 1.\n',DRM_Node_Data(i));

        fprintf(fid,'*Cload, amplitude=AMP-Node%d_y\n',i);
        fprintf(fid,'PART-1-1.%d, 2, 1.\n',DRM_Node_Data(i)); 

        fprintf(fid,'*Cload, amplitude=AMP-Node%d_z\n',i);
        fprintf(fid,'PART-1-1.%d, 3, 1.\n',DRM_Node_Data(i));
end
fclose(fid);

% writing the files that ask ABAQUS to including the files for equivalent
% nodal forces
OutputName = sprintf('%s/IncludeInfo.txt',FileFolderName);
fid = fopen(OutputName,'wt');
fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo.txt\n');
fprintf(fid,'*INCLUDE, INPUT=ForcesInfo.txt\n');
fclose(fid); 



