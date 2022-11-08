clc; clear; close all;
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/05/2019
%
% This code generates required DRM information for doing scattered type
% analysis.
% The output is equivalent nodal forces along the DRM interface
%
% before running this code:
% 1) update shear wave velocity and El_nusson's ratio of the homogeneous soil
%    layer in line 29
% 2) update angle of the inclined SV wave in line 39
% 3) update the coordinates of the first node that will be excited with 
%    format of [x, y] in line 44. In most cases, this node belongs to 
%    one of two nodes on the bottom line of the DRM interface
% 4) update the elevation of the ground surface in line 48. Here we assume 
%    y-axis represents the vertical direction.
% 5) update rayleigh damping coefficients in line 58-62
% 6) update the incident wave function in line 314-329
% 7) update the name of the instance that DRM nodes belong to in line 402
%
%--------------------------------------------------------------------------
% this part should be given by the user
%--------------------------------------------------------------------------
% shear wave velocity, Poisson's ratio and densityof the homogeneous soil 
% layer
Vs = 200; El_nu = 0.3; El_rho = 2000;

% Angle of the inclined SV wave
% angle(1) is the angle between the incident wave and z-axis. So if the
% wave is vertically propogated, angle(1) = 0. Please note this angle
% should not exceed the critical angle. See documentation for the values of
% the critical angle.
% angle(2) is the angle between the projected line of the incident wave on
% XoY plane and x-axis. So if the wave oscillates in x-axis, angle(2) = 0,
% and if the wave oscillates in y-axis, angle(2) = pi/2
angle = 0;

% Node_0: the coordinates of the first node that will be excited with 
%         format of [x, y, z]. In most cases, this node belongs to one of 
%         four corner nodes on the bottom surface of the DRM interface
Node_0 = [-25 -25];

% Surface_z: the elevation of the ground surface. Here we assume z-axis 
%            represents the vertical direction.
Surface_z = 0;

%--------------------------------------------------------------------------
% computing some material proterties for the homogeneous soil
El_mu = El_rho*Vs^2;
El_E = El_mu*2*(1+El_nu);
El_lambda = El_E *El_nu/(1+El_nu)/(1-2*El_nu);

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
NEN = 4;

%--------------------------------------------------------------------------
% reading mesh info from the input file
%--------------------------------------------------------------------------
fid = fopen('Job-DRM-verification-test1-input.inp');
cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true ); 
cac = cac{1};

% reading node coordinates info from the input file    
NodeXYCoordinate = str2num(char(cac(18:1908)));

% reading element connectivity info from the input file    
ElementConnectivity = str2num(char(cac(1910:3709)));

% reading inside DRM interface nodes' number info from the input file    
NodeData_DRM_inner = str2num(char(cac(3938:3943)));
NodeData_DRM_inner = reshape(NodeData_DRM_inner',[size(NodeData_DRM_inner,1)*size(NodeData_DRM_inner,2),1]);
NodeData_DRM_inner = [NodeData_DRM_inner; str2num(char(cac(3944)))'];

% reading outside DRM interface nodes' number info from the input file    
NodeData_DRM_outer = str2num(char(cac(3930:3935)));
NodeData_DRM_outer = reshape(NodeData_DRM_outer',[size(NodeData_DRM_outer,1)*size(NodeData_DRM_outer,2),1]);
NodeData_DRM_outer = [NodeData_DRM_outer; str2num(char(cac(3936)))'];

% reading DRM interface elements' number info from the input file    
% ElementData_DRM = str2num(char(cac(40330:40351)));
% ElementData_DRM = reshape(ElementData_DRM',[size(ElementData_DRM,1)*size(ElementData_DRM,2),1]);
% ElementData_DRM = [ElementData_DRM; str2num(char(cac(40352)))'];

ElementData_DRM = str2num(char(cac(3928)));
ElementData_DRM = ElementData_DRM(1):ElementData_DRM(3):ElementData_DRM(2);
    
% total number of nodes along the DRM interface
nnodesDRM = length((NodeData_DRM_inner))+length((NodeData_DRM_outer));

% total number of DOFs for the DRM nodes, i.e., DOFs = 3 for each node 
nGlDRMDof = 2*nnodesDRM;

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
    
    ElementXY = DRM_Node_Coords(DRM_Node_Index,2:3);
    
    ElementDof3GlDofu1  =                 DRM_Node_Index;
    ElementDof3GlDofu2  =     nnodesDRM + DRM_Node_Index;
        
    %Mu1u1
    M_local_u1u1 = zeros(NEN);
    
    %Mu2u2
    M_local_u2u2 = zeros(NEN);    
       
    %--------------------------------------------------------------------------
    % C matrix
    %--------------------------------------------------------------------------
    %Cu1u1
    C_local_u1u1 = zeros(NEN); 
    
    %Cu1u2
    C_local_u1u2 = zeros(NEN);
    
    %Cu2u1
    C_local_u2u1 = zeros(NEN); 
    
    %Cu2u2
    C_local_u2u2 = zeros(NEN);  
    
    %--------------------------------------------------------------------------
    % K matrix
    %--------------------------------------------------------------------------
    %Ku1u1
    K_local_u1u1 = zeros(NEN); 
    
    %Ku1u2
    K_local_u1u2 = zeros(NEN);
    
    %Ku2u1
    K_local_u2u1 = zeros(NEN); 
    
    %Ku2u2
    K_local_u2u2 = zeros(NEN);  

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
            
            El_mu = El_rho*Vs^2;
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
    % C matrix (Rayleigh damping)
    %--------------------------------------------------------------
    C_local_u1u1 = aC(1)*M_local_u1u1 + aC(2)*K_local_u1u1;
    C_local_u1u2 = aC(2)*K_local_u1u2;
    C_local_u2u1 = aC(2)*K_local_u2u1;
    C_local_u2u2 = aC(1)*M_local_u2u2 + aC(2)*K_local_u2u2;
    % M matrix
    %----------------------------------------------------------------------
    %Mu1u1
    M(ElementDof3GlDofu1,ElementDof3GlDofu1) = ...
        M(ElementDof3GlDofu1,ElementDof3GlDofu1) + M_local_u1u1;
    
    %Mu2u2
    M(ElementDof3GlDofu2,ElementDof3GlDofu2) = ...
        M(ElementDof3GlDofu2,ElementDof3GlDofu2) + M_local_u2u2;
    
    %----------------------------------------------------------------------
    % C matrix
    %----------------------------------------------------------------------
    % Cu1u1
    C(ElementDof3GlDofu1,ElementDof3GlDofu1) = ...
        C(ElementDof3GlDofu1,ElementDof3GlDofu1) + C_local_u1u1;
    
    % Cu1u2
    C(ElementDof3GlDofu1,ElementDof3GlDofu2) = ...
        C(ElementDof3GlDofu1,ElementDof3GlDofu2) + C_local_u1u2;
    
    % Cu2u1
    C(ElementDof3GlDofu2,ElementDof3GlDofu1) = ...
        C(ElementDof3GlDofu2,ElementDof3GlDofu1) + C_local_u2u1;
    
    % Cu2u2
    C(ElementDof3GlDofu2,ElementDof3GlDofu2) = ...
        C(ElementDof3GlDofu2,ElementDof3GlDofu2) + C_local_u2u2;   
    
    %----------------------------------------------------------------------
    % K matrix
    %----------------------------------------------------------------------
    % Ku1u1
    K(ElementDof3GlDofu1,ElementDof3GlDofu1) = ...
        K(ElementDof3GlDofu1,ElementDof3GlDofu1) + K_local_u1u1;
    
    % Ku1u2
    K(ElementDof3GlDofu1,ElementDof3GlDofu2) = ...
        K(ElementDof3GlDofu1,ElementDof3GlDofu2) + K_local_u1u2;
    
    % Ku2u1
    K(ElementDof3GlDofu2,ElementDof3GlDofu1) = ...
        K(ElementDof3GlDofu2,ElementDof3GlDofu1) + K_local_u2u1;
    
    % Ku2u2
    K(ElementDof3GlDofu2,ElementDof3GlDofu2) = ...
        K(ElementDof3GlDofu2,ElementDof3GlDofu2) + K_local_u2u2;   
     
end

disp('done with assembly')

% %-------------------------------------------------------------------------- 
% finding and sorting the inside and outside DRM nodes  
DRM_Inner_Node_Seq = find(ismember(DRM_Node_Coords(:,1), NodeData_DRM_inner));
DRM_Outer_Node_Seq = find(ismember(DRM_Node_Coords(:,1), NodeData_DRM_outer));

% creating the DOFs for inside and outside DRM nodes 
InnerDof = [DRM_Inner_Node_Seq; DRM_Inner_Node_Seq+nnodesDRM];
OuterDof = [DRM_Outer_Node_Seq; DRM_Outer_Node_Seq+nnodesDRM];

%--------------------------------------------------------------------------
% initialization of the force vector
Mbe = M(InnerDof,OuterDof);
Meb = M(OuterDof,InnerDof);
Cbe = C(InnerDof,OuterDof);
Ceb = C(OuterDof,InnerDof);
Kbe = K(InnerDof,OuterDof);
Keb = K(OuterDof,InnerDof);

%-------------------------------------------------------------------------- 
% defining the function of the incident wave

% defining the global variables
global dt Disp Velo Acce

% central frequency of the Ricker Pulse
fricker = 5; T = 1/fricker;

% amplitude of the Ricker Pulse 
Aricker = 1e-4;

% time that the Ricker Pulse reaches the peak value
t0 = 0.3;

% some time info for the input motion
t1 = 1.0;dt = 0.001;t = (0:dt:t1)';n = length(t);

% loading the input motion for Ricker pulse
[Acce,Velo,Disp] = ricker_new(fricker,Aricker,t,t0,dt);

%Disp = [Aricker*sin(2*pi*fricker*[0:dt:T/2]) zeros(1,n)];

% Disp = normpdf(t,t0,0.02);
% Disp = Aricker/max(Disp)*Disp;
% Velo = [0; diff(Disp)/dt];
% Acce = [0; diff(Velo)/dt];
%Velo = [2*pi*fricker*Aricker*cos(2*pi*fricker*[0:dt:T/2]) zeros(1,n)];
%Acce = [-(2*pi*fricker)^2*Aricker*sin(2*pi*fricker*[0:dt:T/2]) zeros(1,n)];

%--------------------------------------------------------------------------
% initialization of acceleration, velocity and displacement response    
Ueddot = zeros(nGlDRMDof,n);
Uedot = zeros(nGlDRMDof,n);
Ue = zeros(nGlDRMDof,n);

% computing the response for all DRM nodes
for j = 1:nnodesDRM 
        
    coords = DRM_Node_Coords(j,2:end);
%     [U,V] = fun_GetResponse_layer_2D(t0,t,coords,angle,Node_0,Surface_z,Layer_info,Velo_info_X,Velo_info_Y);
    
    [U,V] = fun_GetResponse_2D(t,coords,angle,Vs,Poi,Node_0,Surface_z);
        
    Ueddot(j,:) = U(3,:);
    Uedot(j,:) = U(2,:);
    Ue(j,:) = U(1,:);
    Ueddot(j+nnodesDRM,:) = V(3,:);
    Uedot(j+nnodesDRM,:) = V(2,:);
    Ue(j+nnodesDRM,:) = V(1,:);
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
mkdir('Files need to be imported into ABAQUS');

scale_factor = 1;

OutputName = sprintf('%s/EqNodalForceInfo.txt',FileFolderName);
fid = fopen(OutputName,'wt');
for i = 1:nnodesDRM
        fprintf(fid,'*Amplitude, name=AMP-Node%d_x, TIME=TOTAL TIME\n',i);
        fprintf(fid,'%.2f, %.1f\n',[[t(1:scale_factor:end)]; [Feff(i,1:scale_factor:end)]]);
        fprintf(fid,'*Amplitude, name=AMP-Node%d_y, TIME=TOTAL TIME\n',i);
        fprintf(fid,'%.2f, %.1f\n',[[t(1:scale_factor:end)]; [Feff(i+nnodesDRM,1:scale_factor:end)]]);
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
end
fclose(fid);

% writing the files that ask ABAQUS to including the files for equivalent
% nodal forces
OutputName = sprintf('%s/IncludeInfo.txt',FileFolderName);
fid = fopen(OutputName,'wt');
fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo.txt\n');
fprintf(fid,'*INCLUDE, INPUT=ForcesInfo.txt\n');
fclose(fid); 







