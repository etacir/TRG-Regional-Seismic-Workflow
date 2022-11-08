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
% 1) update shear wave velocity and Poisson's ratio of the homogeneous soil
%    layer in line 29
% 2) update angle of the inclined SV wave in line 39
% 3) update the coordinates of the first node that will be excited with 
%    format of [x, y, z] in line 44. In most cases, this node belongs to 
%    one of four nodes on the bottom surface of the DRM interface
% 4) update the elevation of the ground surface in line 48. Here we assume 
%    z-axis represents the vertical direction.
% 5) update rayleigh damping coefficients in line 58-62
% 6) update the incident wave function in line 431-446
% 7) update the name of the instance that DRM nodes belong to in line 516
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
angle = [0 0];

% Node_0: the coordinates of the first node that will be excited with 
%         format of [x, y, z]. In most cases, this node belongs to one of 
%         four corner nodes on the bottom surface of the DRM interface
Node_0 = [-25,-25,-25];

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
NEN = 8;

%--------------------------------------------------------------------------
% reading mesh info from the input file
%--------------------------------------------------------------------------
fid = fopen('Job-DRM-verification-3D-test-input.inp');
cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true ); 
cac = cac{1};

% reading node coordinates info from the input file    
NodeXYCoordinate = str2num(char(cac(19:115369)));
NodeXYCoordinate(:,4) = NodeXYCoordinate(:,4)-30;


% reading element connectivity info from the input file    
ElementConnectivity = str2num(char(cac(115371:223370)));

% reading inside DRM interface nodes' number info from the input file    
NodeData_DRM_inner = str2num(char(cac(238749:239186)));
NodeData_DRM_inner = reshape(NodeData_DRM_inner',[size(NodeData_DRM_inner,1)*size(NodeData_DRM_inner,2),1]);
NodeData_DRM_inner = [NodeData_DRM_inner; str2num(char(cac(239187)))'];

% reading outside DRM interface nodes' number info from the input file    
NodeData_DRM_outer = str2num(char(cac(238272:238746)));
NodeData_DRM_outer = reshape(NodeData_DRM_outer',[size(NodeData_DRM_outer,1)*size(NodeData_DRM_outer,2),1]);
NodeData_DRM_outer = [NodeData_DRM_outer; str2num(char(cac(238747)))'];

% reading DRM interface elements' number info from the input file    
ElementData_DRM = str2num(char(cac(237820:238269)));
ElementData_DRM = reshape(ElementData_DRM',[size(ElementData_DRM,1)*size(ElementData_DRM,2),1]);
ElementData_DRM = [ElementData_DRM; str2num(char(cac(238270)))'];

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

%--------------------------------------------------------------------------
% initialization of acceleration, velocity and displacement response    
Ueddot = zeros(nGlDRMDof,n);
Uedot = zeros(nGlDRMDof,n);
Ue = zeros(nGlDRMDof,n);

% computing the response for all DRM nodes
for i = 1:n-1
    disp(i)
    for j = 1:nnodesDRM 
        
        coords = DRM_Node_Coords(j,2:4);
        time = t(i);
        [U,V,W] = fun_GetResponse_3D(time,coords,angle,Vs,El_nu,Node_0,Surface_z);
        
        Ueddot(j,i) = U(3);
        Uedot(j,i) = U(2);
        Ue(j,i) = U(1);
        Ueddot(j+nnodesDRM,i) = V(3);
        Uedot(j+nnodesDRM,i) = V(2);
        Ue(j+nnodesDRM,i) = V(1);
        Ueddot(j+nnodesDRM*2,i) = W(3);
        Uedot(j+nnodesDRM*2,i) = W(2);
        Ue(j+nnodesDRM*2,i) = W(1);
    end       
end

% computing the equivalent nodal forces for all DRM nodes
Feff = zeros(nGlDRMDof,n);
for i = 1:n
    Feff(InnerDof,i) = -Mbe*Ueddot(OuterDof,i) - Cbe*Uedot(OuterDof,i) - Kbe*Ue(OuterDof,i);
    Feff(OuterDof,i) = Meb*Ueddot(InnerDof,i) + Ceb*Uedot(InnerDof,i) + Keb*Ue(InnerDof,i);
end 

%--------------------------------------------------------------------------
% writing the files that include the equivalent nodal forces for all DRM
% nodes in 3 directions (Fx, Fy, Fz)
mkdir('Files need to be imported into ABAQUS 3D');

for i = 1:nnodesDRM
    OutputName = sprintf('Files need to be imported into ABAQUS 3D/EqNodalForceInfo%d_x.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*Amplitude, name=AMP-Node%d_x\n',i);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i,j));
    end
    fclose(fid);
    
    OutputName = sprintf('Files need to be imported into ABAQUS 3D/EqNodalForceInfo%d_y.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*Amplitude, name=AMP-Node%d_y\n',i);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i+nnodesDRM,j));
    end
    fclose(fid);   
    
    OutputName = sprintf('Files need to be imported into ABAQUS 3D/EqNodalForceInfo%d_z.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*Amplitude, name=AMP-Node%d_z\n',i);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i+nnodesDRM*2,j));
    end
    fclose(fid);
end

% writing the files that creates node sets for each DRM node
fid = fopen('Files need to be imported into ABAQUS 3D/NodeSetsInfo.txt','wt');
for i = 1:nnodesDRM
    fprintf(fid,'*Nset, nset=NodeSet%d, internal, instance=PART-1-1\n',DRM_Node_Data(i));
    fprintf(fid,' %d,\n',DRM_Node_Data(i));
end
fclose(fid);

% writing the files that teach ABAQUS to understand which file belongs to
% which node
fid = fopen('Files need to be imported into ABAQUS 3D/ForcesInfo.txt','wt');
for i = 1:nnodesDRM
    fprintf(fid,'*Cload, amplitude=AMP-Node%d_x\n',i);
    fprintf(fid,'NodeSet%d, 1, 1.\n',DRM_Node_Data(i));
    
    fprintf(fid,'*Cload, amplitude=AMP-Node%d_y\n',i);
    fprintf(fid,'NodeSet%d, 2, 1.\n',DRM_Node_Data(i)); 

    fprintf(fid,'*Cload, amplitude=AMP-Node%d_z\n',i);
    fprintf(fid,'NodeSet%d, 3, 1.\n',DRM_Node_Data(i));
end
fclose(fid);

% writing the files that ask ABAQUS to including the files for equivalent
% nodal forces
fid = fopen('Files need to be imported into ABAQUS 3D/IncludeInfo.txt','wt');
for i = 1:nnodesDRM 
    fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo%d_x.txt\n',i);
    fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo%d_y.txt\n',i);
    fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo%d_z.txt\n',i);

end
fprintf(fid,'*INCLUDE, INPUT=NodeSetsInfo.txt\n');
fprintf(fid,'*INCLUDE, INPUT=ForcesInfo.txt\n');
fclose(fid);




