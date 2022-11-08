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
% 1) update shear wave velocity and El_nusson's ratio of the homogeneous soil
%    layer in line 29
% 2) update angle of the inclined SV wave in line 39
% 3) update the coordinates of the first node that will be excited with 
%    format of [x, y, z] in line 44. In most cases, this node belongs to 
%    one of four nodes on the bottom surface of the DRM interface
% 4) update the elevation of the ground surface in line 48. Here we assume 
%    z-axis represents the vertical direction.
% 5) update element size in line 51
% 6) update rayleigh damping coefficients in line 61-65
% 7) update the incident wave function in line 429-447
% 8) update material properties within gaussian integration
% 9) update Layer info
% 10) update angle, surface, Node_0
% 11) Run GenerateAnalyticalSolu_SV_2D_Layer.m first to get analytical
% solution
% 12) update the results name (.mat file)
%
%--------------------------------------------------------------------------
% this part should be given by the user
%--------------------------------------------------------------------------
% shear wave velocity, Poisson's ratio and densityof the homogeneous soil 
% layer
Vs = 2500; El_nu = 0.3; El_rho = 2230;

% element size
El_size = 20;

%--------------------------------------------------------------------------
% computing some material proterties for the homogeneous soil
El_mu = El_rho*Vs^2;
El_E = El_mu*2*(1+El_nu);
El_lambda = El_E *El_nu/(1+El_nu)/(1-2*El_nu);

%--------------------------------------------------------------------------
% defining rayleigh damping coefficients
xi1 = 0.2; xi2 = 0.2;
f1  = 1 ; f2  = 10;
AC  =1/4/pi*[1/f1 4*pi^2*f1;1/f2 4*pi^2*f2];
aC = AC\[xi1;xi2];
aC = [0;0.00106]; %undamped system

%--------------------------------------------------------------------------
% Gauss Quadrature order
QuadOrder = 3;
% number of nodes per elemens
NEN = 4;

%--------------------------------------------------------------------------
% reading mesh info from the input file
%--------------------------------------------------------------------------
fid = fopen('Job-2D-Layer-valley-input.inp');
cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true ); 
cac = cac{1};

% reading node coordinates info from the input file    
NodeXYCoordinate = str2num(char(cac(10:15467)));
NodeXYCoordinate(:,3) = NodeXYCoordinate(:,3);

% reading element connectivity info from the input file    
ElementConnectivity = str2num(char(cac(15469:30636)));

% reading inside DRM interface nodes' number info from the input file    
% NodeData_DRM_inner = str2num(char(cac(46511:46534)));
% NodeData_DRM_inner = reshape(NodeData_DRM_inner',[size(NodeData_DRM_inner,1)*size(NodeData_DRM_inner,2),1]);
% NodeData_DRM_inner = [NodeData_DRM_inner; str2num(char(cac(46535)))'];

NodeData_DRM_inner = str2num(char(cac(32375)));
NodeData_DRM_inner = (NodeData_DRM_inner(1):NodeData_DRM_inner(2))';

% reading outside DRM interface nodes' number info from the input file    
NodeData_DRM_outer = str2num(char(cac(32377:32396)));
NodeData_DRM_outer = reshape(NodeData_DRM_outer',[size(NodeData_DRM_outer,1)*size(NodeData_DRM_outer,2),1]);
NodeData_DRM_outer = [NodeData_DRM_outer; str2num(char(cac(32397)))'];

% reading DRM interface elements' number info from the input file    
% ElementData_DRM = str2num(char(cac(46485:46508)));
% ElementData_DRM = reshape(ElementData_DRM',[size(ElementData_DRM,1)*size(ElementData_DRM,2),1]);
% ElementData_DRM = [ElementData_DRM; str2num(char(cac(46509)))'];

ElementData_DRM = str2num(char(cac(32373)));
ElementData_DRM = (ElementData_DRM(1):ElementData_DRM(2))';


% ElementData_DRM = str2num(char(cac(3928)));
% ElementData_DRM = ElementData_DRM(1):ElementData_DRM(3):ElementData_DRM(2);
    
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
            
            if x2 >= -1000
                Vs = 1500;
                El_rho = 2050;
                El_nu = 0.35;
            else
                Vs = 2500;
                El_rho = 2230;
                El_nu = 0.3;
            end
            
            
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

% defining the function of the incident wave

Layer_info = [1000 2050 1500 0.35;
               400 2230 2500  0.3];

           
angle = 0;
Node_0 = [-1900 -1400];
Surface_z = 0;

dt = 0.01;
t0 = 0:dt:100;
t = 0:dt:15;
n = length(t);

[Acce,Velo,Disp] = ricker_new(1,1e-1,t0,1,dt);

[ut,wt] = fun_GenerateAnalyticalSolu_SV_2D_Layer(t0,dt,Velo,Layer_info,angle);

Velo_info_X = ut;
Velo_info_Y = wt;


%--------------------------------------------------------------------------
% initialization of acceleration, velocity and displacement response    
Ueddot = zeros(nGlDRMDof,n);
Uedot = zeros(nGlDRMDof,n);
Ue = zeros(nGlDRMDof,n);

% computing the response for all DRM nodes

for j = 1:nnodesDRM 
        
    coords = DRM_Node_Coords(j,2:end);
    [U,V] = fun_GetResponse_layer_2D(t0,t,coords,angle,Node_0,Surface_z,Layer_info,Velo_info_X,Velo_info_Y);
    %[U,V] = fun_GetResponse(time,coords,0,25);
        
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

%--------------------------------------------------------------------------
% writing the files that include the equivalent nodal forces for all DRM
% nodes in 3 directions (Fx, Fy, Fz)
mkdir('Files need to be imported into ABAQUS 2D Layer');

for i = 1:nnodesDRM
    OutputName = sprintf('Files need to be imported into ABAQUS 2D Layer/EqNodalForceInfo%d_x.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*Amplitude, name=AMP-Node%d_x\n',i);
    fprintf(fid,'%.4f, %.1f\n',[[t(1:1:end)]; [Feff(i,1:1:end)]]);
    fclose(fid);
    
    OutputName = sprintf('Files need to be imported into ABAQUS 2D Layer/EqNodalForceInfo%d_y.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*Amplitude, name=AMP-Node%d_y\n',i);
    fprintf(fid,'%.4f, %.1f\n',[[t(1:1:end)]; [Feff(i+nnodesDRM,1:1:end)]]);
    fclose(fid);      
end

% writing the files that teach ABAQUS to understand which file belongs to
% which node
fid = fopen('Files need to be imported into ABAQUS 2D Layer/ForcesInfo.txt','wt');
for i = 1:nnodesDRM
    fprintf(fid,'*Cload, amplitude=AMP-Node%d_x\n',i);
    fprintf(fid,'Part-1-1.%d, 1, 1.\n',DRM_Node_Data(i));
    
    fprintf(fid,'*Cload, amplitude=AMP-Node%d_y\n',i);
    fprintf(fid,'Part-1-1.%d, 2, 1.\n',DRM_Node_Data(i)); 
end
fclose(fid);

% writing the files that ask ABAQUS to including the files for equivalent
% nodal forces
fid = fopen('Files need to be imported into ABAQUS 2D Layer/IncludeInfo.txt','wt');
for i = 1:nnodesDRM 
    fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo%d_x.txt\n',i);
    fprintf(fid,'*INCLUDE, INPUT=EqNodalForceInfo%d_y.txt\n',i);

end
fprintf(fid,'*INCLUDE, INPUT=ForcesInfo.txt\n');
fclose(fid);








