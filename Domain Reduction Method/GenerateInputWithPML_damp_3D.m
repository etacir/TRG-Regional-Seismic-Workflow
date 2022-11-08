clc;clear;close all;

fid = fopen('Job-DRM-vefification-3D-input.inp');
cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true ); 
cac = cac{1};
fclose(fid);

ElementConnectivity = str2num(char(cac(67637:130136)));
PML_ElementConnectivity = ElementConnectivity;
PML_ElementConnectivity(:,1) = PML_ElementConnectivity(:,1) + 5e5;

OutputName = sprintf('PML_Info_3D.txt');
fid = fopen(OutputName,'wt');


PML_Element_Set = str2num(char(cac(135437)));
PML_Element_Index = PML_Element_Set(1):PML_Element_Set(3):PML_Element_Set(2);

OutputName = sprintf('Files need to be imported into ABAQUS 3D Layer/PML_Info_3D.txt');
fid = fopen(OutputName,'wt');
fprintf(fid, '*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=12, COORDINATES=3, VARIABLES=360, unsym\n');
fprintf(fid, '1,2,3,21,22,23,24,25,26\n');
fprintf(fid,'*Element, type=U1, elset=el-PML-type-1\n');
for j = 1:length(PML_Element_Index)
        fprintf(fid,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',...
                PML_ElementConnectivity(PML_Element_Index(j),1),...
                PML_ElementConnectivity(PML_Element_Index(j),2),...
                PML_ElementConnectivity(PML_Element_Index(j),3),...
                PML_ElementConnectivity(PML_Element_Index(j),4),...
                PML_ElementConnectivity(PML_Element_Index(j),5),...
                PML_ElementConnectivity(PML_Element_Index(j),6),...
                PML_ElementConnectivity(PML_Element_Index(j),7),...
                PML_ElementConnectivity(PML_Element_Index(j),8),...
                PML_ElementConnectivity(PML_Element_Index(j),9));
end
fprintf(fid,'*UEL PROPERTY, elset=el-PML-type-1\n');
fprintf(fid,'2.08e+08, 0.3, 2000.0, 2., 25., 2., 0.00000001, 100.,\n');
fprintf(fid,'100., 100., 0.0, 0.0\n');
fclose(fid);



    