clc;clear;close all;

fid = fopen('Job-DRM-vefification-1.inp');
cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true ); 
cac = cac{1};

ElementConnectivity = str2num(char(cac(22377:44426)));
PML_ElementConnectivity = ElementConnectivity;
PML_ElementConnectivity(:,1) = PML_ElementConnectivity(:,1) + 1e5;

PML_Element_Set = str2num(char(cac(46592)));
PML_Element_Index = PML_Element_Set(1):PML_Element_Set(3):PML_Element_Set(2);

OutputName = sprintf('Files need to be imported into ABAQUS/PML_Info.txt');
fid = fopen(OutputName,'wt');
fprintf(fid, '*USER ELEMENT, NODES=4, TYPE=U1, PROPERTIES=11, COORDINATES=2, VARIABLES=80, unsym\n');
fprintf(fid, '1,2,21,22,23\n');
fprintf(fid,'*Element, type=U1, elset=el-PML-type-1\n');
for j = 1:length(PML_Element_Index)
    fprintf(fid,'%d, %d, %d, %d, %d\n',...
            PML_ElementConnectivity(PML_Element_Index(j),1),...
            PML_ElementConnectivity(PML_Element_Index(j),2),...
            PML_ElementConnectivity(PML_Element_Index(j),3),...
            PML_ElementConnectivity(PML_Element_Index(j),4),...
            PML_ElementConnectivity(PML_Element_Index(j),5));
end
fprintf(fid,'*UEL PROPERTY, elset=el-PML-type-1\n');
fprintf(fid,'2.08e+08, 0.3, 2000.0, 5., 5., 2., 0.00000001, 100.,\n');
fprintf(fid,'100., 0.0, 0.0\n');
fclose(fid);
    