% **********************************************************************
% *                                                                    *
% *    WHAT IT DOES:                                                   *
% *                                                                    *
% *    This function computes shape functions and their gradients for  *
% *    lineat quadrilateral element in parent coordinates.             *
% *                                                                    *
% *                                                                    *
% **********************************************************************
% *                                                                    *
% *   cee235b: Finite Element Analysis of Structures                   *
% *                                                                    *
% *                                                                    *
% *   History: <dd/mm/yyyy> {who} [what]                               *
% *                                                                    *
% *   <02,22,2014> {E. Esmaeilzadeh}                                   *
% *   [ original coding ]                                              *
% *                                                                    *
% **********************************************************************
function [N,B] = fun_ShapeFunction_3D(xi,eta,mu)

% N = [N1 N2 N3 N4 N5 N6 N7 N8] 
% shape functions for linear quadrilateral element
N1 = 1/8 * (1-xi)*(1-eta)*(1-mu);

N2 = 1/8 * (1+xi)*(1-eta)*(1-mu);

N3 = 1/8 * (1+xi)*(1+eta)*(1-mu);

N4 = 1/8 * (1-xi)*(1+eta)*(1-mu);

N5 = 1/8 * (1-xi)*(1-eta)*(1+mu);

N6 = 1/8 * (1+xi)*(1-eta)*(1+mu);

N7 = 1/8 * (1+xi)*(1+eta)*(1+mu);

N8 = 1/8 * (1-xi)*(1+eta)*(1+mu);

% B = [B1 B2 B3 B4 B5 B6 B7 B8]
% gradient of shape functions for linear quadrilateral element
B1 = [-1/8*(1-eta)*(1-mu);-1/8*(1-xi)*(1-mu);-1/8*(1-xi)*(1-eta)];

B2 = [+1/8*(1-eta)*(1-mu);-1/8*(1+xi)*(1-mu);-1/8*(1+xi)*(1-eta)];

B3 = [+1/8*(1+eta)*(1-mu);+1/8*(1+xi)*(1-mu);-1/8*(1+xi)*(1+eta)];

B4 = [-1/8*(1+eta)*(1-mu);+1/8*(1-xi)*(1-mu);-1/8*(1-xi)*(1+eta)];

B5 = [-1/8*(1-eta)*(1+mu);-1/8*(1-xi)*(1+mu);+1/8*(1-xi)*(1-eta)];

B6 = [+1/8*(1-eta)*(1+mu);-1/8*(1+xi)*(1+mu);+1/8*(1+xi)*(1-eta)];

B7 = [+1/8*(1+eta)*(1+mu);+1/8*(1+xi)*(1+mu);+1/8*(1+xi)*(1+eta)];

B8 = [-1/8*(1+eta)*(1+mu);+1/8*(1-xi)*(1+mu);+1/8*(1-xi)*(1+eta)];

N = [N1, N2, N3, N4, N5, N6, N7, N8]';

B = [B1, B2, B3, B4, B5, B6, B7, B8]';



