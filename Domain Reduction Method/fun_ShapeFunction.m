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
function [N,B] = fun_ShapeFunction(xi,eta)

% N = [N1 N2 N3 N4] 
% shape functions for linear quadrilateral element
N1 = 1/4 * (1-xi)*(1-eta);

N2 = 1/4 * (1+xi)*(1-eta);

N3 = 1/4 * (1+xi)*(1+eta);

N4 = 1/4 * (1-xi)*(1+eta);

% B = [B1 B2 B3 B4]
% gradient of shape functions for linear quadrilateral element
B1 = [-1/4*(1-eta);-1/4*(1-xi)];

B2 = [+1/4*(1-eta);-1/4*(1+xi)];

B3 = [+1/4*(1+eta);+1/4*(1+xi)];

B4 = [-1/4*(1+eta);+1/4*(1-xi)];

N = [N1, N2, N3, N4]';

B = [B1, B2, B3, B4]';



