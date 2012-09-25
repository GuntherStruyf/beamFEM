%
%% Copyright  (C)  2012  Gunther Struyf <https://github.com/GuntherStruyf>
%% Version: 1.0
%% Author: Gunther Struyf <https://github.com/GuntherStruyf>
%% URL: https://github.com/GuntherStruyf

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%
function [x, u] = Nx2u(u_i, theta_i, L, Npoints)
%% [x, u] = Nx2u(u1, theta1, u2, theta2, Npoints)
%Convert the parametric solution [u_i theta_i] to real values.
%The used form functions are:
%	N1(xi) =  2*xi^3 - 3*xi^2 + 1
%	N2(xi) = L*(xi^3 - 2*xi^2 + xi)
%	N3(xi) = -2*xi^3 + 3*xi^2
%	N4(xi) = L*(xi^3 -   xi^2)
%
%INPUT:
%	u_i			nNode x 1 vector
%					Displacement at nodes.
% 	theta_i		nNode x 1 vector
%					Rotation at nodes.
% 	L			(nNode-1) x 1 vector
%					Length of each element.
% 	Npoints		Number of interpolating points in each element.
%

	nNode = numel(u_i);
	nElem = nNode-1;

	if nNode~=numel(theta_i)
		error('BeamFEM:ArgumentMismatch','u and theta must be of same length');
	elseif nElem~=numel(L)
		if isscalar(L)
			L = L(ones(nElem,1));
		else
			error('BeamFEM:ArgumentMismatch','length of L does not match with length of u and theta');
		end
	end

	xi 	= linspace(0,1,Npoints+1);
	Nu1	= 2*xi.^3 -3*xi.^2+1;
	Nt1b=  (xi.^3 -2*xi.^2+xi);
	Nu2	=-2*xi.^3 +3*xi.^2;
	Nt2b=  (xi.^3   -xi.^2);

	x = NaN(1,Npoints*nElem+1);
	u = NaN(1,Npoints*nElem+1);

	x0 = 0;
	for ii = 1:nElem
		idxs = (ii-1)*Npoints + (1:Npoints+1);

		x(idxs) = x0 + xi*L(ii);
		u(idxs) =	Nu1	* u_i(ii)		+ ...
					Nu2	* u_i(ii+1)	+ ...
					Nt1b*L(ii)	* theta_i(ii) +...
					Nt2b*L(ii)	* theta_i(ii+1);

		x0 = x0+L(ii);
	end
end