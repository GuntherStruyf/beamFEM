%% BEAMFEM Finite elements method applied to a 1-dimensional beam
%
%INPUTS:
%	xnode	Nnx1 vector of x-positions at which there are nodes.
%
%	F		Nnx1 vector of given external forces acting on the
%			corresponding nodes. If there is no (known) external force at
%			a node, F should be NaN (0N is still a known force!)
%
%	M		Nnx1 vector of given external moments acting on a corresponding
%			nodes. If there is no (known) external force at a node, F
%			should be NaN.
%
%	u		Nnx1 vector of given displacements. No displacement -> NaN.
%			A node's displacement cannot be set in conjunction with a 
%			force F in the same node!
%
%	theta	Nnx1 vector of given rotations. No rotation -> NaN.
%			A node's rotation cannot be set in conjunction with a 
%			moment M in the same node!
%
%	qelem	(Nn-1)x1 vector indicating the value of the distributed load
%			per element.
%			If no distributed force is applied, 0 is also a valid value.
%
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
function [xnode, unode, tnode, x, u, Qs, Mb, stress ] = beamFEM( xnode, F,M,u,theta,qelem, Jxx, E, max_dist_from_neutral)
	%% check input
	if any(~isnan(F) & ~isnan(u))
		error('BeamFEM:ArgumentMismatch','F and u are not complementary');
	elseif any(~isnan(M) & ~isnan(theta))
		error('BeamFEM:ArgumentMismatch','M and theta are not complemenary');
	end

	%% initialize output
	% xnode
	Qs=0;
	Mb=0;
	stress=0;

	%% STEP 1: Initialization

	% initialization of the stiffness matrix K, net external force F and the
	% parameters d to solve for (displacement & rotation in each node).
	%    nElem = number of elements
	%    nNode = number of nodes
	%	 nnz_K = number of non-zero elements in K
	%    K     = total stiffness matrix
	%    F     = net external force (right hand side of the equation)
	%    d     = solution parameters

	nNode	= numel(xnode);
	nElem	= nNode-1;
	nGauss	= 2;% number of Gauss quadrature points

	nDOF	= 2*nNode;
	nnz_K	= 12*nNode+4;
	K		= spalloc( nDOF, nDOF, nnz_K );
	f		= zeros( nDOF, 1 );
	d		= zeros( nDOF, 1 );

	L		= diff(xnode);

	if numel(E)~=nElem
		if isscalar(E)
			E = E(ones(nElem,1));
		else
			error('BeamFEM:ArgumentMismatch','E is of different size as the number of elements');
		end
	end
	if numel(qelem)~=nElem
		if isscalar(qelem)
			qelem = qelem(ones(nElem,1));
		else
			error('BeamFEM:ArgumentMismatch','qelem is of different size as the number of elements');
		end
	end

	EI		= E.*Jxx;


	%% Calculation of integration points
	%    xi  = normalized x-coordinates of Gaussian points ( 0 <= xi  <=1 )
	%    H   = weighting factors of the Gaussian points
	switch nGauss
		case 1, xi = 1 / 2;                                     H = [ 1 ];
		case 2, xi = [ 1-sqrt(1/3), 1+sqrt(1/3) ] / 2;          H = [ 1, 1 ] / 2;
		case 3, xi = [ 1-sqrt(3/5), 1, 1+sqrt(3/5) ] / 2;       H = [ 5, 8, 5 ] / 18;
		otherwise, [XI,H]=Gauss1D(nGauss);  xi=(1+XI)/2; H=H/2;
	end


	% ----------------------------------------------------------------------------- %
	%% SETUP TOTAL STIFFNESS MATRIX K
	% ----------------------------------------------------------------------------- %

	% loop over all elements
	for pp = 1:nElem

		% initalization of the element stiffness matrix Ke_glo in a global reference coordinate system
		Ke_glo = zeros( 4,4 );

		% loop over Gaussian points
		for gg = 1:nGauss
			% compute the 2nd order derivative of the form function N
			% to the normalized coordinates xi for Gaussian point xi(gg)
			% N1(xi) =  2*xi^3 - 3*xi^2 + 1
			% N2(xi) = L*(xi^3 - 2*xi^2 + xi)
			% N3(xi) = -2*xi^3 + 3*xi^2
			% N4(xi) = L*(xi^3 -   xi^2)
			% see Course Numeric Modelling for Mechanics, p. 130
			dN2dxi2 = [12*xi(gg)-6  ;  L(pp)*(6*xi(gg)-4) ;  -12*xi(gg)+6  ;   L(pp)*(6*xi(gg)-2)];

			% Compute the B-matrix for Gaussian point xi(gg)
			% see Course Numeric Modelling for Mechanics, Exercices p. 19-20
			B = dN2dxi2./L(pp)^2;
			% Integrate B*EI*B'*L by Gauss-Legendre approximation
			I_Galerkin = H(gg)*(B*EI(pp)*B'*L(pp));
			Ke_glo = Ke_glo + I_Galerkin;
		end

		% assembly of total stiffness matrix
		%   idx = row numbers of element degrees of freedom in d
		idx = (pp-1)*2+1;
		K(idx:idx+3,idx:idx+3) = K(idx:idx+3,idx:idx+3)   +   Ke_glo;
	end

	% visualization of the total stiffness matrix
	% figure;spy(K);

	% ----------------------------------------------------------------------------- %
	%% SETUP OF RIGHT HAND SIDE OF STIFFNESS EQUATION
	% ----------------------------------------------------------------------------- %

	% The unknown reaction forces F_c are initially set to 0. They have no influence
	% on the calculation of beam deformation, because the corresponding rows fall out
	% during partitioning.
	%    idx = true if the row corresponds with a given external inputs.

	F_bi		= ~isnan(F);
	M_bi 		= ~isnan(M);
	indexF		= find(F_bi)*2-1;
	indexM		= find(M_bi)*2;
	f(indexF)	= F(F_bi);
	f(indexM)	= M(M_bi);

	% add offset for dealing with distributed loads
	for pp=1:nElem
		idxs = (pp-1)*2+(1:4);
		f(idxs) = f(idxs) + qelem(pp)*L(pp)/12*[6	L(pp)	6	-L(pp)]';
	end

	% ----------------------------------------------------------------------------- %
	%% PARTITIONING
	% ----------------------------------------------------------------------------- %

	% index_c	= degrees of freedom in d, that are constrained
	% index_f	= degrees of freedom in d, that are unknown
	% K_ff		= part of the stiffness matrix associated with the unknown degrees of freedom
	% K_cf		= part of the stiffness matrix needed for computation of reactive forces
	% f_f		= part of the right hand side of the equation with given load forces

	index_c = find(~isnan([u(:)  theta(:)]'));
	index_f = setdiff( 1:nDOF, index_c );
	K_ff	= K(index_f,index_f);
	K_cf	= K(index_c,index_f);
	f_f		= f(index_f,1);

	K_fc	= K(index_f,index_c);

	% ----------------------------------------------------------------------------- %
	%% SOLVING THE STIFFNESS EQUATION
	% ----------------------------------------------------------------------------- %

	% Solve the system of equations to the unknown degrees of freedom
	d(index_f) 	= K_ff\(f_f-K_fc*d(index_c));

	% Compute the reaction forces
	f = K*d;

	[x,u]=Nx2u(d(1:2:end), d(2:2:end), L, 50);


	unode = d(1:2:end);
	tnode = d(2:2:end);
	return;
	
	% ----------------------------------------------------------------------------- %
	%% POSTPROCESSING
	% ----------------------------------------------------------------------------- %
	% loop over the elements
	for pp = 1:nElem
		% set up the element stiffness matrix Ke_loc in the local reference coordinate system

		Ke_loc = E(pp)*A(pp)/L(pp) * [	1	0	-1	0;
										0	0	0	0;
										-1	0	1	0;
										0	0	0	0];

		% set up the transformation matrix
		% theta(pp) is already calculated ^^
		Am = [	cos(theta(pp))	-sin(theta(pp));
				sin(theta(pp))	cos(theta(pp))];
		G(1:2,1:2) = Am;
		G(3:4,3:4) = Am;

		p1id 	= Element(pp,3);
		p2id 	= Element(pp,4);
		dglo	= d([p1id*2-1	p1id*2	p2id*2-1	p2id*2]);
		dloc	= G' *dglo;
		Nloc 	= Ke_loc * dloc;
		RodForce(pp,1:2) = Nloc(1:2);
	end

	% ----------------------------------------------------------------------------- %
	% EXPORT VAN RESULTAAT
	% ----------------------------------------------------------------------------- %

	% maximale vertikale verplaatsing en hoekverdraaiing
	[C,I]=max(abs(d(1:2:end)));
	maxu = d(I*2-1)

	[C,I]= max(abs(d(2:2:end)));
	maxtheta = d(I*2)*180/(2*pi)

	% reactiekrachten
	%$$$

	% visualiseer vervormingspatroon
	EEM_balk_plot( L_tot, nElem, d );
end

