%% Copyright  (C)  2012  Gunther Struyf <https://github.com/GuntherStruyf>
%% Version: 1.0
%% Author: Gunther Struyf <https://github.com/GuntherStruyf>
%% URL: https://github.com/GuntherStruyf
%
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
%%
clear
close all
clc

%% Configuration

% environment:
g	= 9.81;	% [m/s^2]			% gratitational constant
m	= 10;	% [kg]				% point mass load of 10kg
G	= 10*g;	% [N]				% gravitational force resulting from the point load
qdistributed = 0; % [N/m]		% distributed force over the beam's length

% beam characteristics
E	= 70e3;	% [N/mm^2 = MPa]	% Elasticity modulus (Young's modulus)
E	= E*1e6;% [N/m^2 = Pa]
[Jxx yd] = beam_characteristics('I',0.04,0.05,0.005);
L	= 2;	% [m]				% Length of the beam
xF	= 1.5;	% [m]				% position where the force will be exerted

%% Preparation: filling input arguments of beamFEM
xnode	= [	0		xF		L	];
F		= [	NaN		G		NaN	];
u		= [	0		NaN		0	];
M		= [	0		0		0	];
theta	= [	NaN		NaN		NaN	];

%% Execution
[xnode, unode, tnode, x, u, Qs, Mb, stress ] = beamFEM(xnode, F, M, u, theta, qdistributed, Jxx, E, yd);

%% Visualization

figure
plot(x,-u);

hold on;
plot(xlim,[0 0],'k--');
xlabel('x [m]');
ylabel('deflection [m]');


