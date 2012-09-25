%% BEAM_CHARACTERISTICS 1-dimensional beam characteristics
%
%INPUTS:
%	beamtype
%			'I'	for an I-beam.
%			'L'	for an L-beam.
%			'T'	for an t-beam.
%			'U'	for an U-beam.
%			Other types will raise an error.
%
%	B		Width of the beam [m].
%	H		Height of the beam [m].
%	e		Thickness of the beam sections [m].
%			see http://en.wikipedia.org/wiki/Second_moment_of_area#Second_moment_of_area_for_various_cross_sections
%			for more information
%
%OUTPUTS:
%	Jxx		Second moment of area for the specified beam [m^4].
%	yd_max	Distance from the neutral line to both top and bottom of the beam [m].
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
function [ Jxx, yd_max ] = beam_characteristics( beamtype, B,H,e )
	% neutral axis = yg: yg*Atot = sum(Ai*yi)

	switch lower(beamtype)
		case {'t'}
			b	= [B-e	e];
			h	= [e	H];
			yi	= [B-e/2	H/2];
			sgn	= [1	1];
		case {'l'}
			b	= [B-e	e];
			h	= [e	H];
			yi	= [e/2	H/2];
			sgn	= [1	1];
		case 'i'
			b	= [B	B-e];
			h	= [H	H-2*e];
			yi	= [H/2	H/2];
			sgn	= [1	-1];
		case 'u'
			b	= [B-2*e e];
			h	= [e	H];
			yi	= [e/2	H/2];
			sgn	= [1	2];
		otherwise
			error('Unsupported beamtype, only L,T,I and U are valid beam types.');
	end

	yg	= neutral_axis(b,h,yi,sgn);
	Jxx = Jgen(b,h,yi-yg,sgn);

	yd_max = [-yg H-yg];

end

function yg = neutral_axis(b,h,yi,sgn)
	A			= b.*h.*sgn;
	empty_idx	= A==0;
	if all(empty_idx)
		error('Area of zero surface area has no neutral axis.');
	end
		
	A(empty_idx)	= [];
	yi(empty_idx)	= [];
	
	if all(yi==yi(1))
		yg=yi(1);
	else
		Atot	= sum(A);
		yg		= sum(A.*yi)/Atot;
	end
end
function Jxx = Jgen(b,h,y,sgn)
	A		= b.*h;
	Jloc	= b.*h.^3/12;
	Jxx		= sum(sgn.*(Jloc + y.^2.*A));
end

