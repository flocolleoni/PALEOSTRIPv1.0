%-----------------------------------------------------------------------------------
% Copyright (C) 2021 PALEOSTRIP Authors
% 
% This file is part of PALEOSTRIP.
% 
% PALEOSTRIP is free software; you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
% 
% PALEOSTRIP is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License
% along with PALEOSTRIP; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%-----------------------------------------------------------------------------------
%
% Interpolate input var (sea level and dynamic topography 2D maps) onto:
%
% [1] An expanded data mesh (rectangular regular grid based on original pologyn
% coordinates) which is expanded by 30% outside of all edges. Loads are 
% interpolated with natural neighbour algorithm and linearily extrapolated where 
% data are missing. Horizontal resolution from original input sediment data is 
% preserved (and must be regular in x and y).
%
% [2] Newly created expanded mesh is placed at the center of a second
% mesh, twice as large as the first expanded mesh, with similar horizontal
% resolution
%
%-----------------------------------------------------------------------------------
function [var_interp]=mesh3D_inputvar(X,Y,Xsl,Ysl,var1)
% CREATE 3D MESH FOR 3D FLEXURE

% Define the 2D-grid

    % GET UNIQUE set of X and Y
    Xorig=unique(X); 
    Yorig=unique(Y);    
    
    % get dx and dy - Needs to be regular spacing
    dx=abs(Xorig(2)-Xorig(1));
    dy=abs(Yorig(2)-Yorig(1));    
    
    % get Xmin-Xmax and Ymin-Ymax
    Xmin=min(Xorig); Xmax=max(Xorig);
    Ymin=min(Yorig); Ymax=max(Yorig);   
    
 % Expand data to fill a regular rectangular domain
    
    [Xorig2d Yorig2d] = meshgrid(Xmin:dx:Xmax,Ymin:dy:Ymax);
     %Interpolate EET on original data grid

    F = scatteredInterpolant(Xsl,Ysl,var1,'natural','linear');
    Q_interp=F(Xorig2d,Yorig2d);    

    
% Fill the nodes with similar coordinates with thickness values (THICK)
    % get index of original coordinates in the new 2D-grid
    [~, index] = ismember([X Y] , [Xorig2d(:) Yorig2d(:)], 'rows' );
    

    %FLEX3D WANTS 1D VECTOR AND NOT 2D MATRIX
    var_interp(:)=Q_interp(index(:)); %Fill the new vector with value with identical coordinates
      
    
end         
