%-----------------------------------------------------------------------------------
% Copyright (C) 2020 PALEOSTRIP Authors
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
function [Load_1d,NX,NY,dx,dy,X2d,Y2d,index]=mesh3D(X,Y,Q)
% CREATE 3D MESH FOR 3D FLEXURE (TAFI coupling)

% Define the 2D-grid

    % GET UNIQUE set of X and Y
    Xorig=unique(X); 
    Yorig=unique(Y);

    % get dx and dy - Needs to be regular spacing
    dx=Xorig(2)-Xorig(1);
    dy=Yorig(2)-Yorig(1);

    % get Xmin-Xmax and Ymin-Ymax
    Xmin=min(Xorig); Xmax=max(Xorig);
    Ymin=min(Yorig); Ymax=max(Yorig);

    % create 2D-grid
    Ynew=Ymin:dy:Ymax ; NY= length(Ynew);
    Xnew=Xmin:dx:Xmax ; NX= length(Xnew);

   
    
    [X2d Y2d] = meshgrid(Xmin:dx:Xmax,Ymin:dy:Ymax);
                
% Fill the nodes with similar coordinates with thickness values (THICK)
    
    % get index of similar coordinates in the new 2D-grid
    [~, index] = ismember([X Y] , [X2d(:) Y2d(:)], 'rows' );

    Load_2d=zeros(NX,NY);
    Load_1d= Load_2d(:);          %Get 2D Matrix to 1D

    %TAFI WANTS 1D VECTOR AND NOT 2D MATRIX
    
    Load_1d(index(:))=Q(:); %Fill the new vector with value with identical coordinates
    %Load_2d = reshape(Load_1d, [NY NX]); % Transform into 2D matrix to pass to TAFI     
     
%    figure
%    mesh(X2d,Y2d,Load_2d)
%    hold on
%    plot3(X,Y,Q,'o')
    
end         
