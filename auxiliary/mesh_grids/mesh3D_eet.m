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
% Interpolate input spatially variable Effective Elastic Thickness file
% onto:
%
% [1] An expanded data mesh (rectangular regular grid based on original pologyn
% coordinates) which is expanded by 30% outside of all edges. EET input data are interpolated 
% with natural neighbour algorithm and linearily extrapolated where data are missing. 
% Horizontal resolution from original input sediment data is preserved 
% (and must be regular in x and y).
%
% [2] Newly created expanded mesh is placed at the center of a second
% mesh, twice as large as the first expanded mesh, with similar horizontal
% resolution
%
%-----------------------------------------------------------------------------------

function [EET_1d,EET_out]=mesh3D_eet(X,Y,Xeet,Yeet,EET)
% CREATE 3D MESH FOR 3D FLEXURE

 % ORIGINAL GRID DIMENSIONS
 
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
    Yorig2d=Ymin:dy:Ymax ; NY= length(Yorig2d);
    Xorig2d=Xmin:dx:Xmax ; NX= length(Xorig2d);
    
    Xmin2=Xmin-(ceil(NX*0.3)*dx); Xmax2=Xmax+(ceil(NX*0.3)*dx);
    Ymin2=Ymin-(ceil(NY*0.3)*dy); Ymax2=Ymax+(ceil(NY*0.3)*dy); 
     
    [Xorig2d Yorig2d] = meshgrid(Xmin2:dx:Xmax2,Ymin2:dy:Ymax2);   

 % Interpolate EET on original sediment data grid
    F = scatteredInterpolant(Xeet,Yeet,EET,'natural','linear');%'natural','linear'
    EET_interp=F(Xorig2d,Yorig2d); %    EET_interp=F(X,Y);

    % get index of original coordinates in the new 2D-grid
    [~, index] = ismember([X Y] , [Xorig2d(:) Yorig2d(:)], 'rows' );
    EET_out(1:length(X)) = EET_interp(index(:)); % For plotting initial conditions
    
% Create 2D-grid twice as large as the expanded meshgrid domain (30%) for flexural isostasy
     NX_orig2d = length(Xorig2d);
     NY_orig2d = length(Yorig2d);

     Xemin=Xmin-(ceil(NX_orig2d*1)*dx); Xemax=Xmax+(ceil(NX_orig2d*1)*dx);
     Yemin=Ymin-(ceil(NY_orig2d*1)*dy); Yemax=Ymax+(ceil(NY_orig2d*1)*dy);
     
     [X2d Y2d] = meshgrid(Xemin:dx:Xemax,Yemin:dy:Yemax);
      NY2= length(Y2d);
      NX2= length(X2d);
      
     % get index of similar coordinates in the new 2D-grid
     [~, index_2d] = ismember([Xorig2d(:) Yorig2d(:)] , [X2d(:) Y2d(:)], 'rows' );

     
% OUTPUT FOR FLEX3D (takes vector not matrix)
    EET_2d=zeros(NX2,NY2);
    EET_1d= EET_2d(:);

    EET_1d(index_2d(:))= EET_interp(:); %Fill the new vector with value with identical coordinates
    EET_2d = reshape(EET_1d, [NY2 NX2]); % Transform into 2D matrix
    
   % Populate the 1d vector with right dimensions order
    k=1;
    for i=1:NY2
        for j=1:NX2
            EET_1d(k)=EET_2d(i,j);
            k=k+1;
        end
    end
    
end         
