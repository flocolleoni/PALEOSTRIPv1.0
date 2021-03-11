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
% Interpolates points between two prescribed coordinates from Extract panel
% of the GUI or from user-input files cointaining a transect of X,Y
% coordinates
%
% Interpolation is linear and uses scatteredInterpolant
%
% Number of points within the interpolated transect: 500 
% (except for extracting well for which number of points = 1)
%
% This number can be modified below (numnpoints)
%-----------------------------------------------------------------------------------

function [interpvar,Xnew,Ynew]=extract_pts(opt_extract,X,Y,Z,DataPoints)


Xmin=min(X); Xmax=max(X);
Ymin=min(Y); Ymax=max(Y);

% number of points within the interpolated transect 
%(except for well for which number of points = 1)
numnpoints = 500;

if (opt_extract==1) % Create profile from start-end points in GUI Table
        
    X0=DataPoints{1,1}; 
    Y0=DataPoints{1,2};
    X1=DataPoints{2,1}; 
    Y1=DataPoints{2,2};

    if (X0==X1 && Y0==Y1) % If only one points required    
        Xnew=X0; Ynew=Y0;
        numnpoints = 1;
    elseif (X0~=X1 && Y0==Y1)
        Xinc=(X1-X0)/numnpoints;
        Xnew=X0:Xinc:X1; Ynew=Y0;
        
    elseif (X0==X1 && Y0~=Y1)
        Yinc=(Y1-Y0)/numnpoints;
        Xnew=X0; Ynew=Y0:Yinc:Y1;
        
    elseif (X0~=X1 && Y0~=Y1)
        Xinc=(X1-X0)/numnpoints; Yinc=(Y1-Y0)/numnpoints;
        Xnew=X0:Xinc:X1; Ynew=Y0:Yinc:Y1;
    
    end

else  % Create profile from points read in user-provided file
     Xnew=DataPoints(:,1);
     Ynew=DataPoints(:,2);
   
    numnpoints=length(DataPoints(:,1));
end

% Interpolate data
    
    interpvar=zeros(numnpoints,1);

    F = scatteredInterpolant(X,Y,Z,'linear','none');
    interpvar=F(Xnew,Ynew);


end

