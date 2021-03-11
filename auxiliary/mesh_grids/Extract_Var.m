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
% Extract a transect of the selected variable between two prescribed points
% in the Table from the GUI or from user-input files cointaining a transect of X,Y
% coordinates
%
%-----------------------------------------------------------------------------------

function [extract_prof,Xpts,Ypts]=Extract_Var(opt_extract,plot_extracted,DataPoints,X,Y,profiles,isopachs,density,...
    porosity,tsubs_corr,isostasy,betaF,tectonic)


if (plot_extracted == 2) %'Horizons depths'
                   VAR=profiles;
                   IMAX=length(profiles);
            elseif (plot_extracted == 3) %,'Isopachs'
                   VAR=isopachs;
                   IMAX=length(isopachs); 
            elseif (plot_extracted == 4) %,'Density'
                   VAR=density;
                   IMAX=length(density) ;
            elseif (plot_extracted == 5) %,'Porosity'
                   VAR=porosity;
                   IMAX=length(porosity) ;
            elseif (plot_extracted == 6) %,'Thermal Subsidence'
                   VAR=-tsubs_corr;
                   IMAX=size(tsubs_corr,2);
            elseif (plot_extracted == 7) %,'Isostatic correction'
                   VAR=isostasy;
                   IMAX=length(isostasy);
            elseif (plot_extracted == 8) %,'Beta Factor'
                   VAR=betaF';
                   IMAX=size(betaF,1);
            elseif (plot_extracted == 9) %,'Tectonic'
                   VAR=tectonic';
                   IMAX=size(tectonic,2); 
end  
            
            if (plot_extracted < 6 || plot_extracted == 7 )
            for ii=1:IMAX
                 for jj=IMAX-ii+1:-1:1
                     to_be_extracted = VAR(ii).interface(:,jj);
                     [extracted,Xpts,Ypts]=extract_pts(opt_extract,X,Y,to_be_extracted,DataPoints);
                     extract_prof(ii).interface(:,jj)=extracted;
                 end
            end
             
            else % IN THE CASE OF BETA FACTOR & THERMAL SUBSIDENCE
            
                 for ii=1:IMAX
                     to_be_extracted = VAR(:,ii);
                     [extracted,Xpts,Ypts]=extract_pts(opt_extract,X,Y,to_be_extracted,DataPoints);
                     extract_prof(:,ii)=extracted;
                 end    
                
                
            
            end % END PLOT_EXTRACTED MENU
end

