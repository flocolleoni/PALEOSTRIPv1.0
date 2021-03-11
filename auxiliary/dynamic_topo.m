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
% SEALEVEL CHANGES
%
% Retrieve sea level changes relative to present (dSL)
%
%
% Several options are available from the GUI interface:
%
% [1] Spatially uniform and constant prescribed dSL
%
% [2] Interpolated dSL from eustatic reconstructions from Haq et al.: values
% are retrieved at the ages of each decompacted layers (input parameters file)
%
% [3] Interpolated dSL (eustatic reconstruction) from Miller et al (2020):
% values are retrieved at the ages of each decompacted layers (input parameters file)
%
% [4] User-based provided file: may contain a timeseries (time dSL) or 2D
% horizontal maps of GIA-based sea level change
%
%
%  *****  VALUES PROVIDED NEED TO BE IN METERS *****
%  
%-------------------------------------------------------------------------------------

    function [dyntopo_corr,outputstring]=dynamic_topo(Np,N,Xcoord,Ycoord,ages,opt_dyn,dtopo_cte,dtopo_file,dtopo_dir,grid_opt)
       dyntopo_corr=zeros(Np,N);
       outputstring=' ';



       if (opt_dyn == 1) % Spatially uniform constant prescribed value
          for t=1:N
             dyntopo_corr(:,t) = dtopo_cte;
          end     
       elseif (opt_dyn == 2) % User-provided 1D timeseries file
           
             [DT, outputstring]=read_data_1d(dtopo_file,dtopo_dir,ages);
             for t=1:N
               dyntopo_corr(:,t)= DT(t);
             end
       elseif (opt_dyn == 3) % User-provided 2d/3D file

          [DT, outputstring]=read_dtopo(grid_opt,Np,N,Xcoord,Ycoord,ages,dtopo_file,dtopo_dir);
          for t=1:N
             dyntopo_corr(:,t)= DT(:,t);
          end    
       elseif (opt_dyn == 4) % No dyn. topo. corrections
         for t=1:N
            dyntopo_corr(:,t)= 0.;
         end
       end

       
end %END FUNCTION