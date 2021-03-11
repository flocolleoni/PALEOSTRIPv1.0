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

    function [SL,outputstring]=sealevel(Np,N,opt_sl,ds,ages,Xcoord,Ycoord,SLfile,SLdir)
       SL=zeros(Np,N);
       outputstring=' ';

       data_folder = 'data';
       addpath(genpath('./read_data'));

       for t=1:N
           if (opt_sl == 1) % Spatially uniform constant prescribed value
                 SL(:,t) = ds;
                 
           elseif (opt_sl == 2) % interpolated from timeseries
                 file_sl='haq_sealevel_curve.dat';
                 [dSL, outputstring]=read_data_sl_1d(file_sl,data_folder,ages);

                 SL(:,t) = dSL(t);
                 
           elseif (opt_sl == 3)% interpolated from timeseries
                 file_sl='Miller_et_al_2020_0-66Ma.dat';

                 [dSL, outputstring]=read_data_sl_1d(file_sl,data_folder,ages);

                 SL(:,t) = dSL(t);
                 
           elseif (opt_sl == 4) % User-provided file

                 [dSL, outputstring]=read_data_sl_1d(SLfile,SLdir,ages);
                 
                 SL(:,t)= dSL(t);
                  
           end    
       end
       
       % Spatially variable sea level maps
       % Note that there is no possibility to provide time-evolving 2D maps
       % for now
       if (opt_sl == 5) % User-provided file

         [dSL, outputstring]=read_data_sl_2d(Xcoord,Ycoord,SLfile,SLdir);

         for i=1:N
            SL(:,t)= dSL(:);
         end

       elseif (opt_sl == 6) % No sea level corrections
         for i=1:N
            SL(:,t)= 0.;
         end
       end
       
end %END FUNCTION