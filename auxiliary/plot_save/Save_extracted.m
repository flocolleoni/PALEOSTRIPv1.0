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
%%% SAVE EXTRACTED VARIABLES TO FILE %%%
% --- Executes on button press in pushbutton_save_extract.
%
% Uses "grid_opt" to save 1D well, 2D transect to ASCII files
%
% At each timestep, the ASCII file is written from top to bottom as follows
% time1: [X Y Z1 Z2 Z3...Zn]
% time2: [X Y Z2 Z3...Zn]
% time3: [X Y Z3...Zn]
%
% and so on. 
% Where Z1 corresponds to the uppermost layer and Zn corresponds to the
% basement
%
% Decompacted and corrected horizon depths corresponds to:
% depths of horizons Z1 to Zn at time 1
% depths of horizons Z2 to Zn at time 2
% .
% .
% .
% and so on.
%-----------------------------------------------------------------------------------
function Save_extracted(plot_extracted,destinationfolder,extractid,extract_prof,Xpts,Ypts,profiles)

  
% RUNID String
ext='.dat';
    

if (plot_extracted < 7 || plot_extracted == 8 )
       IMAX=length(extract_prof);
else % BETA FACTOR AND THERMAL SUBSIDENCE
       IMAX=size(extract_prof,2);
end

% SELECT VARIABLE TO BE SAVED
if (plot_extracted == 2) %'Horizons depths'
             string_run = strcat('_backstripped_depths_extracted_',extractid,ext);
elseif (plot_extracted == 3) %,'Isopacks'
             string_run = strcat('_backstripped_isopachs_extracted_',extractid,ext);
elseif (plot_extracted == 4) %'Density'
             string_run = strcat('_backstripped_density_extracted_',extractid,ext);
elseif (plot_extracted == 5) %'porosity'
             string_run = strcat('_backstripped_porosity_extracted_',extractid,ext);
elseif (plot_extracted == 6) %'THERMAL SUBSIDENCE'
                string_run = strcat('_thermal_subs_correction_extracted_',extractid,ext);
elseif (plot_extracted == 7) %'ISOSTASY'
             string_run = strcat('_Isostatic_correction_extracted_',extractid,ext);
elseif (plot_extracted == 8) %'Beta Factor'        
             string_run = strcat('Beta_Factor_extracted_',extractid,ext);
end
 


if (plot_extracted < 6 || plot_extracted == 7 )    
    for ii=length(extract_prof):-1:1 % for each time 
             
             filename=[profiles(ii).layername string_run];
             fpath=fullfile(destinationfolder,filename);
             [fid]=check_save_to_file(fpath,destinationfolder,filename);
            
             to_be_written=fliplr(extract_prof(ii).interface);
             formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
             fprintf(fid,formatstring,[Xpts, Ypts, to_be_written]');

             fclose(fid);
    end
       
    
elseif (plot_extracted == 6) % THERMAL SUBSIDENCE

        for ii=1:size(extract_prof,2)
             to_be_written=fliplr(extract_prof(:,ii));
             filename=[profiles(ii).layername string_run];
             fpath=fullfile(destinationfolder,filename);
            
             [fid]=check_save_to_file(fpath,destinationfolder,filename);         

             formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
             fprintf(fid,formatstring,[Xpts, Ypts, to_be_written]');
             fclose(fid);
        end
          
elseif (plot_extracted == 8) %'Beta Factor'        
             to_be_written=extract_prof;
             filename=string_run; 
             fpath=fullfile(destinationfolder,filename);
            
             [fid]=check_save_to_file(fpath,destinationfolder,filename);         

            %CHOOSE GRID DIMENSION
             formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
             fprintf(fid,formatstring,[Xpts,Ypts, to_be_written]');
             fclose(fid);
           
end
     

end