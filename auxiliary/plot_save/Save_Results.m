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
%%% SAVE VARIABLES TO FILE %%%
% --- Executes on button press in TpushbuttonSaveResults.
%
% Uses "grid_opt" to save 1D well, 2D transect or 3D maps to ASCII files
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

function Save_Results(M,save_menu,opt_grid,runid,destinationfolder,profiles,isopachs,density,porosity,tsubs_corr,isostasy,tectonic) 
 % File extension
 ext='.dat';

 if (save_menu == 1)
        errordlg('Select a variable to be saved', 'Error', 'modal');
        return 
 end

if (save_menu == 2) %'Horizons depths'
       VAR=profiles;
       IMAX=length(profiles);
       string_run = strcat('_backstripped_depths_',runid,ext);

elseif (save_menu == 3) %,'Isopacks'
       VAR=isopachs;
       IMAX=length(isopachs);
       string_run = strcat('_backstripped_isopachs_',runid,ext);

elseif (save_menu == 4) %,'Density'
       VAR=density;
       IMAX=length(density) ;
       string_run = strcat('_backstripped_density_',runid,ext);

elseif (save_menu == 5) %,'Porosity'
       VAR=porosity;
       IMAX=length(porosity) ;
       string_run = strcat('_backstripped_porosity_',runid,ext);

elseif (save_menu == 6) %,'Thermal Subsidence'
       VAR=tsubs_corr;
       IMAX=size(tsubs_corr,2);
       string_run = strcat('_thermal_subs_correction_',runid,ext);

elseif (save_menu == 7) %,'Isostatic correction'
       VAR=isostasy;
       IMAX=length(isostasy);
       string_run = strcat('_isostatic_correction_',runid,ext);
       
elseif (save_menu == 8) %,'Tectonic subsidence'
       VAR=tectonic;
       IMAX=size(tectonic,2);
       string_run = strcat('_tectonic_subsidence_',runid,ext);
       
end 

          
% SAVE VARIABLE TO FILE
for ii=1:IMAX % for each time
   if (save_menu < 6 || save_menu == 7 )
            
            if (save_menu == 2) % horizons depth
             filename=[profiles(ii).layername string_run];
            else % all other variables
             filename=[profiles(ii+1).layername '_' profiles(ii).layername string_run];  
            end
            
            fpath=fullfile(destinationfolder,filename);
            
            [fid]=check_save_to_file(fpath,destinationfolder,filename);
            
            to_be_written=fliplr(VAR(ii).interface);

            %CHOOSE GRID DIMENSION
             if (opt_grid == 1)
                 fprintf(fid,[profiles(ii).layername '  and deeper \n']); % labels
                 formatstring=[repmat(' %4.4f',1,size(to_be_written(1,:),2)) '\n'];
                 fprintf(fid,formatstring,[to_be_written(1,:)]');
                 fclose(fid);

             elseif (opt_grid == 2)
                 fprintf(fid,[' X (m) ' profiles(ii).layername '  and deeper \n']); % labels
                 formatstring=['%6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 fprintf(fid,formatstring,[M(:,1), to_be_written]');
                 fclose(fid);

             elseif(opt_grid == 3)
                 fprintf(fid,[' X (m) Y(m) ' profiles(ii).layername '  and deeper \n']); % labels
                 formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 fprintf(fid,formatstring,[M(:,1), M(:,2), to_be_written]');
                 fclose(fid);
             end
            
  elseif (save_menu == 6) %'THERMAL SUBSIDENCE'

             to_be_written=fliplr(VAR(:,ii));
             filename=['Present-day_' profiles(ii+1).layername string_run];
             fpath=fullfile(destinationfolder,filename);
            
             [fid]=check_save_to_file(fpath,destinationfolder,filename);         

            %CHOOSE GRID DIMENSION
             if (opt_grid == 1)
                 fprintf(fid,[profiles(ii).layername '\n']); % labels
                 formatstring=[repmat(' %4.4f',1,size(to_be_written(1,:),2)) '\n'];
                 fprintf(fid,formatstring,[to_be_written(1,:)]');
                 fclose(fid);

             elseif (opt_grid == 2)
                 fprintf(fid,[' X (m)  ' profiles(ii).layername '\n']); % labels
                 formatstring=['%6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 fprintf(fid,formatstring,[M(:,1), to_be_written]');
                 fclose(fid);

             elseif(opt_grid == 3)
                 fprintf(fid,[' X (m)   Y (m)  ' profiles(ii).layername '\n']); % labels
                 formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 fprintf(fid,formatstring,[M(:,1), M(:,2), to_be_written]');
                 fclose(fid);
             end
             
  elseif (save_menu == 8) %'TECTONIC SUBSIDENCE HISTORY'

             to_be_written=fliplr(VAR(:,ii));
             filename=['Basement_at_' num2str(profiles(ii).age) '_Ma' string_run];
             fpath=fullfile(destinationfolder,filename);
            
             [fid]=check_save_to_file(fpath,destinationfolder,filename);         

            %CHOOSE GRID DIMENSION
             if (opt_grid == 1)
                 fprintf(fid,[profiles(ii).layername '\n']); % labels
                 formatstring=[repmat(' %4.4f',1,size(to_be_written(1,:),2)) '\n'];
                 fprintf(fid,formatstring,[to_be_written(1,:)]');
                 fclose(fid);

             elseif (opt_grid == 2)
                 fprintf(fid,[' X (m)  ' profiles(ii).layername '\n']); % labels
                 formatstring=['%6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 fprintf(fid,formatstring,[M(:,1), to_be_written]');
                 fclose(fid);

             elseif(opt_grid == 3)
                 fprintf(fid,[' X (m)   Y (m)  ' profiles(ii).layername '\n']); % labels
                 formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 fprintf(fid,formatstring,[M(:,1), M(:,2), to_be_written]');
                 fclose(fid);
             end
  end
end