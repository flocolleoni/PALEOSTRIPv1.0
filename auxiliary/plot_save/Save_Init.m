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
function Save_Init(M,plot_init,opt_grid,runid,destinationfolder,profiles,depth,EET_init,SL_init,DynTopo_init,betaF) 
% File extension
 ext='.dat';

 if (plot_init == 1)
        errordlg('Select a variable to be saved', 'Error', 'modal');
        return 
 end

if (plot_init == 2) %'Horizons depths'
       VAR=depth;
       IMAX=size(depth,2);
       string_run = strcat('_initial_compacted_depths_',runid,ext);

elseif (plot_init == 3) %,'EET'
       VAR=EET_init;
       IMAX=1;
       string_run = strcat('Initial_EET_',runid,ext);

elseif (plot_init == 4) %,'SL'
       VAR=SL_init;
       IMAX=size(SL_init,2);
       string_run = strcat('_initial_Sea_Level_',runid,ext);

elseif (plot_init == 5) %,'BETA FACTOR'
       VAR=betaF';
       IMAX=1;
       string_run = strcat('Initial_Beta_factor_',runid,ext);
 elseif (plot_init == 6) %,'Dynamic topography'
        VAR=DynTopo_init;
        IMAX=size(DynTopo_init,2);
        string_run = strcat('_initial_Dynamic_Topography_',runid,ext);
end 

          
% SAVE VARIABLE TO FILE
if (plot_init == 2 || plot_init == 4 || plot_init == 6)
          for ii=1:IMAX % for each time
            if (plot_init == 2)
                filename=[profiles(IMAX-ii+1).layername string_run];
            else
                filename=[profiles(ii+1).layername string_run];
            end

            
            fpath=fullfile(destinationfolder,filename);
            
            [fid]=check_save_to_file(fpath,destinationfolder,filename);
            
            to_be_written=VAR(:,ii);

            %CHOOSE GRID DIMENSION
             if (opt_grid == 1)
                 fprintf(fid,[profiles(ii).layername '  and deeper \n']); % labels
                 formatstring=['%4.4f' '\n'];
                 fprintf(fid,formatstring,[to_be_written(1,:)]');
                 fclose(fid);

             elseif (opt_grid == 2)
                 fprintf(fid,[' X (m) ' profiles(ii).layername '  and deeper \n']); % labels
                 formatstring=['%6.6f %4.4f' '\n'];
                 fprintf(fid,formatstring,[M(:,1), to_be_written]');
                 fclose(fid);

             elseif(opt_grid == 3)
                 fprintf(fid,[' X (m) Y(m) ' profiles(ii).layername '  and deeper \n']); % labels
                 %formatstring=['%6.6f %6.6f' repmat(' %4.4f',1,size(to_be_written,2)) '\n'];
                 formatstring=['%6.6f %6.6f %4.4f' '\n'];
                 fprintf(fid,formatstring,[M(:,1), M(:,2), to_be_written]');
                 fclose(fid);
             end
          end
elseif (plot_init == 3) %'Elastic thickness'        
             to_be_written=VAR(:);
             filename=string_run; 
             fpath=fullfile(destinationfolder,filename);
            
             [fid]=check_save_to_file(fpath,destinationfolder,filename);         

            %CHOOSE GRID DIMENSION
             if (opt_grid == 1)
                 fprintf(fid,'EET (km) \n'); % labels
                 formatstring=['%4.4f' '\n'];
                 fprintf(fid,formatstring,[to_be_written(1,:)]');
                 fclose(fid);

             elseif (opt_grid == 2)
                 fprintf(fid,[' X (m)  EET (km) \n']); % labels
                 formatstring=['%6.6f %4.4f' '\n'];
                 fprintf(fid,formatstring,[M(:,1), to_be_written]');
                 fclose(fid);

             elseif(opt_grid == 3)
                 fprintf(fid,[' X (m)  Y(m)  EET(km) \n']); % labels
                 formatstring=['%6.6f %6.6f %4.4f' '\n'];
                 fprintf(fid,formatstring,[M(:,1), M(:,2), to_be_written]');
                 fclose(fid);
             end
             
elseif (plot_init == 5) %'Beta Factor'        
             to_be_written=VAR(:);
             filename=string_run; 
             fpath=fullfile(destinationfolder,filename);
            
             [fid]=check_save_to_file(fpath,destinationfolder,filename);         

            %CHOOSE GRID DIMENSION
             if (opt_grid == 1)
                 fprintf(fid,'Beta \n'); % labels
                 formatstring=['%4.4f' '\n'];
                 fprintf(fid,formatstring,[to_be_written(1,:)]');
                 fclose(fid);

             elseif (opt_grid == 2)
                 fprintf(fid,[' X (m)  Beta \n']); % labels
                 formatstring=['%6.6f %4.4f' '\n'];
                 fprintf(fid,formatstring,[M(:,1), to_be_written]');
                 fclose(fid);

             elseif(opt_grid == 3)
                 fprintf(fid,[' X (m)  Y(m)  Beta \n']); % labels
                 formatstring=['%6.6f %6.6f %4.4f' '\n'];
                 fprintf(fid,formatstring,[M(:,1), M(:,2), to_be_written]');
                 fclose(fid);
             end
               
 end
