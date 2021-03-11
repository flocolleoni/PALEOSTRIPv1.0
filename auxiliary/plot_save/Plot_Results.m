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
% Plot computed variables selected from the "Plot GUI" and related options:
%
% - plot min and max automatically calculated or prescribed from the GUI
% - for 3D backtracking: plot with plain map view or real 3D visual effect
%
%-----------------------------------------------------------------------------------

function Plot_Results(M,grid_opt,plot_menu,maps_plot,plot_min,plot_max,plot_levels,profiles,profiles_min,profiles_max,isopachs,isopachs_min,isopachs_max,...
                density,density_min,density_max,porosity,porosity_min,porosity_max,tsubs_corr,tsubs_min,...
                tsubs_max,isostasy,isostasy_min,isostasy_max,tectonic,tectonic_min,tectonic_max)


% Get input grid coordinates
 X=M(:,1);
 Y=M(:,2);
 
 iscollinear=rank([X-X(1), Y-Y(1)])<2; %true if the points are collinear or there is a single point
 
 if (plot_menu == 1)
        errordlg('Select a variable to be plotted', 'Error', 'modal');
        return 
 end

% Get variables dimensions
if (plot_menu == 2) %'Horizons depths'
       IMAX=length(profiles);
       ZMAX = profiles_max;
       ZMIN = profiles_min;
elseif (plot_menu == 3) %,'Isopachs'
       IMAX=length(isopachs); 
       ZMAX = isopachs_max;
       ZMIN = isopachs_min;
elseif (plot_menu == 4) %,'Density'
       IMAX=length(density) ;
       ZMAX = density_max;
       ZMIN = density_min;
elseif (plot_menu == 5) %,'Porosity'
       IMAX=length(porosity) ;
       ZMAX = porosity_max;
       ZMIN = porosity_min;
elseif (plot_menu == 6) %,'Thermal Subsidence'
       IMAX=size(tsubs_corr,2);
       ZMAX = tsubs_max;
       ZMIN = tsubs_min;  
elseif (plot_menu == 7) %,'Isostatic correction'
       IMAX=length(isostasy);
       ZMAX = isostasy_max;
       ZMIN = isostasy_min;
elseif (plot_menu == 8) % Tectonic subsidence of basement
       IMAX=size(tectonic,2);
       ZMAX = tectonic_max;
       ZMIN = tectonic_min  ;    
end   


%Select plot min and max as prescribed by users in the GUI:
% Overwrites calculated ZMIN and ZMAX
if plot_levels == 1
    ZMIN = plot_min;
    ZMAX = plot_max;
end

%INITIALISE FIGURES ARRAY
 vhnd=[];
 try
     warning off
 if ~iscollinear
     tri=delaunay(X,Y);
 

%--------------------------------------------------------------------------
% 3D plots (plain maps or with 3D visual)
%-------------------------------------------------------------------------- 
for ii=1:IMAX
 hnd=figure; vhnd=[vhnd hnd];

%---------------------------------------
% HORIZONS DEPTH
%---------------------------------------      
    if (plot_menu == 2) 
     P=size(profiles(ii).interface,2);

     to_be_plotted=profiles(ii).interface(:,P);
     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');
     if(maps_plot == 0) % Flat 2D maps
       material dull;
     else               % Real 3D maps
       camlight(-61,43);
     end

     colormap(flipud(hsv));
     caxis([ZMIN ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Elevation (m)';

     % Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');
     
     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end             

     title ([num2str(profiles(ii).age) ' Ma   ' upper(profiles(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')

%---------------------------------------
% ISOPACHS
%--------------------------------------- 
  elseif (plot_menu == 3) %,'isopachs'

     P=size(isopachs(ii).interface,2);
     to_be_plotted=isopachs(ii).interface(:,P);

     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');

     if(maps_plot == 1) % Real 3D maps
       camlight(-61,43);
     else               % 2D flat maps
       material dull;
     end

     colormap parula;
     caxis([ZMIN ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Sediment thickness (m)';

     % Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end  

     title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')

%---------------------------------------
% DENSITY
%--------------------------------------- 
 elseif (plot_menu == 4) %,'Density'

     P=size(density(ii).interface,2);
     to_be_plotted=density(ii).interface(:,P);

     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');
             
     if(maps_plot == 1) % 3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end

     colormap jet;
     caxis([1800 ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Density (kg/m3)';

     % Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D maps flat view
     if(maps_plot == 0) 
       view(2);
     end  

     title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')

%---------------------------------------
% POROSITY
%--------------------------------------- 
 elseif (plot_menu == 5) %,'Porosity'

     P=size(porosity(ii).interface,2);
     to_be_plotted=porosity(ii).interface(:,P);

     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');

     if(maps_plot == 1) % 3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end

     colormap jet;
     caxis([ZMIN ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Porosity (%)';

     % Contours on plot
      clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D maps flat view
     if(maps_plot == 0) 
       view(2);
     end 

     title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')
         
%---------------------------------------
% ISOSTASY
%---------------------------------------        
 elseif (plot_menu == 7) %,'Isostasy'

     P=size(isostasy(ii).interface,2);
     to_be_plotted=isostasy(ii).interface(:,P);

     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);
     %clevels=[ZMIN:25:ZMAX];

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');

     if(maps_plot == 1) %3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end

     colormap jet;
     caxis([ZMIN ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Isostasy (m)';

     % Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end 

     title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')
   
              
     
%---------------------------------------
% THERMAL SUBSIDENCE
%---------------------------------------  
elseif (plot_menu == 6) %,'Thermal Subsidence'

     to_be_plotted=tsubs_corr(:,ii);
     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');

     if(maps_plot == 1) % 3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end
     
     colormap jet;
     caxis([ZMIN ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Thermal subsidence correction (m)';

     % Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end 

     title (upper([profiles(ii).layername]), 'Interpreter', 'none')

%---------------------------------------
% TECTONIC SUBSIDENCE HISTORY
%---------------------------------------    
elseif (plot_menu == 8) %,'Tectonic Subsidence of Basement' 
     
     to_be_plotted=tectonic(:,ii);
     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)');
     ylabel('Y cart. (m)');
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');
     
     if(maps_plot == 1) % 3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end
     
     colormap jet;
     caxis([ZMIN ZMAX]);
     cb=colorbar;
     cb.Label.String = 'Tectonic subsidence (m)';

     %Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');
     
     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end 

     title (upper([profiles(ii).layername]), 'Interpreter', 'none')

 end %END IF PLOT_MENU        
end

         
%------------------------------------------------------------------------
%  PLOT: 2D profiles
%------------------------------------------------------------------------
   
else % iscollinear
    cm=colormap(parula(IMAX));
    
  for ii=1:IMAX
    hnd=figure; vhnd=[vhnd hnd];

%---------------------------------------
% HORIZONS DEPTH
%--------------------------------------- 
    if (plot_menu == 2)
         for jj=size(profiles(ii).interface,2):-1:1

           if (grid_opt == 2)
             h=plot(X, profiles(ii).interface(:,jj)); hold on
           elseif (grid_opt == 1)
             h=plot(profiles(ii).interface(50:60,jj)); hold on
           end

            % Set axis limits
              axis([min(X) max(X) ZMIN ZMAX]);
            
            % Set profiles lines properties            
            if (jj == 1)
              set(h,'LineWidth',3,'Color','black');
            else
              set(h,'LineWidth',2,'Color',cm(mod(jj-1,length(cm))+1,:));
            end
         end

         title ([num2str(profiles(ii).age) ' Ma   ' upper(profiles(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')

%---------------------------------------
% ISOPACH
%---------------------------------------          
    elseif (plot_menu == 3)
         for jj=size(isopachs(ii).interface,2):-1:1

           if (grid_opt == 2)
            h=plot(X,isopachs(ii).interface(:,jj)); hold on
           elseif (grid_opt == 1)
            h=plot(isopachs(ii).interface(50:60,jj)); hold on
           end             
             

            % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);

            if (jj == 1)
              set(h,'LineWidth',3,'Color','black');
            else
              set(h,'LineWidth',2,'Color',cm(mod(jj-1,length(cm))+1,:));
            end
         end

         title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')

%---------------------------------------
% DENSITY
%--------------------------------------- 
    elseif (plot_menu == 4)
         for jj=size(density(ii).interface,2):-1:1

           if (grid_opt == 2)
            h=plot(X,density(ii).interface(:,jj)); hold on
           elseif (grid_opt == 1)
            h=plot(density(ii).interface(50:60,jj)); hold on
           end                

            % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);

            if (jj == 1)
              set(h,'LineWidth',3,'Color','black');
            else
              set(h,'LineWidth',2,'Color',cm(mod(jj-1,length(cm))+1,:));
            end
         end

         title ([num2str(density(ii).age) ' Ma   ' upper(density(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')

%---------------------------------------
% POROSITY
%--------------------------------------- 
    elseif (plot_menu == 5)
         for jj=size(porosity(ii).interface,2):-1:1

           if (grid_opt == 2)
            h=plot(X,porosity(ii).interface(:,jj)); hold on
           elseif (grid_opt == 1)
            h=plot(porosity(ii).interface(50:60,jj)); hold on
           end 

            % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);

            if (jj == 1)
              set(h,'LineWidth',3,'Color','black');
            else
              set(h,'LineWidth',2,'Color',cm(mod(jj-1,length(cm))+1,:));
            end
         end

         title ([num2str(porosity(ii).age) ' Ma   ' upper(porosity(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')

%---------------------------------------
% THERMAL SUBSIDENCE
%--------------------------------------- 
    elseif (plot_menu == 6)

           if (grid_opt == 2)
               h=plot(X,tsubs_corr(:,ii)); hold on
           elseif (grid_opt == 1)
               h=plot(tsubs_corr(50:60,ii)); hold on
           end         

        % Set axis limits
        axis([min(X) max(X) ZMIN ZMAX]);

        if (ii == 1)
          set(h,'LineWidth',3,'Color','black');
        else
          set(h,'LineWidth',2,'Color',cm(mod(ii-1,length(cm))+1,:));
        end

        title ([num2str(profiles(ii+1).age) ' Ma   ' upper(profiles(ii+1).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')

%---------------------------------------
% ISOSTASY
%---------------------------------------          
    elseif (plot_menu == 7)
        for jj=size(isostasy(ii).interface,2):-1:1
            
           if (grid_opt == 2)
              h=plot(X,isostasy(ii).interface(:,jj)); hold on
           elseif (grid_opt == 1)
              h=plot(isostasy(ii).interface(50:60,jj)); hold on
           end  

           % Set axis limits
           axis([min(X) max(X) ZMIN ZMAX]);

           if (jj == 1)
             set(h,'LineWidth',3,'Color','black');
           else
             set(h,'LineWidth',2,'Color',cm(mod(jj-1,length(cm))+1,:));
           end
        end

        title ([num2str(isostasy(ii).age) ' Ma   ' upper(isostasy(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none') 


%---------------------------------------
% Tectonic subsidence
%---------------------------------------
    elseif (plot_menu == 8) 
        for ii=1:IMAX-ii+1
            
           if (grid_opt == 2)
              h=plot(X,tectonic(:,ii)); hold on
           elseif (grid_opt == 1)
              h=plot(tectonic(50:60,ii)); hold on
           end  
           
            
             % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);
              
            if (ii == 1)
              set(h,'LineWidth',3,'Color','black');
            else
              set(h,'LineWidth',2,'Color',cm(mod(ii-1,length(cm))+1,:));
            end
        end % end for loop
        
        title ([num2str(profiles(ii).age) ' Ma   ' upper(profiles(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')

    end % END PLOT MENU OPTION    
    
     grid on
 end 

 
 end  % PLOT 2D 
 warning on;
 catch
     close(vhnd)
     errordlg('The data cannot be visualized. Possible causes: NaNs and/or too large values due to wrong parameters.','', 'modal')
     warning on;     
 end


end