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
% Plot initial input data selected from the "Plot GUI" and related options:
%
% - plot min and max automatically calculated or prescribed from the GUI
% - for 3D backtracking: plot with plain map view or real 3D visual effect
%
%-----------------------------------------------------------------------------------
function Plot_Init(M,grid_opt,plot_init,maps_plot,plot_min,plot_max,plot_levels,profiles,depth,EET_init,SL_init,DynTopo_init,betaF)



% Get coordinates
 X=M(:,1);
 Y=M(:,2);
 
 iscollinear=rank([X-X(1), Y-Y(1)])<2; %true if the points are collinear or there is a single point
 
 if (plot_init == 1)
        errordlg('Select a variable to be plotted', 'Error', 'modal');
        return 
 end

if (plot_init == 2) %'Initial compacted present-day horizons depths'
       IMAX=size(depth,2);
       ZMAX = max(-depth,[],'all');
       ZMIN = min(-depth,[],'all');
elseif (plot_init == 3) %,'Input Effective Elastic Thickness'
       IMAX=1;
       ZMAX = max(EET_init,[],'all');
       ZMIN = min(EET_init,[],'all');
elseif (plot_init == 4) %,'Input Sea level change at sediment layer age'
       IMAX=size(SL_init,2);
       ZMAX = max(SL_init,[],'all');
       ZMIN = min(SL_init,[],'all');
elseif (plot_init == 5) %,'Input Sea level change at sediment layer age'
       IMAX=1;
       ZMAX = max(betaF,[],'all');
       ZMIN = min(betaF,[],'all');
elseif (plot_init == 6) %,'Input dynamic topography at sediment layer age'
       IMAX=size(DynTopo_init,2);
       ZMAX = max(DynTopo_init,[],'all');
       ZMIN = min(DynTopo_init,[],'all');
end   


%Select plot min and max as prescribed by users in the GUI:
% Overwrites calculated ZMIN and ZMAX
if (plot_levels == 1)
    ZMIN = plot_min;
    ZMAX = plot_max;
end

%INITIALISE FIGURE ARRAY
 vhnd=[];
 try
     warning off
 if ~iscollinear
     tri=delaunay(X,Y);
 
    
%--------------------------------------------------------------------------
% 3D plots
%-------------------------------------------------------------------------- 
for ii=1:IMAX
 hnd=figure; vhnd=[vhnd hnd];

%---------------------------------------
% HORIZONS DEPTH
%---------------------------------------      
    if (plot_init == 2) 
     P=size(depth,2);

     to_be_plotted=-depth(:,ii);
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
     caxis([ZMIN ZMAX]); %-4000 --> ZMIN
     cb=colorbar;
     cb.Label.String = 'Depth (m)';

     % Contours on plot
     clevels=round(ZMIN):round(abs(ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');
     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end             

     title ([num2str(profiles(IMAX-ii+1).age) 'Ma -' profiles(IMAX-ii+1).layername], 'Interpreter', 'none')

%---------------------------------------
% SEA LEVEL CHANGE
%--------------------------------------- 
  elseif (plot_init == 4) %,'SL_init'

     P=size(SL_init,2);
     to_be_plotted=SL_init(:,ii);

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

     colormap parula;

     cb=colorbar;
     cb.Label.String = 'Sea Level correction (m)';
     if (ZMIN ~= ZMAX)
        caxis([ZMIN ZMAX]); % In case EET is not constant
     end
     % Contours on plot
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end  

     title (upper([profiles(ii+1).layername]), 'Interpreter', 'none')

%---------------------------------------
% DYNAMIC TOPOGRAPHY
%--------------------------------------- 
  elseif (plot_init == 6) %,'Dynamic topography'
 
    P=size(DynTopo_init,2);
     to_be_plotted=DynTopo_init(:,ii);

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

     colormap parula;

     cb=colorbar;
     cb.Label.String = 'Dyn. Topo. correction (m)';
     if (ZMIN ~= ZMAX)
        caxis([ZMIN ZMAX]); % In case EET is not constant
     end
     % Contours on plot;
     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black');

     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end  

     title (upper([profiles(ii+1).layername]), 'Interpreter', 'none')

    end 
end
%---------------------------------------
% EFFECTIVE ELASTIC THICKNESS
%---------------------------------------    
 if (plot_init == 3) %,'EET'
     to_be_plotted=EET_init;

     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);

     figure
     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)'); hold on
     ylabel('Y cart. (m)'); hold on
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');
     
     if(maps_plot == 1) % 3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end
     
     colormap jet;
     cb=colorbar; hold on
     cb.Label.String = 'EET (m)';
     
     if (ZMIN ~= ZMAX)
        caxis([ZMIN ZMAX]); % In case EET is not constant
     end

     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black'); hold on

     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end 

     title ('Effective Elastic Thickness', 'Interpreter', 'none')

%---------------------------------------
% BETA FACTOR
%--------------------------------------- 
 elseif (plot_init == 5) % Beta factor

     to_be_plotted=betaF;

     [var3D,X2d,Y2d]=mesh3D_plot(X,Y,to_be_plotted);
     
     figure
     h=mesh(X2d,Y2d,var3D); hold on
     xlabel('X cart. (m)'); hold on
     ylabel('Y cart. (m)'); hold on
     set(h,'edgealpha',0,'facecolor','interp','facelighting','gouraud');

     if(maps_plot == 1) % 3D maps
       camlight(-61,43);
     else               % 2D maps
       material dull;
     end

     colormap jet;
     cb=colorbar; hold on
     cb.Label.String = 'Beta factor (-)';
     if (ZMIN ~= ZMAX)
        caxis([ZMIN ZMAX]); % In case Beta Factor is not constant
     end

     clevels=round(ZMIN):round((ZMAX-ZMIN)/30):round(ZMAX);
     [C,h]=contour3(X2d,Y2d,var3D,clevels,'black'); hold on
     
     % 2D flat maps
     if(maps_plot == 0) 
       view(2);
     end 

     title ('Beta factor', 'Interpreter', 'none')
     
 end %END IF plot_init        
        

         
%------------------------------------------------------------------------
%  PLOT: 2D profiles
%------------------------------------------------------------------------
   
else % iscollinear
    cm=colormap(parula(IMAX));
  
  figure
  for ii=1:IMAX
    %hnd=figure; vhnd=[vhnd hnd];

%---------------------------------------
% HORIZONS DEPTH
%--------------------------------------- 
    if (plot_init == 2)
        
           if (grid_opt == 2)
            h=plot(X, -depth(:,ii)); hold on
           elseif (grid_opt == 1)
            h=plot(-depth(50:60,ii)); hold on
           end 
            
            % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);
            
            cmap = parula(IMAX);
            RGB1 = cmap(ii,:);
            % Set profiles lines properties            
            if (ii == 1)
              set(h,'LineWidth',3,'Color',RGB1);
            else
              set(h,'LineWidth',2,'Color',RGB1);
            end
            
            % Create legend
            for i= 1:IMAX
               leg(i)=arrayfun(@num2str,profiles(IMAX-i+1).age,'UniformOutput',false);
            end
            
            for i=1:IMAX
               RGB2 = cmap(i,:);
               hh(i) = plot(NaN,'Color', RGB2);
            end
            legend(hh,leg,'Location','EastOutside');
            [hleg,att] = legend('show');
            title(hleg,'Age (Ma)')
            
         title (upper('Present-day - compacted depths (m)'), 'Interpreter', 'none')
%---------------------------------------
% SEA LEVEL CHANGE
%---------------------------------------          
    elseif (plot_init == 4)

           if (grid_opt == 2)
            h=plot(X,SL_init(:,ii)); hold on
           elseif (grid_opt == 1)
            h=plot(SL_init(50:60,ii)); hold on
           end 

            % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);
            
            cmap = parula(IMAX);
            RGB1 = cmap(ii,:);
            if (ii == 1)
              set(h,'LineWidth',3,'Color',RGB1);
            else
              set(h,'LineWidth',3,'Color',RGB1);
            end
            
            % Create legend
            for i= 1:IMAX
               leg(i)=arrayfun(@num2str,profiles(IMAX-i+2).age,'UniformOutput',false);
            end
            
            for i=1:IMAX
               RGB2 = cmap(i,:);
               hh(i) = plot(NaN,'Color', RGB2);
            end
            legend(hh,leg,'Location','EastOutside');
            [hleg,att] = legend('show');
            title(hleg,'Age (Ma)')
            
            
         title (upper('Sea Level Correction (meters)'), 'Interpreter', 'none')

%---------------------------------------
% DYNAMIC TOPOGRAPHY
%--------------------------------------- 
    elseif (plot_init == 6)
        
           if (grid_opt == 2)
            h=plot(X,DynTopo_init(:,ii)); hold on
           elseif (grid_opt == 1)
            h=plot(DynTopo_init(50:60,ii)); hold on
           end         

            % Set axis limits
            axis([min(X) max(X) ZMIN ZMAX]);
            
            cmap = parula(IMAX);
            RGB1 = cmap(ii,:);
            if (ii == 1)
              set(h,'LineWidth',3,'Color',RGB1);
            else
              set(h,'LineWidth',3,'Color',RGB1);
            end
            
            % Create legend
            for i= 1:IMAX
               leg(i)=arrayfun(@num2str,profiles(IMAX-i+2).age,'UniformOutput',false);
            end
            
            for i=1:IMAX
               RGB2 = cmap(i,:);
               hh(i) = plot(NaN,'Color', RGB2);
            end
            legend(hh,leg,'Location','EastOutside');
            [hleg,att] = legend('show');
            title(hleg,'Age (Ma)')
            
            
         title (upper('Dynamic Topography Correction (meters)'), 'Interpreter', 'none')
     end
  end

%---------------------------------------
% EFFECTIVE ELASTIC THICKNESS
%---------------------------------------         
if (plot_init == 3)
    
       if (grid_opt == 2)
        h=plot(X,EET_init); hold on
       elseif (grid_opt == 1)
        h=plot(EET_init(50:60)); hold on
       end 

        % Set axis limits
        if (ZMIN ~= ZMAX) % In the EET is not constant
           axis([min(X) max(X) ZMIN ZMAX]);
        end
        
        set(h,'LineWidth',3,'Color','black');

        title (upper('Effective Elastic Thickness (km)'), 'Interpreter', 'none')

%---------------------------------------
% BETA FACTOR
%---------------------------------------         
elseif (plot_init == 5)
    if (grid_opt == 2)
    h=plot(X,betaF); hold on
    elseif (grid_opt == 1)
    h=plot(betaF(50:60)); hold on
    end 
       

    % Set axis limits
    if (ZMIN ~= ZMAX) % In the case Beta Factor is not constant
       axis([min(X) max(X) ZMIN ZMAX]);
    end

    set(h,'LineWidth',3,'Color','black');

    title (upper('Beta Factor'), 'Interpreter', 'none')

end
 
 
 end  % PLOT 2D 
 warning on;
 catch
     close(vhnd)
     errordlg('The data cannot be visualized. Possible causes: NaNs and/or too large values due to wrong parameters.','', 'modal')
     warning on;     
 end


end