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
% Plot extracted transects (extended horizontal coordinates in case of 1D well, for 
% visualization) selected from the "Plot GUI"
%
%-----------------------------------------------------------------------------------
function Plot_Extracted(M,plot_extracted,extract_prof,Xpts,profiles) %#ok

            

if (plot_extracted < 6 || plot_extracted == 7 )
       IMAX=length(extract_prof);
else % BETA FACTOR AND THERMAL SUBSIDENCE
       IMAX=size(extract_prof,2);
       IMAX=size(extract_prof,2);
end


vhnd=[];
cm=colormap(parula(IMAX));

if (plot_extracted < 6 || plot_extracted == 7 )
 for ii=1:IMAX
     hnd=figure; vhnd=[vhnd hnd];
     
     % Plot extracted
     for jj=size(extract_prof(ii).interface,2):-1:1
         to_be_plotted = extract_prof(ii).interface(2:end-1,jj);

         h=plot(Xpts(2:end-1),to_be_plotted); hold on
         set(h,'LineWidth',2,'Color',cm(mod(jj-1,length(cm))+1,:));
         xlabel('X cart. (m)');
     end
     
    % Get plot title 
    if (plot_extracted == 2) %'Horizons depths'
     title (['Horizons depths - ' num2str(profiles(ii).age) ' Ma   ' upper(profiles(ii).layername) '      ' '(layers removed: ' num2str(ii-1) ')'   ], 'Interpreter', 'none')
    elseif (plot_extracted == 3) %,'Isopacks'
      title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')
    elseif (plot_extracted == 4) %,'Density'
      title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')
    elseif (plot_extracted == 5) %,'Porosity'
      title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')
    elseif (plot_extracted == 7) %,'Isostatic correction'
      title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')
    end

    grid on
    hold off ;
 end

else  % For Thermal subsidence and Beta Factor

  for ii=1:IMAX
     hnd=figure; vhnd=[vhnd hnd];

     % Plot variables
     to_be_plotted = extract_prof(:,ii);
     h=plot(Xpts,to_be_plotted); hold on

     if(size(extract_prof,2) > 1)
        set(h,'LineWidth',2,'Color',cm(mod(ii-1,length(cm))+1,:));
     else
        set(h,'LineWidth',2);
     end

     xlabel('X cart. (m)');

    % Get titles
    if (plot_extracted == 6) %,'Thermal Subsidence'   
     title (upper([profiles(ii+1).layername '      to     ' profiles(ii).layername]), 'Interpreter', 'none')    

    elseif (plot_extracted == 8) %,'Beta Factor'
     title ('Beta factor', 'Interpreter', 'none')
    end

    grid on
    hold off;
  end  % END FOR LOOP  

end % END SELECTED VARIABLES


end