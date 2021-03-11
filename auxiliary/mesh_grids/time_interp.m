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
% Interpolate prescribed Ages through the GUI for input dynamic topography time slices 
% on input sediment horizon ages
%
%-----------------------------------------------------------------------------------
function [dt_DT]=time_interp(NT,N,dynAges,tHages,DT_corr)

% Check time order
if ( dynAges(2) > dynAges(1) )
  ll_dynAges_backward = 1;
else
  ll_dynAges_backward = 0;
end

%

% Interpolate DT on tHages array

  for t=N:-1:1
      %Interpolate on tHages  
        if ( ll_dynAges_backward == 1)
          kinf=1;
          for k=2:NT-1
            if ( dynAges(k) <= tHages(t) ) 
                kinf = k;
            end
                       
          end
          ksup = kinf + 1;
        else
          ksup=NT;
          for k=NT-1:-1:2
            if ( dynAges(k) >= tHages(t) ) 
                ksup = k;
            end
          end
          kinf = ksup - 1;
        end
        
        
       dt_DT(:,t) =  (DT_corr(:,kinf).*(dynAges(ksup)-tHages(t))+ (tHages(t)-dynAges(kinf)).* DT_corr(:,ksup))./(dynAges(ksup)-dynAges(kinf));            
        

  end
    