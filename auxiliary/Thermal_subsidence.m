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
% THERMAL SUBSIDENCE
%
% Based on the 1D model of McKenzie et al (1978)
% Calculate syn-rift and post-rift subsidence
% Allows only for one phase of rifting (implementation of several riftings is in progress)
%
%
% Different ways to prescribe the streching factor (beta factor: Sfactor) from the GUI:
%
% [1] Spatially uniform Beta factor: only B1 value is used
%
% [2] Linear spatial interpolation following X direction: interpolate between B1 and B2
%
% [3] Step-wise interpolation: uses B1 if X < Xlim and B2 otherwise
%     Xlim is defined with the 'coeff' var form the GUI and then as follows:
%     - if 'coeff' is set at 0.5: it means that the grid domain limit is at
%     half domain based on Xmin and Xmax coordinates
%
%     - if 'coeff' is set at 0.33: it means that the grid domain limits is at
%     one third based on Xmin and Xmax coordinates
% and so on
%
%
% [4] User-provided file of 2d maps of spatially variable Beta factor 
%------------------------------------------------------------------------------------------- 

function [ts_corr,Sfactor]=Thermal_subsidence(NN,Np,S,Xcoord,Ycoord,opt_subs,opt_BF,user_beta,Blimit,B1,B2,...
         filebeta,dirbeta,tlitho,tmoho,alpha,kappa,rift,Rhom,Rhoc,Rhow,ages)

     
%----------------------------------------------------       
%--- DEFINE STRECHING FACTOR BASED ON GUI OPTIONS
%----------------------------------------------------       
ts_corr=zeros(Np,NN);   % Total subsidence correction relative to present            

if (opt_subs == 1) % THERMAL SUBSIDENCE ON
    
     if (opt_BF == 1) % Spatially uniform (prescribed from interface)
         for k=1:Np
           Sfactor(k) = B1;
         end

     % opt_BF == 2 correponds to the first choice in the
     % menu on the GUI: "Beta factor Interpolation" and retrun the interface ERROR message:
     % "Please select an interpolation methods"

     elseif (opt_BF == 3) % Step_wise function (with limit set at 'coeff' from GUI, and using B1,B2 prescribed from GUI)

         Xlim = min(Xcoord)+(max(Xcoord) - min(Xcoord))*Blimit;
         %Ylim = min(Ycoord)+(max(Ycoord) - min(Ycoord))*Blimit;

         for k=1:Np
             if (Xcoord(k) < Xlim)
                 Sfactor(k)=B1;
             else
                 Sfactor(k)=B2; 
             end
         end

     elseif (opt_BF == 4) % Linear interpolation(East to West, two Beta values prescribed from GUI)
         for k=1:Np  
           Sfactor(k) = B1 + (B2-B1)*(Xcoord(k)-Xcoord(1))/(Xcoord(Np)-Xcoord(1));
         end

     end        

     if (user_beta == 1) % User-provided 2D maps of beta factors   
         
         [Beta]=read_beta(filebeta,dirbeta,Xcoord);
         Sfactor(:) = Beta(:) ;
     end

%----------------------------------------------------                    
%--- COMPUTE THERMAL SUBSIDENCE AT HORIZONS AGE
%----------------------------------------------------       

%--- SYN-RIFT SUBSIDENCE
     tlitho = tlitho*1000. ;  % Initial thickness of lithosphere (km) Busetti et al.
     tmoho = tmoho*1000. ;    % Initial thickness of the crust (km) Decesari et al.
     PI = 3.14159265359 ;
     ma2s = 1e6*3600*24*365 ; % convert Myrs into seconds
     Tm = 1330 ;              % temperature of the asthenosphere (°C)
     N=NN-1;
     
     
     for k=1:Np
       Ys(k) = tlitho* ((Rhom - Rhoc) * (tmoho/tlitho)*( 1 - (alpha*Tm*tmoho)/(2*tlitho) ) - (alpha*Tm*Rhom)/2 ) *... 
           (1- 1/Sfactor(k))/ ( Rhom * (1 - alpha*Tm) - Rhow);
     end
                 
%--- POST-RIFT THERMAL SUBSIDENCE

     E0=zeros(NN);
     S=zeros(Np,NN);         % Post_rift thermal subsidence
     tsubs=zeros(Np,NN);     % Total subsidence

     % Retrieve Ages (Ma) of each horizon, including present-day
     hzages=[ages{:}];       % Age of the different layers from input parameters file
     thzages=[hzages 0]';    % Transpose ages for calculation and add present-day (0 Ma)
     
         for t=1:NN % NN = N+1 because it includes present-day
           for k=1:Np
             E0(t) = 4*tlitho * Rhom * alpha * Tm/(PI.^2 * (Rhom - Rhow));
             tau = (tlitho).^2 / (PI.^2*kappa) ;

             S(k,t) = E0(t)* (Sfactor(k)/PI)*sin(PI/Sfactor(k)).* (1 - exp( -(rift-thzages(t))*ma2s/tau) ); %(rift-ages(i)) --> time ellapsed since rifting

             % TOTAL SUBSIDENCE = SYN-RIFT SUBS + POST-RIFT SUBS
             tsubs(k,t) = Ys(k) + S(k,t);
           end
         end


%--- CALCULATION OF THERMAL SUBSIDENCE CORRECTION
% Remode present-day value and add interpolated subsidence at horizon ages
  for t=1:N
         
     if (thzages(t) < rift) % Post rift
       for k=1:Np
          ts_corr(k,t) = tsubs(k,N-t+1)-tsubs(k,NN);
       end
       
     else % If basement older than end of rift
         'basement is older than rift ages'
        for k=1:Np
          ts_corr(k,t) = Ys(k)-tsubs(k,NN);
       end
     end
  end


else %THERMAL SUBSIDENCE OFF
     for t=1:NN
            ts_corr(:,t) = 0.;
            Sfactor(:) = 0.;
     end

end %END OPT

end