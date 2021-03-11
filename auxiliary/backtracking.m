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
%
% Backtracking is an implementation in Matlab language of the BASIC program
% reported at the end of chapter 8 of "Principles of Sedimentary Basin
% Analysis"
% [ Parallel version for faster computation]
% 
%   
%     OUTPUTS:  
%     out.Nzt contains the decompacted depth of the top surface of the
%     layers in different times
%
%    out.Nzb         % backstripped depths (in meters)
%    out.Dt          % decompacted thicknesses (in meters)
%    out.T           % Initial thicknesses
%    out.Rhoa        % Averaged density
%    out.Rhoa_units  % Density of each units
%    out.Iso_units   % Isostasy of each units
%    out.Por         % Porosity of each units
%    out.tsubs_corr  % Thermal subsidence
%    out.beta_factor % Streching factor
%    out.Tect        % Tectonic Subsidence
%    
%    % INITIAL DATA AND CONDITIONS
%    out.EET         % Elastic thickness
%    out.SL          % Sea level change
%    out.dyntopo     % Initial dynamic topography
%--------------------------------------------------------------------------         

function [out, outputstring]=backtracking(S,params,grid_opt,iso_opt,EET_opt,EETfile,EETdir,EET,YM,v,opt_sl,dS,...
    SLfile,SLdir,opt_subs,rift,alpha,user_beta,filebeta,dirbeta,opt_BF,B1,B2,Blimit,kappa,tlitho,tmoho,...
    opt_dyn,dtopo_cte,dtopo_file,dtopo_dir,Rhom,Rhoc,Rhow)

    % DATA FOLDER: contains sea level curves and other data needed for
    % backstripping


    %Get coordinates
    Xcoord=S(:,1);
    Ycoord=S(:,2);
    
    % convert depth from m to km
    S(:,3:end)=S(:,3:end)/1000.;
    

    N=length(params); % Number of stratigraphic units
    NN=N+1; % Number of stratigraphic horizon (includes bathymetry)
    Surpor = [params.Surpor]; % Estimated surface porosity
    C = [params.C]; % Porosity-depth coefficient
    Rhos = [params.Rhos]; % Sediment grain density in kg / m^3
    %L = {params.Name}; % Name of stratigraphic unit
    ages = {params.AgeBase};
    
    if NN~=(size(S,2)-2) %N~=(size(S,2)-2)
        out=[];
        outputstring='Number of stratigraphic units not consistent with the data.';
        return
    end
       
    Np = size(S,1); %number of points


%--------------------------------------------
% Get bottom (Zb) and top (Zt) of each layer
%--------------------------------------------     
    
    Zb=zeros(Np,NN); % bottom
    Zt=zeros(Np,NN); % top
    T=zeros(Np,N);  % thickness
    ZZb=zeros(Np,N); % Variable to update thickness for each time step of main loop
    ZZt=zeros(Np,N); % Variable to update thickness for each time step of main loop

    
%--- ATTRIBUTE TOP AND BOTTOM OF PRESENT_DAY LAYERS
%--- from bottom to top j=1 ---> basement ; j=N+1 ---> bathymetry (Zt(N))
%--- !! Depth of each layers is in KM !!
    for j =1:N
         Zb(:,j)=S(:,2+j);
         Zt(:,j)=S(:,3+j); 
    end
    
%--- CALCULATE PRESENT-DAY LAYER THICKNESS (T)
   for j = 1:N
         T(:,j)= Zb(:,j) - Zt(:,j) ;
   end
   
%--- PREPARE VARIABLES FOR DECOMPACTION
%--- ZZt and ZZb will be iteratively updated between each time steps
%--- For the first time step, ZZt and ZZb correspond to present-day depths
%--- For the following time steps ZZt and ZZb will correspond to decompacted new depths corrected from isostasy

   Ratio = zeros(N);
   for j = 1:N
         ZZt(:,j) = Zt(:,j)-Zt(:,N); %remove present bathymetry to move up all levels so the first one = 0
         ZZb(:,j) = Zb(:,j)-Zt(:,N); %remove present bathymetry to move up all levels to be consistent with ZZt
         Ratio(j) = Surpor(j)/C(j);  % Ratio between surf. poros. and decompaction coeff.
   end
   
    
%---------------------------------------------------------------------------------------
% THERMAL SUBSIDENCE, DYNAMIC TOPOGRAPHY AND SEALEVEL CORRECTION ON INDIVIDUAL LAYERS
%
% Calculation occurs outside the main backtracking loop since those
% corrections do not depend on layers but on time only. They are applied to newly
% decompacted layers at each time steps at the end of the backtracking process
%---------------------------------------------------------------------------------------

%--- CALCULATE WATER LOAD DUE TO SEA LEVEL CHANGES
   % Sea level changes
   [SL,outputstring]=sealevel(Np,N,opt_sl,dS,ages,Xcoord,Ycoord,SLfile,SLdir);%sealevel() ; % to compute sealevel correction and add later to tectonic subsidence of horizons

    %"Water load due to sealevel change"
   [dSL_iso]=Iso_water_load(Np,N,Xcoord,Ycoord,opt_sl,SL,iso_opt,grid_opt,EET_opt,EETfile,EETdir,EET,YM,v,Rhom,Rhow,Rhoc);

%--- CALCULATE THERMAL SUBSIDENCE   
   
   [ts_corr,Sfactor]=Thermal_subsidence(NN,Np,S,Xcoord,Ycoord,opt_subs,opt_BF,user_beta,Blimit,B1,B2,...
       filebeta,dirbeta,tlitho,tmoho,alpha,kappa,rift,Rhom,Rhoc,Rhow,ages) ;  
      
%--- CALCULATE DYNAMIC TOPOGRAPHY
   [dyntopo_corr,outputstring]=dynamic_topo(Np,N,Xcoord,Ycoord,ages,opt_dyn,dtopo_cte,dtopo_file,dtopo_dir,grid_opt);
    


%--------------------------------------------------------------------------------------
% PERFORM BACKTRACKING
%--------------------------------------------------------------------------------------
    Nzb=[];
    Nzt=[];
    Nzb=zeros(Np,N,N);          % Decompacted bottom depths
    Nzt=zeros(Np,N,N);          % Decompacted top depths
    Nzb_iso=zeros(Np,N,N);      % Decompacted top depths with isostatic correction
    Nzt_iso=zeros(Np,N,N);      % Decompacted bottom depths with isostatic correction    
    Nzb_corr=zeros(Np,N,N);     % Decompacted top depths with isostasy, thermal subsidence, sea level and dynamic topography
    Nzt_corr=zeros(Np,N,N);     % Decompacted top depths with isostasy, thermal subsidence, sea level and dynamic topography
        
    Dt=zeros(Np,N,N);           % Decompacted thicknesses
    Dt_cum=zeros(Np,N);         % Cumulative decompacted thicknesses
    Thick_cum=zeros(Np,N);      % Cumulative initial thicknesses
    WT_corr = zeros(Np,N,N);    % Iterative water thickness increment of each layer

    Por=zeros(Np,N,N);          % Porosity of newly decompacted layers
    Rhoa_units=zeros(Np,N);     % Density of newly decompacted layers
    Rhoa=zeros(Np,N);           % Average density of newly decompacted total sediment thickness
    Iso_units=zeros(Np,N,N);    % Isostatic depth correction resulting from each newly decompacted layer
    Tect=zeros(Np,N);           % Tectonic subsidence of the basement at each time step   
      
    WC(:)=Zt(:,N);              % Present-day water depth (bathymetry)


%------------ MAIN LOOP FOR BACKTRACKING -----------------------------------------------------------------%   
for t=1:N 

%---------------------------
% 1- DECOMPACTION
%---------------------------

%--- CALL DECOMPACTION
Decompaction();

%--- CALCULATE DECOMPACTED THICKNESSES
   for j = N-t+1:-1:1
      for k=1:Np
          Dt(k,t,j) = Nzb(k,t,j) - Nzt(k,t,j);
      end
   end

%--- CALCULATE CUMULATIVE DECOMPACTED THICKNESSES   
   temp1 = 0;
   for j = N-t+1:-1:1
      temp1 = temp1 + Dt(:,t,j);
   end
   Dt_cum(:,t) = temp1;   
   
%--- CALCULATE WATER THICKNESS INCREMENT IN NEWLY DECOMPACTED LAYERS   
   for j = N-t+1:-1:1
       WT_corr(:,t,j)= Dt(:,t,j)-(ZZb(:,j) -ZZt(:,j));
   end
  

%--------------------------------------------
% 2- POROSITY OF NEWLY DECOMPACTED LAYERS
%--------------------------------------------    
   Porosity();

%-------------------------------------------
% 3- DENSITY OF NEWLY DECOMPACTED LAYERS
%-------------------------------------------    
   Density_units();
        
   Density();
   
%-------------------------------------------
% 4- ISOSTATY OF NEWLY DECOMPACTED LAYERS
%-------------------------------------------     
   [Iso_layers,EET_out]=Isostasy_units(t,Np,N,Xcoord,Ycoord,Dt,Rhoa_units,iso_opt,grid_opt,EET_opt,EETfile,EETdir,EET,YM,v,Rhom,Rhow,Rhoc) ;% isostasy correction after removing a single layer at once
   Iso_units(:,t,:)=Iso_layers(:,t,:);


%---------------------------------------------------------------------
% 5- CORRECT DECOMPACTED DEPTHS WITH ITERATIVE ISOSTATIC CORRECTION
%
% From top to bottom, the upper unit is removed (Iso_units(:,t,N-t+1)) 
% and the underlying layers are isostatically relaxed
%
% The process is repeated iteratively from one time step to another
% and decompaction increment (water thickness) resulting from porosity
% changes are accounted for (WT_corr).
%---------------------------------------------------------------------


   for j= N-t+1:-1:1 

           if t < 2 % First time sept (only one layer removed)
              if j < 2
                 Nzb_iso(:,t,j) = Nzb(:,t,j)- Iso_units(:,t,N-t+1); % basement
                 Nzt_iso(:,t,j) = Nzt(:,t,j)- Iso_units(:,t,N-t+1); % basement
              else
                 Nzb_iso(:,t,j) = Nzb(:,t,j)- Iso_units(:,t,N-t+1)-WT_corr(:,t,j-1); % overlying layers
                 Nzt_iso(:,t,j) = Nzt(:,t,j)- Iso_units(:,t,N-t+1)-WT_corr(:,t,j-1); % overlying layers
              end
           else 
               if j < 2
                 Nzb_iso(:,t,j) =  Nzb_iso(:,t-1,j)- Iso_units(:,t,N-t+1); % basement
                 Nzt_iso(:,t,j) =  Nzt_iso(:,t-1,j)- Iso_units(:,t,N-t+1); % basement
               else
                 Nzb_iso(:,t,j) =  Nzb_iso(:,t-1,j)- Iso_units(:,t,N-t+1)-WT_corr(:,t,j-1); % overlying layers
                 Nzt_iso(:,t,j) =  Nzt_iso(:,t-1,j)- Iso_units(:,t,N-t+1)-WT_corr(:,t,j-1); % overlying layers                  
               end

           end
   end
 
   
%--- ADDITIONAL CORRECTIONS: Thermal subsidence (ts_corr), Sea Level, Dynamic Topography
% Equation for backtracking:  
% WD(t) = Z(t) - (rho_m-rho_savg)/(rho_m-rho_w)* S(i) + dSL(t)*rho_m/(rho_m-rho_w) + dDYN_topo(t)

    for j= N-t+1:-1:1
        Nzt_corr(:,t,j) = Nzt_iso(:,t,j) + ts_corr(:,t)./1000.  + dSL_iso(:,t)./1000. + dyntopo_corr(:,t)./1000. + WC(:);
        Nzb_corr(:,t,j) = Nzb_iso(:,t,j) + ts_corr(:,t)./1000.  + dSL_iso(:,t)./1000. + dyntopo_corr(:,t)./1000. + WC(:);    
    end    
    
  % Correct for Dt = 0 --> unconformities
  for i=1:t
   for j= N-t+1:-1:2 
     for k=1:Np
       if (Nzb_corr(k,t,j) > Nzb_corr(k,t,j-1))
          Nzb_corr(k,t,j) =  Nzb_corr(k,t,j-1);
       end 
     end 
   end
  end
  
%--- TOTAL TECTONIC SUBSIDENCE
% Corresponds to the decompacted and corrected basement (j=1)
   Tect(:,t) = -Nzb_corr(:,t,1);
%   figure
%   plot(Tect(:,t));hold on 
    
    
%--- UPDATE: ZZt and ZZb vriables for next iteration of decompaction
    for j= N-t+1:-1:1
        ZZt(:,j)= Nzt_iso(:,t,j);
        ZZb(:,j)= Nzb_iso(:,t,j);     
    end
   
    
end %%END MAIN LOOP ON BACKTRACKING 



%-------------------------------------------------------------------
%%% OUTPUTS TO BE SAVED OR PLOTTED %%% 
%-------------------------------------------------------------------
   
   out.Nzb         = -Nzb_corr.*1000.;  % backstripped depths (in meters)
   out.Dt          = Dt*1000;           % decompacted thicknesses (in meters)
   out.T           = T*1000.;           % Initial thicknesses
   out.Rhoa        = Rhoa;              % Averaged density
   out.Rhoa_units  = Rhoa_units;        % Density of each units
   out.Iso_units   = -Iso_units*1000.;  % Isostasy of each units
   out.Por         = Por*100;           % Porosity of each units
   out.tsubs_corr  = ts_corr;          % Thermal subsidence
   out.beta_factor = Sfactor;       % Streching factor
   out.Tect        = Tect.*1000.;              % Tectonic Subsidence
   
   % INITIAL DATA AND CONDITIONS
   out.EET         = EET_out;              % Elastic thickness
   out.SL          = SL;               % Sea level change
   out.dyntopo     = dyntopo_corr;           % Initial dynamic topography
   
   outputstring='';
  
   
%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    function Decompaction()
       eps = 0.000001;  
       for j= N-t+1:-1:1
           
             if j == N-t+1
                Nzt(:,t,j) = 0. ; % Put first level at 0 m.
             else
                Nzt(:,t,j) = Nzb(:,t,j+1);
             end
             
             for k=1: Np
                  tZb(k) = 1.;
                  delta(k) = 2*eps;
                  while (delta(k) > eps)
                       Nzb(k,t,j) = Nzt(k,t,j) + ( ZZb(k,j) - ZZt(k,j) ) - Ratio(j) * ( exp(-C(j).*ZZt(k,j)) - exp(-C(j).*ZZb(k,j))) + ... 
                           Ratio(j)* ( exp(-C(j).*Nzt(k,t,j)) - exp(-C(j).*tZb(k)) );

                       delta(k) = abs(Nzb(k,t,j)-tZb(k))./tZb(k);
                       tZb(k) = Nzb(k,t,j);
                  end % end while      
             end % end k loop
       end % end j loop
    end
%----------------------------------------------------------------------------
% POROSITY FOR EACH UNITS
%----------------------------------------------------------------------------

    function Porosity()
    
      for j=  N-t+1:-1:1
           Por(:,t,j) = Surpor(j);
        for k=1:Np
         %  Por(k,i,j) = Ratio(j)*exp(-C(j)*Nzt(k,i,j))*(1 - exp(-C(j)*Dt(k,i,j)));
          if (Dt(k,t,j) > 0)
            Por(k,t,j) = Ratio(j)*(exp(-C(j)*Nzt(k,t,j))- exp(-C(j)*Nzb(k,t,j)))./Dt(k,t,j);
          end
        end
      end          
          
    end

%----------------------------------------------------------------------------
% DENSITY FOR EACH LAYERS
%----------------------------------------------------------------------------

    function Density_units()

        for j=N-t+1:-1:1
            for k=1:Np
                 Rhoa_units(k,t,j) = Rhos(j) + (Rhow - Rhos(j)).*Por(k,t,j);
            end
        end
        
    end

%----------------------------------------------------------------------------
% AVERAGE DENSITY OVER TOTAL SEDIMENT THICKNESS
%----------------------------------------------------------------------------

function Density()

    % Call to single unit density:
    Density_units()

    rtemp = 0.;      
    for j= N-t+1:-1:1      
          rtemp = rtemp + Rhoa_units(:,t,j).*Dt(:,t,j);
    end

    Rhoa(:,t) = rtemp./Dt_cum(:,t) ;

end % END FUNCTION DENSITY


end % END BACKTRACKING