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
% ISOSTASY: WATER LOAD 
%
% Two methods are implemented:
% - Airy local isostasy
% - Flexural isostasy
%
% In the case of flexural isostasy:
% - 2D: Chapman and Cardozo () - Finite difference grid and variable Te
% - 3D: N. Cardozo () - Finite difference grid and variable Te
%
% Flexural isostasy requires calculations: 
% - on distances (m) and not on cartesian coordinates thus the initial vector or matrix is
% transformed to get distance within the routines flex2d, flex3dv.
%
% - on a broader domain to avoid loads interferences at edges of the domain:
% matrix and vectors are interpolated on REGULAR grid broader than the
% original one and reinterpolated on original coordinates after flexure has
% been calculated.
%
%----------------------------------------------------------------------------     
function [dSL_iso]=Iso_water_load(Np,N,Xcoord,Ycoord,opt_sl,SL,iso_opt,grid_opt,EET_opt,EETfile,EETdir,EET,YM,v,Rhom,Rhow,Rhoc) 
        

    
  dSL_iso=zeros(Np,N);
    
  if (opt_sl < 6 ) % If sea level change is switch on
        
%----------------------------------------------------------------------
% AIRY LOCAL COMPENSATION
%----------------------------------------------------------------------        
      if (iso_opt == 1)
         for t=1:N           
          dSL_iso(:,t) = -SL(:,t).*Rhom/(Rhom-Rhow) ;
         end   


%----------------------------------------------------------------------
% 2D FLEXURE
%----------------------------------------------------------------------
      elseif (iso_opt == 2) % Flexure         
          if (grid_opt == 2) % Flexure 2D
              
            %%%%%%% Get EET vector %%%%%%%%%%
            dx_grid=ceil(abs(Xcoord(2)-Xcoord(1))/1000);
            xpos=1:dx_grid:Np*dx_grid;
            xnew=1:dx_grid:Np*dx_grid*3; % Expand domain for isostasy
            NX=length(xnew);
            
            % Check on Te option
            if (EET_opt == 1) % Variable Te                                                 
                [EET_1d]=read_EET_3d(EETfile,EETdir);
                EET_out(1:Np) = EET_1d(:);
            else % Spatially uniform Te
                EET_1d(1:Np)=EET;
            end
            
            % Expand EET to the broader domain
            EET_new(1:Np*3) = 0.; 
            EET_new(Np+1:Np*2)= EET_1d(1:Np);
            EET_new(Np-ceil(Np*0.9):Np) = EET_1d(1); 
            EET_new(Np*2+1:Np*2+ceil(Np*0.9)) = EET_1d(Np);             
                          
         g = 9.81; % Gravity constant
         for t=1:N
                % Shift values in the middle of the array to avoid edges artifacts
                % set to 0 the entire array
                % populate the array at the middle with decompacted depths
                dSL_new(1:Np*3,t) = 0.; 
                dSL_new(Np+1:Np*2,t)= SL(:,t);
                dSL_new(Np-round((Np)*0.4):Np,t) = SL(1,t); 
                dSL_new(Np*2+1:Np*2+round(Np*0.4),t) = SL(Np,t);             

                % Create X vector using GUI spacing value
                x_vector=0:dx_grid:xnew(NX);
                
                % Calculate input load on X vector
                EET_interp(:) = interp1(xnew,EET_new(:),x_vector,'linear','extrap');
                dSL_interp(:) = interp1(xnew,dSL_new(:,t),x_vector,'linear','extrap');                           
                q_sl(:)=dSL_interp(:).*(Rhom).*g;

                %CALCULATE 2D-FLEXURE - Infinite plate and distributed load            
                [SL_1d]=flex2d(q_sl,x_vector,Rhom,Rhow,EET_interp,YM,v,dx_grid);
                 sl_interp = interp1(x_vector,SL_1d,xnew,'linear',0);
                 sl_total= sl_interp(Np+1:Np*2); % Retrieve values
                 
                 
             end
                dSL_iso(:,t)=-sl_total(1:Np);
         
%----------------------------------------------------------------------
% 3D FLEXURE
%----------------------------------------------------------------------          
          
          elseif (grid_opt == 3) % Flexure 3D

              %%%%%%%% Effective Elastic Thickness %%%%%%%%%%%%%%%%%
              %Interpolate variable EET on the data grid & CREATE MESH and get 1D vector from 2D new matrix
              if (EET_opt == 1) % Variable EET (N. Cardozo - fled3dv code)
                [EET, Xeet,Yeet,outputstring]=read_EET_3d(EETfile,EETdir);
                [EET_1d,EET_out]=mesh3D_eet(Xcoord,Ycoord,Xeet,Yeet,EET);
                 
              end

             %%%%%%%% Effective Elastic Thickness %%%%%%%%%%%%%%%%%
              for t=1:N                      
                    
                 % if (EET_opt == 1) % Variable EET - Finite difference method
                    q_sl(:)=-SL(:,t).*(Rhom)*9.81;  %Calculate load magnitude (N/m2)

                    %CREATE MESH and get 1D vector from 2D new matrix
                    [SL_1d,mask_1d,NX,NY,dx,dy,X2d,Y2d,index]=mesh3D_load(Xcoord,Ycoord,q_sl);
                    
                    if (EET_opt == 2) % Uniform EET
                        EET_1d(1:length(mask_1d)) = mask_1d(1:end).*EET;
                    end
                    
                    %COMPUTE FLEXURE  - Infinite plate and distributed load
                    [SL_2d]=flex3d(SL_1d,EET_1d,YM,v,Rhom,Rhow,dx,NX,NY);
                  
                    SL_w1d = -SL_2d(:);
  
                    sl_total = SL_w1d(index(:)); % Extract values only at original Xcoord and Ycoord

                    dSL_iso(1:Np,t)= -sl_total(1:Np);

              end % FOR LOOP
          end % END GRID OPTION
      end % END ISOSTASY OPTION
      
  else % NO ISOSTASY
      for t=1:N
         dSL_iso(:,t) = 0. ;
      end
  end % END SEA LEVEL OPTION

 end %END FUNCTION