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
% ISOSTASY: SINGLE SEDIMENTARY UNITS 
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
    function [Iso_layers,EET_out]=Isostasy_units(t,Np,N,Xcoord,Ycoord,Dt,Rhoa_units,iso_opt,grid_opt,EET_opt,EETfile,EETdir,EET,YM,v,Rhom,Rhow,Rhoc)  
        
    Iso_layers=zeros(Np,N,N);
    Iso_corr=zeros(Np,N,N);
    %EET_1d=zeros(Np);

      if (iso_opt == 1) % AIRY LOCAL COMPENSATION

            for j=N-t+1:-1:1
                  Iso_layers(:,t,j) = Dt(:,t,j).*(Rhom-Rhoa_units(:,t,j))./(Rhom - Rhow) ;
            end
            EET_out(1:Np) = 0.;
            
%---------------------------------------------------------------------------------------
% FLEXURE 2D
%---------------------------------------------------------------------------------------
      elseif (iso_opt == 2) % Flexure
          if (grid_opt == 2) % Flexure 2D

            % Get Effective Elastic Thickness vector
            dx_grid=ceil(abs(Xcoord(2)-Xcoord(1))/1000);
            xpos=1:dx_grid:Np*dx_grid;
            xnew=1:dx_grid:Np*dx_grid*3; % Expand domain for isostasy
            NX=length(xnew);
            % Check on Te option
            if (EET_opt == 1) % Variable Te                                                 
                [EET_1d]=read_EET_data(EETfile,EETdir);
                EET_out(1:Np) = EET_1d(:);
            else % Spatially uniform Te
                EET_1d(1:Np)=EET;
                EET_out(1:Np) = EET;
            end
            
            % Expand EET to the broader domain
            EET_new(1:Np*3) = 0.; 
            EET_new(Np+1:Np*2)= EET_1d(1:Np);
            EET_new(Np-ceil(Np*0.9):Np) = EET_1d(1); 
            EET_new(Np*2+1:Np*2+ceil(Np*0.9)) = EET_1d(Np);
            
            %figure
            
            g = 9.81;
            for j=N-t+1:-1:1
                
               % Shift values in the middle of the array to avoid edges artifacts
               % set to 0 the entire array
               % populate the array at the middle with decompacted depths
               Dt_new(1:Np*3,t,j) = 0.; 
               Dt_new(Np+1:Np*2,t,j)= Dt(:,t,j);
               Dt_new(Np-ceil(Np*0.9):Np,t,j) = Dt(1,t,j); 
               Dt_new(Np*2+1:Np*2+ceil(Np*0.9),t,j) = Dt(Np,t,j);


               Rhoa_units_new(1:Np*3,t,j) = 0.;
               Rhoa_units_new(Np+1:Np*2,t,j)= Rhoa_units(:,t,j); 
               Rhoa_units_new(Np-ceil(Np*0.9):Np,t,j) = Rhoa_units(1,t,j); 
               Rhoa_units_new(Np*2+1:Np*2+ceil(Np*0.9),t,j) = Rhoa_units(Np,t,j); 

               %Interpolate onto isostasy grid at spacing dx (from Isostasy GUI)
               x_vector=0:dx_grid:xnew(NX);
               EET_interp(:) = interp1(xnew,EET_new(:),x_vector,'linear','extrap');
               Dt_interp(:) = interp1(xnew,Dt_new(:,t,j),x_vector,'linear','extrap');
               Rhoa_units_interp(:) = interp1(xnew,Rhoa_units_new(:,t,j),x_vector,'linear','extrap');


               % Calculate surface load
               h(:) = Dt_interp(:)*1e3;
               q(:)=h(:).*(Rhom-Rhoa_units_interp(:)).*g;

               %CALCULATE 2D-FLEXURE - Infinite plate and distributed load
               %Jay Chapman/Cardozo
               [W_1d]=flex2d(q,x_vector,Rhom,Rhow,EET_interp,YM,v,dx_grid);
               w_interp = interp1(x_vector,W_1d,xnew,'linear',0);
               w_total= w_interp(Np+1:Np*2); % Retrieve values on original array

               Iso_layers(:,t,j)= w_total(1:Np)/1000.; 
                
            end

                     
%----------------------------------------------------------------------
% 3D FLEXURE
%----------------------------------------------------------------------
          elseif (grid_opt == 3) % Flexure 3D
              
              %%%%%%%% Effective Elastic Thickness %%%%%%%%%%%%%%%%%
              %Interpolate variable EET on the data grid & CREATE MESH and get 1D vector from 2D new matrix
              if (EET_opt == 1) % Variable EET 
                [EET, Xeet,Yeet,outputstring]=read_EET_data(EETfile,EETdir);
                [EET_1d,EET_out]=mesh3D_eet(Xcoord,Ycoord,Xeet,Yeet,EET);                
              end

               for j= N-t+1:-1:1 
                    h = Dt(:,t,j)*1e3;
                    q(:)=h.*(Rhom-Rhoa_units(:,t,j)).*9.81;  %Calculate load magnitude (N/m2)
                    

                    %CREATE MESH and get 1D vector from 2D new matrix
                    [Load_1d,mask_1d,NX,NY,dx,dy,X2d,Y2d,index]=mesh3D_load(Xcoord,Ycoord,q);
                    
                    if (EET_opt == 2) % Uniform EET
                        EET_1d(1:length(mask_1d)) = mask_1d(1:end).*EET;
                        EET_out(1:Np) = EET_1d(index(:));
                    end
                    %Calculate 3D flexure from Cardozo et al.
                    [W_2d]=flex3d(Load_1d,EET_1d,YM,v,Rhom,Rhow,dx,NX,NY);
                    W_1d = -W_2d(:);

                  end   

                     w_total = W_1d(index(:)); % Extract values only at original Xcoord and Ycoord
                     Iso_layers(:,t,j)= -w_total(1:Np)./1000;

               end %GRID OPT

                      
      elseif (iso_opt == 3) % NO ISOSTASY
            for j=N-t+1:-1:1
                  Iso_corr(:,t,j) = 0. ;
            end
            
            EET_out(1:Np) = 0.;
      end % END ISOSTASY OPTION
      
 end %END FUNCTION