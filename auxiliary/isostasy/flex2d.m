%-----------------------------------------------------------------------------------
% Copyright (C) 2015 Jay Chapman
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
% Flex2D script for variable flexural rigitidy in x direction
% Centered Finite Difference Technique
% uses flexural equation from Turcotte and Schubert, 1982
% d2/dx2(D(x)(d2w/dx2))+ (pm-pc)gw = q(x)
% written by Jay Chapman, Nov. 2015
%
% Adapted for PALEOSTRIP by Colleoni et al. (2021)
%
%------------------------------------------------------------------------------------
function [w]=flex2d(Load_1d,Xcoord,Rhom,Rhow,EET_interp,YM,PR,dx)


%Grid size, number of nodes in x direction
sizex = length(Xcoord);

%Node spacing (meters)
delta = dx*1e3;

% profile distance in km
x = [0:delta:delta*(sizex-1)];


% Young Modulus in Pascals
E = YM;
% Poisson's ratio
v = PR;
% Gravity acceleration in m/s^2
g = 9.81;
% Density of the mantle kg/m^3
pm = Rhom;
% Density of the material filling resultant depression kg/m^3
pfill = Rhow;

%initialize vector for elastic thickness
te = ones(1,sizex);  
te = EET_interp.*1e3.*te;


% Compute flexural rigidity
D = te'.^3.0*E/(12.0*(1.0-v^2.0));

% Compute load
tload = 1*Load_1d';


% --------------------------------------------------------------------
% These are variables that go into making a sparse matrix for the inversion
% "sizeind" is the number of non-zero values in the sparse matrix
% "rowi" and "coli" point the position of the non-zero element "values"
sizeind = ((sizex-4)*5)+4; 
rowi = zeros(sizeind,1);
coli = zeros(sizeind,1);


values = zeros(sizeind,1);



% These are the outer non-zero terms
count = 1;
for i = [1 2 sizex-1 sizex]
       %sets the diagonal on first two and last two rows in DM =1
       %corresponds to the first and last two entries in the sizex vector
        rowi(count)=i;
        coli(count)=i;
        values(count)=1.0;
        count = count+1;
end

%Nodes
%  a-b-c-d-e
% nodec is center node

% Inner nodes
for i=3:sizex-2
        % Node indexes for computing finite differences
        % current node
        nodec = i;

        % surrounding nodes
        nodea = nodec-2;
        nodeb = nodec-1;
        noded = nodec+1;
        nodee = nodec+2;
        
        % Compute constants and derivatives in D
        A = D(nodec);
        B = (pm-pfill)*g; 
                      
        %Nodes
      %  a-b-c-d-e
    %k=  1-2-3-4-5
        
    %coefficients for 4th derivative, 2nd order accruacy
    % 1 -4 6 -4 1
               
     % indexes and values of design matrix
        for k=1:5
            rowi(count)=nodec;
            if k == 1
                coli(count)=nodea;
                values(count)=A/(delta^4.0);
            elseif k == 2
                coli(count)=nodeb;
                values(count)=-4*A/(delta^4.0);
            elseif k == 3
                coli(count)=nodec;
                values(count)=6*A/(delta^4.0)+B;
            elseif k == 4
                coli(count)=noded;
                values(count)=-4*A/(delta^4.0);
            elseif k == 5
                coli(count)=nodee;
                values(count)=A/(delta^4.0);
            
            end
            count = count+1;
        end
    end
    

% --------------------------------------------------
% Design Matrix (DM)
%size of sparse matrix is sizex by sizex
DM = sparse(rowi,coli,values,sizex,sizex);

% --------------------------------------------------
% Direct inversion to calculate deflections (w)
w = DM\tload;

% -------------------------------------------------
% Uncomment to write deflections out to txt file
% fid = fopen('deflection.txt','wt');
% for i=1:sizex
%         fprintf(fid,'%f\n',w(i));
%     w_1d(:) = w(:);
% end        
% fclose(fid);

end
  
    