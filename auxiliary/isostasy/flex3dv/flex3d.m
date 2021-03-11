%-----------------------------------------------------------------------------------
% Copyright (C) Nestor Cardozo 2009
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
% Deflection of an elastic plate of variable thickness
% Eq. 3.83 of Ventsel and Krauthammer, 2001: Thin plates and shells
% Solution by centered finite differences
%
% Adapted for PALEOSTRIP by Colleoni et al. (2021) with the agreement of
% Nestor Cardozo
%
%------------------------------------------------------------------------------------
function [w_2d]=flex3d(tload,te,E,v,pm,pfill,dx,nx,ny)


% ---------------------------------------------------
% MODIFY THE PARAMETERS BELOW TO SUIT YOUR PROBLEM:
% As they are now, they correspond to test1.m -> test6.m
% ---------------------------------------------------
% Number of points in x
pointsx = nx ;
% Number of points in y
pointsy = ny ;
% Grid size or distance between points in meters
delta = dx ;
% Gravity in m/s^2
g = 9.8;


% -------------------------------------------------------
% Support of foundation
delp = (pm-pfill)*g;

% Compute flexural rigidity
D = (te.*1000).^3.0*E/(12.0*(1.0-v^2.0));

% --------------------------------------------------------------------
% Make row, column, and values vectors for non-zero entries of design
% matrix
% --------------------------------------------------------------------
sizeind= pointsx*4+(pointsy-4)*4+(pointsx-4)*(pointsy-4)*13;
rowi = zeros(sizeind,1);
coli = zeros(sizeind,1);
values = zeros(sizeind,1);
% --------------------------------------------------------------------
%  Indexes and values of non zero entries of design matrix
% --------------------------------------------------------------------

% NOTE: The program assumes that in the two columns and rows most proximal
% to the boundaries of the plate, the displacement is zero. This will work 
% if the loads are far from the boundaries: i.e. CONSTRUCT YOUR LOADS
% SUCH AS THAT THEY ARE IN THE CENTER OF THE PLATE, AWAY FROM THE
% BOUNDARIES

count = 1;
% Two outermost rows
for i=1:2
    for j=1:pointsx
        % current node
        nodek = (i-1)*pointsx+j;
        rowi(count)=nodek;
        coli(count)=nodek;
        values(count)=1.0;
        count = count+1;
    end
end
for i=pointsy-1:pointsy
    for j=1:pointsx
        % current node
        nodek = (i-1)*pointsx+j;
        rowi(count)=nodek;
        coli(count)=nodek;
        values(count)=1.0;
        count = count+1;
    end
end
% Two outermost columns
for j=1:2
    for i=3:pointsy-2
        % current node
        nodek = (i-1)*pointsx+j;
        rowi(count)=nodek;
        coli(count)=nodek;
        values(count)=1.0;
        count = count+1;
    end
end
for j=pointsx-1:pointsx
    for i=3:pointsy-2
        % current node
        nodek = (i-1)*pointsx+j;
        rowi(count)=nodek;
        coli(count)=nodek;
        values(count)=1.0;
        count = count+1;
    end
end
% Inner nodes
for i=3:pointsy-2
    for j=3:pointsx-2
        % Node indexes for computing finite differences
        % current node
        nodek = (i-1)*pointsx+j;
        % surrounding nodes
        noden = nodek-(pointsx*2);
        nodeh = nodek-pointsx-1;
        noded = nodek-pointsx;
        nodeg = nodek-pointsx+1;
        nodei = nodek-2;
        nodea = nodek-1;
        nodec = nodek+1;
        nodem = nodek+2;
        nodee = nodek+pointsx-1;
        nodeb = nodek+pointsx;
        nodef = nodek+pointsx+1;
        nodel = nodek+(pointsx*2);
        % Compute constants and derivatives in D
        A = D(nodek);
        B = 2.0*((D(nodec)-D(nodea))/(2.0*delta));
        C = 2.0*((D(nodeb)-D(noded))/(2.0*delta));
        DD = (D(nodea)+D(nodeb)+D(nodec)+D(noded)-4.0*D(nodek))/(delta^2.0);
        E = -(1.0-v)*((D(nodea)-2.0*D(nodek)+D(nodec))/(delta^2.0));
        F = 2.0*(1.0-v)*((D(nodeh)-D(nodee)+D(nodef)-D(nodeg))/(4.0*delta^2.0));
        G = -(1.0-v)*((D(noded)-2.0*D(nodek)+D(nodeb))/(delta^2.0));
        H = delp;
        % indexes and values of design matrix
        for k=1:13
            rowi(count)=nodek;
            if k == 1
                coli(count)=noden;
                values(count)=A/(delta^4.0)-C/(2.0*delta^3.0);
            elseif k == 2
                coli(count)=nodeh;
                values(count)=2.0*A/(delta^4.0)-B/(2.0*delta^3.0)-C/(2.0*delta^3.0)+F/(4.0*delta^2.0);
            elseif k == 3
                coli(count)=noded;
                values(count)=-8.0*A/(delta^4.0)+2.0*C/(delta^3.0)+DD/(delta^2.0)+E/(delta^2.0);
            elseif k == 4
                coli(count)=nodeg;
                values(count)=2.0*A/(delta^4.0)+B/(2.0*delta^3.0)-C/(2.0*delta^3.0)-F/(4.0*delta^2.0);
            elseif k == 5
                coli(count)=nodei;
                values(count)=A/(delta^4.0)-B/(2.0*delta^3.0);
            elseif k == 6
                coli(count)=nodea;
                values(count)=-8.0*A/(delta^4.0)+2.0*B/(delta^3.0)+DD/(delta^2.0)+G/(delta^2.0);
            elseif k == 7
                coli(count)=nodek;
                values(count)=20.0*A/(delta^4.0)-4.0*DD/(delta^2.0)-2.0*E/(delta^2.0)-2.0*G/(delta^2.0)+H;
            elseif k == 8
                coli(count)=nodec;
                values(count)=-8.0*A/(delta^4.0)-2.0*B/(delta^3.0)+DD/(delta^2.0)+G/(delta^2.0);
            elseif k == 9
                coli(count)=nodem;
                values(count)=A/(delta^4.0)+B/(2.0*delta^3.0);
            elseif k == 10
                coli(count)=nodee;
                values(count)=2.0*A/(delta^4.0)-B/(2.0*delta^3.0)+C/(2.0*delta^3.0)-F/(4.0*delta^2.0);
            elseif k == 11
                coli(count)=nodeb;
                values(count)=-8.0*A/(delta^4.0)-2.0*C/(delta^3.0)+DD/(delta^2.0)+E/(delta^2.0);
            elseif k == 12
                coli(count)=nodef;
                values(count)=2.0*A/(delta^4.0)+B/(2.0*delta^3.0)+C/(2.0*delta^3.0)+F/(4.0*delta^2.0);
            elseif k == 13
                coli(count)=nodel;
                values(count)=A/(delta^4.0)+C/(2.0*delta^3.0);
            end
            count = count+1;
        end
    end
end

% --------------------------------------------------
% Construct design matrix as a sparse matrix
% --------------------------------------------------
totnod = pointsx*pointsy;
DM = sparse(rowi,coli,values,totnod,totnod);

% --------------------------------------------------
% Compute deflections
% --------------------------------------------------
w = DM\tload;

% -------------------------------------------------
% Write deflections
% -------------------------------------------------
% %fid = fopen('deflection.txt','wt');
count = 1;
for i=1:pointsy
    for j=1:pointsx
       % % fprintf(fid,'%f\n',w(count));
       if(w(count)> 2000 || w(count) < -2000)
           w(count)=0.;
       end
        w_2d(i,j) = w(count);
        count = count+1;
    end
end        
% %fclose(fid);
end