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
% Read dynamic topography 2D maps and interpolate on original sediment layers
% input grid and on corresponding seidment horizon ages
% This routine consents multiple input dynamic topography maps and thus
% time-evolving 2D dynamic topographic correction.
%
% Input:
%---------
% dtop_file.dat with [X Y DynTopo] input from GUI
% all input files are placed within cell array of filenames, 
% example {'file_1.dat','file_2.dat',...,'file_m.dat'}
%
% DynTopo is provided in m.
%
% Ages (in Ma) of each dtopo_file.dat are MANDATORY and prescribed from the GUI 
%
% NOTE: it is important that input dynamic topo data are provided on cartesian
% coordinates grid and with the same geographical projection and same horizontal 
% resolution than input sediment layers. But grid extent and shape can
% differ
%
%
% Example of formatted input text file:
% X1 [Y1] Z1
% X2 [Y2] Z2
% .
% .
% .
% Xn [Yn] Zn
%
%
% When there is no [Y] column, the second column of M is set to 0 (fictitious
% column)
%
% 
% Dynamic topography maps are interpolated at corresponding Ages of input
% sediment layers and used as correction for decompated sediment horizon
% depths
%
%-----------------------------------------------------------------------------------

function [dt_DT, outputstring]=read_dtopo(grid_opt,Np,N,Xcoord,Ycoord,Hages,dtopo_file,dtopo_dir)


% Get number of files to be treated
D = readcell(fullfile(dtopo_dir,dtopo_file));
dynAges=cell2mat(D(:,1));
dyn_files=D(:,2);
nfiles=size(D,1);


if nargin==1
    dtopo_dir='';
end

outputstring=' ';

delimiterIn = ' ';
format long g

for ii=1:nfiles
    fid=fopen(fullfile(dtopo_dir,dyn_files{ii}));
    if fid <0
        outputstring=['Unable to open file ' fullfile(dtopo_dir,dyn_files{ii}) '.'];
        DT=[];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==3
       if grid_opt==3
         dtdata=importdata(fullfile(dtopo_dir,dyn_files{ii}),delimiterIn);
       else
        fclose(fid);
        outputstring='Input data must have 3 columns - Select Grid dimension: 3D';
        DT=[];
        return         
       end
    elseif number_of_columns==2
       if grid_opt==2
        dtdata=importdata(fullfile(dtopo_dir,dyn_files{ii}),delimiterIn);
        
       else
        fclose(fid);
        outputstring='Input data have 2 columns - Select Grid dimension: 2D';
        DT=[];
        return         
       end                 
    else
        fclose(fid);
        outputstring='Use uniform dyn. topo. correction option';
        DT=[];
        return
    end
    
    
    if ii==1
        if number_of_columns==3
            [dtopo_interp]=mesh3D_inputvar(Xcoord,Ycoord,dtdata(:,1),dtdata(:,2),dtdata(:,3));
            DT=dtopo_interp(:);

        else
            [X,index]=unique(dtdata(:,1));
            orig_data=dtdata(:,2);
            DT=interp1(X,orig_data(index),Xcoord,'linear','extrap');
        end
        
    else
     
    % Populate matrix DT by adding new dyn topo column
    % Interpolate on sediment data grid.
        if number_of_columns==2
            %Interpolate on sediment data grid
            newcolumn = interp1(dtdata(:,1),dtdata(:,2),Xcoord,'linear','extrap');
        elseif number_of_columns==3 % no interpolation for grid data
            [dtopo_interp]=mesh3D_inputvar(Xcoord,Ycoord,dtdata(:,1),dtdata(:,2),dtdata(:,3));
            newcolumn=dtopo_interp(:);
            %newcolumn(1)
        end
        try
            DT=[DT newcolumn];
        catch
            fclose(fid);
                outputstring='Wrong number of rows in the .dat files.';
            DT=[];
        return
        end
        
    end
    
    %
    
end % End Loop read data files


% Interpolate at sediment layer ages

% Sediment layer ages
HAges = [Hages{:}];
tHages = HAges.';

if (grid_opt == 2)
    for i=1:nfiles
       DT_corr(:,i) = DT(:,i) - DT(:,1);% compute correction relative to present-day (2nd column)
    end   

elseif (grid_opt == 3)
    for i=1:nfiles
       DT_corr(:,i) = DT(:,i) - DT(:,1);% compute correction relative to present-day (2nd column)
    end
end

%Interpolate on sediment ages
[dt_DT] = time_interp(nfiles,N,dynAges,tHages,DT_corr);
dt_DT=fliplr(dt_DT); % To get it in chronological order for correction on backtracking

    
end