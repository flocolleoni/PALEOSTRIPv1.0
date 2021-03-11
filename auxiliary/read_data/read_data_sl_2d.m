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
% Read sea level changes 2D map and interpolate on original sediment layers
% input grid
%
% Input:
%---------
% file.dat with [X Y Sealevel]
%
%
% NOTE: it is important that input sea level data are provided on cartesian
% coordinates grid and with the same geographical projection and same horizontal 
% resolution than input sediment layers. But grid extent and shape can
% differ
%
%-----------------------------------------------------------------------------------

function [dSL,outputstring]=read_sl_2d(Xcoord,Ycoord,sl_file,sl_dir,Hages,N)
% Read depth data from text files.


if nargin==1
    sl_file='';
end

outputstring='';

delimiterIn = ' ';
format long g

fid=fopen(fullfile(sl_dir,sl_file));
    if fid <0
        outputstring=['Unable to open file ' fullfile(sl_dir,sl_file) '.'];
        dSL=[];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==3
        datasl=importdata(fullfile(sl_dir,sl_file),delimiterIn);
        %Interpolate SL on original data grid
        [datasl_interp]=mesh3D_inputvar(Xcoord,Ycoord,datasl(:,1),datasl(:,2),datasl(:,3));       
    else
        fclose(fid);
        outputstring='Sea Level maps should contain 3 colums X Y Z';
        datasl=[];
        return
    end
     
    dSL =datasl_interp(:,3);


fclose(fid);

end