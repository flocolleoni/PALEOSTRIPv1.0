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
% Read input EET data from text files.
%
% Input
%--------
% file.dat with [X EET] for 2D and [X Y EET] for 3D
%
% NOTE 1: it is important that input EET data are provided on cartesian
% coordinates grid and with the same geographical projection and same horizontal 
% resolution than input sediment layers. But grid extent and shape can
% differ
%
% NOTE 2: For 1D EET use "uniform" option from the GUI
%
% Output
%--------
% M: matrix of form:
% X1 Y1 Z1_1 Z1_2 Z1_3 ... Z1_m
% .
% .
% .
% Xn Yn Zn_1 Zn_2 Zn_3 ... Zn_m
%
% When there is no [Y] column, the second column of M is set to 0 (fictitious
% column)
%-----------------------------------------------------------------------------------

function [EET, Xeet,Yeet,outputstring]=read_EET_data(fileEET,EETdir)

if nargin==1
    fileEET='';
end

outputstring='';

delimiterIn = ' ';
format long g

fid=fopen(fullfile(EETdir,fileEET));
    if fid <0
        outputstring=['Unable to open file ' fullfile(EETdir,fileEET) '.'];
        EET=[];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==3
        data=importdata(fullfile(EETdir,fileEET),delimiterIn);
        EET=data(:,3);
        Yeet=data(:,2);
        Xeet=data(:,1);
    elseif number_of_columns==2
        data=importdata(fullfile(EETdir,fileEET),delimiterIn);
        EET=data(:,2);
        Yeet=data(:,1);
        Xeet=data(:,1);
    else
        fclose(fid);
        outputstring='Wrong number of column in EET input file. If 1D well selected: use uniform EET opt.';
        EET=[];
        return
    end
    
%     figure
%     plot(EET)

fclose(fid);

end