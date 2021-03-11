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
%-------------------------------------------------------------------------
% READ IMPORTED COORDINATES FROM INPUT FILE TO BE EXTRACTED FROM 3D maps
%
% Files must contain cartesian coordinates on two columns (no headers):
% X Y
%
% This can be for extracting:
% - a well (1D): only one set of (X,Y) is necessary
% - a 2D transect: a list of (X,Y)
%
% The routine automatically detects if 1D or 2D.
%-------------------------------------------------------------------------
function [DataPoints outputstring]=read_imported_coord(filecoord,coorddir)
% Read X Y from text files.


if nargin==1
    filecoord='';
end

outputstring='';

delimiterIn = ' ';
format long g

fid=fopen(fullfile(coorddir,filecoord));
    if fid <0
        %disp(['Unable to open file ' fullfile(foldername,cfilenames{ii}) '.'])
        outputstring=['Unable to open file ' fullfile(coorddir,filecoord) '.'];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==2
        DataPoints=importdata(fullfile(coorddir,filecoord),delimiterIn);
    else
        fclose(fid);
        outputstring='File should contain only 2 columns (X and Y).';
        return
    end

fclose(fid);

end