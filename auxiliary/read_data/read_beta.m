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
% Read spatially variable Beta factor from text files.
%
% Input:
%---------
% file.dat with [X Beta] for 2D and [X Y Beta] for 3D
%
% Output:
%---------
% M: matrix of form:
% X1 Y1 Z1_1 Z1_2 Z1_3 ... Z1_m
% .
% .
% .
% Xn Yn Zn_1 Zn_2 Zn_3 ... Zn_m
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
%-----------------------------------------------------------------------------------
function [Beta, outputstring]=read_beta(filebeta,betadir,Xcoord)


if nargin==1
    filebeta='';
end

outputstring='';

delimiterIn = ' ';
format long g

fid=fopen(fullfile(betadir,filebeta));
    if fid <0
        %disp(['Unable to open file ' fullfile(foldername,cfilenames{ii}) '.'])
        outputstring=['Unable to open file ' fullfile(betadir,filebeta) '.'];
        Beta=[];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==3
        data=importdata(fullfile(betadir,filebeta),delimiterIn);
    elseif number_of_columns==2
        data=importdata(fullfile(betadir,filebeta),delimiterIn);
    else
        fclose(fid);
        outputstring='Wrong number of columns in the text file.';
        Beta=[];
        return
    end
    
    

    if number_of_columns==3    
        Beta=data;
    else
        Beta=[data(:,1) zeros(size(data,1),1) data(:,2)]; %add a fictitious null column
    end
         
        
    % for line data, the coordinates may be different from layer to
    % layer.
    % If this is the case, the remaining layers are interpolated at the
    % coordinates of the first.

    if number_of_columns==2
        newcolumn = interp1(data(:,1),data(:,2),Xcoord,'nearest','extrap');% Data are interpolated on horizons X coordinates
    else % no interpolation for grid data
        newcolumn=data(:,3);
    end
    try
        Beta=[Beta newcolumn]; %#ok grow
    catch
        fclose(fid);
            outputstring='Wrong number of rows in the .dat files.';
        Beta=[];
    return
    end
        
fclose(fid);

end


