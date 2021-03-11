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
% Read sea level changes time series and interpolate at corresponding ages of 
% input sediment horizons.
%
% Input:
%---------
% file.dat with [Time Sealevel]
%
% Time is provided in Ma
% Sealevel changes are provided in m relative to present
%
% Imput sea level time series are interpolated linearily at the time of input
% sediment horizons and return for correction of decompated depths
%-----------------------------------------------------------------------------------

function [dSL, outputstring]=read_data_sl_1d(file_sl,data_folder,Hages)


if nargin==1
    data_folder='';
end


outputstring='';
    fid=fopen(fullfile(data_folder,file_sl));
    if fid <0
        outputstring=['Unable to open file ' fullfile(data_folder,file_sl) '.'];
        dSL=[];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==2
        datasl=textscan(fid,'%f %f');
        fclose(fid) ;
    else
        fclose(fid);
        outputstring='Timeseries of sea level changes must contain two colums (time,sl)';
        datasl=[];
        return
    end
    
    % To remove several occurrence of the same time
    xSL = [datasl{:,1}] ;
    ySL = [datasl{:,2}] ;
    [xSL, index] = unique(xSL) ;
    HAges = [Hages{:}];
    tHAges = HAges.';
   
    %Interpolate sea level from sea level curve on backstripped horizons ages
    dSL = interp1(xSL,ySL(index),tHAges,'linear','extrap');
    
    % Attribute 1st or last values if Ages are younger or older than the
    % 1st or last ages of the SL time series
    for i=1:length(tHAges)
        if (tHAges(i) > xSL(end))
            dSL(i) = ySL(end);
        elseif (tHAges(i) < xSL(1))
            dSL(i) = ySL(end);        
        end
    end

end


