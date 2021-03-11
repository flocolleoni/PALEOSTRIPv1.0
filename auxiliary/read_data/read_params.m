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
% READ LAYERS LITHOLOGICAL PARAMETER FROM INPUT FILE
%
% From Bottom (layer n°1) to top (layer n)
%
% Example of formatted text: 
%
% [empty line]
% NUMBER OF LAYERS =           8
% [empty line]
%   LAYER  POROSITY   DEC CON   MAT DEN  AGE BASE     NAME
%                     (1/KM)   (KG / MC)   (MA)    
% [empty line]
%      1    0.45       0.450     2680    30.000     layer1name
%      2    0.45       0.450     2680    26.000     layer2name
%      3    0.45       0.450     2680    21.000     layer3name
%      4    0.45       0.450     2680    19.500     layer4name
%      5    0.45       0.450     2680    16.000     layer5name
%      6    0.45       0.450     2680     8.000     layer6name
%      7    0.45       0.450     2680     4.000     layer7name
%      8    0.45       0.270     2680     0.600     layer8name
%-------------------------------------------------------------------------

function params=read_params(filename)

if ~exist(filename,'file')
    disp(['file ' filename ' not found'])
    params=[];
    return
end

fid=fopen(filename);

% get number of layers
number_of_layers=[];
while isempty(number_of_layers)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    number_of_layers=sscanf(tline,'NUMBER OF LAYERS = %d');
end

if  isempty(number_of_layers)
    disp('Invalid format')
    params=[];
    fclose(fid);
    return
end

%skip 4 lines
fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);

data=textscan(fid,'%d %f %f %f %f %s',number_of_layers);

if ( size(data{1},1)~=number_of_layers )
    disp('Invalid format')
    params=[];
    fclose(fid);
    return
end

for ii=1:number_of_layers
    params(ii).Surpor = data{2}(ii);
    params(ii).C = data{3}(ii);
    params(ii).Rhos = data{4}(ii);
    params(ii).AgeBase = data{5}(ii);
    params(ii).Name = data{6}{ii};
end
fclose(fid);