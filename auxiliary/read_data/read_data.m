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
% Read input sediment horizons depth data from text files.
%
% Input:
%--------
% Sediment horizon depths.
% cfilenames: cell array of filenames, example {'file_1.dat','file_2.dat',...,'file_m.dat'}
%
% Output
%-----------
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
function [M, outputstring]=read_data(cfilenames,foldername,grid_opt)


if nargin==1
    foldername='';
end

outputstring='';

delimiterIn = ' ';
format long g
for ii=1:length(cfilenames)
    fid=fopen(fullfile(foldername,cfilenames{ii}));
    if fid <0
        %disp(['Unable to open file ' fullfile(foldername,cfilenames{ii}) '.'])
        outputstring=['Unable to open file ' fullfile(foldername,cfilenames{ii}) '.'];
        M=[];
        return
    end
    
    % check number of columns
    tline = fgetl(fid);
    number_of_columns=length(str2num(tline));
    fseek(fid,0,-1); %rewind the file
   
    if number_of_columns==3
       if grid_opt==3
         data=importdata(fullfile(foldername,cfilenames{ii}),delimiterIn);
       else
        fclose(fid);
        outputstring='Input data have 3 columns - Select Grid dimension: 3D';
        M=[];
        return         
       end
    elseif number_of_columns==2
       if grid_opt==2
        data=importdata(fullfile(foldername,cfilenames{ii}),delimiterIn);
       else
        fclose(fid);
        outputstring='Input data have 2 columns - Select Grid dimension: 2D';
        M=[];
        return         
       end        
    elseif number_of_columns==1
       if grid_opt==1
        data1d=importdata(fullfile(foldername,cfilenames{ii}),delimiterIn);
        x=[1:1:100]; % Factice coordinates
        data=zeros(100,2);
        data(:,2)=data1d(:); %Expand the well to 2d uniform values to perform calculations
        data(:,1)=x(:); %Expand the well to 2d uniform values to perform calculations
       else
        fclose(fid);
        outputstring='Input data have 1 column - Select Grid dimension: 1D';
        M=[];
        return         
       end          
    else
        fclose(fid);
        outputstring='Wrong number of columns in the text file.';
        M=[];
        return
    end
    
    
    if ii==1
        if number_of_columns==3    
            M=data;
        else
            [X index]=unique(data(:,1));
            orig_data=data(:,2);
            M=[X zeros(size(X,1),1) orig_data(index)];
        end
        
    else
     
        
        if number_of_columns==1
            newcolumn = data(:,2);
        elseif number_of_columns==2
            newcolumn = data(:,2);
        else % no interpolation for grid data
            newcolumn=data(:,3);
        end
        try
            M=[M newcolumn];
        catch
            fclose(fid);
                outputstring='Wrong number of rows in the .dat files.';
            M=[];
        return
        end
        
    end
    
    %
    
end


% correct for possible inconsistency (a bottom layer that jumps above a top
% layer)
for ii=size(M,2):-1:4
    mask=M(:,ii)<M(:,ii-1);
    M(:,ii-1)=M(:,ii).*~mask+M(:,ii-1).*mask;
end


fclose(fid);

end


