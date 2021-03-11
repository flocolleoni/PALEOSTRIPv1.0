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
% TEST IF FILES TO SAVE DATA ALREADY EXISTS
%
% Ask user if overwritting existing files
%-----------------------------------------------------------------------------------

function [fid]=check_save_to_file(fpath,destinationfolder,filename) 
 overwriteFile = false; % Default to overwriting or creating new file.

 fexist= exist(fpath);
 if (fexist > 0)
      % Ask user if they want to overwrite the file.
      filestring="already exists - Overwrite it?";
      promptMessage = sprintf(fpath,filestring);%sprintf('This file already exists:\n%s\nDo you want to overwrite it?', fpath);
      titleBarCaption = 'Overwrite it?';
      buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
      if strcmpi(buttonText, 'Yes')
         overwriteFile = true; % User want to overwrite.
      end

     if overwriteFile
      % File does not exist yet, or the user wants to overwrite an existing file.
      delete(fullfile(destinationfolder,filename));

       fid=fopen(fullfile(destinationfolder,filename),'wt');
        if fid<0
          errordlg(['Error opening the file' fullfile(destinationfolder,filename) ' for writing.'],'', 'modal')
          return
        end
     else
       [filename,destinationfolder] = uiputfile({'*.dat','Data file (*.dat)';'*.txt','Text file (*.txt)';'*.*','All files (*.*)'},'Save File Name',filename);

       fid=fopen(fullfile(destinationfolder,filename),'wt');
        if fid<0
          errordlg(['Error opening the file' fullfile(destinationfolder,filename) ' for writing.'],'', 'modal')
          return
        end
     end 

 elseif (fexist==0)  
 %OPEN FILE
  fid=fopen(fullfile(destinationfolder,filename),'wt');
  if fid<0
    errordlg(['Error opening the file' fullfile(destinationfolder,filename) ' for writing.'],'', 'modal')
    return
  end         
 end
 
end