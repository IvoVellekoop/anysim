%
% setPath()
% setPath(setNewPathFlag)
%
% Clears the Matlab path and sets it for this project.
% Do not move this function to another folder!
%
% setNewPathFlag: Optional, set to false to reset the original path. (default: true)
%
function setPath(setNewPathFlag)
  persistent user_path_backup;
  
  if nargin < 1
    setNewPathFlag = true;
  end
  
  if setNewPathFlag
    updateOnly = exist('user_path_backup', 'var');
    user_path_backup = path();

    folder = fileparts(mfilename('fullpath'));

    if ~updateOnly
      sprintf('Clearing Matlab path and setting it for this project to %s...\n', folder);
      restoredefaultpath();
    else
      sprintf('Updating Matlab path and setting it for this project to %s...\n', folder);
    end

    addpath(genpath(folder));

    disp('Done setting path. Enter ''resetPath'' to return to the previously set path.');
  else
    if ~isempty(user_path_backup)
      disp('Resetting Matlab path to the previously set path...');
  
      path(user_path_backup);
    else
      disp('No previously set Matlab path known...');
    end

    user_path_backup = [];
  end
end