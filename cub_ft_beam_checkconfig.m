function [ config, fieldtrip_path, subject_path, dataset_path, mrifile_name ] = cub_ft_beam_checkconfig( config )
%[ config, fieldtrip_path, subject_path, dataset_path, mrifile_name ] = cub_ft_beam_checkconfig( config )
%   LM 2-Feb-2016



%try to find the fieldtrip toolbox directory
if ~isfield(config,'fieldtrip_path') || isempty(config.fieldtrip_path),
    thepath = path;
    thepatharray = strread(thepath,'%s','delimiter',':');
    ftpath_idx = strfind(thepatharray,'fieldtrip');
    theemptycells = cellfun(@isempty,ftpath_idx);
    theftarray = thepatharray(find(~theemptycells));
    thelengtharray = cellfun(@length,theftarray);
    ftshort_idx = find(thelengtharray == min(thelengtharray));
    config.fieldtrip_path = theftarray{ftshort_idx};
end

%check fieldtrip path again
if isempty(config.fieldtrip_path),                                          error('Please specify the full path to the fieldtrip toolbox directory as "config.fieldtrip_path".'); end
fieldtrip_path = config.fieldtrip_path;

%check subject path
if ~isfield(config,'subject_path') || isempty(config.subject_path),        	error('Please specify the full path to the subject''s directory as "config.subject_path".'); end
if ~exist(config.subject_path, 'dir'),                                     	error(['Directory "' config.subject_path '" cannot be found. Please check "config.subject_path".']); end
subject_path = config.subject_path;

%check meg dataset name
if ~isfield(config,'dataset_path') || isempty(config.dataset_path),         error('Please specify the dataset folder path as "config.dataset_path"'); end
if ~exist(config.dataset_path, 'dir'),                                      fprintf(['*** Tips ***\nCheck .ds folder name. DO provide full path, not just name of .ds folder.\n']);
                                                                            error(['Directory "' config.dataset_path '" cannot be found.']); end
dataset_path = strtrim(config.dataset_path2);
%%%dataset_path = strtrim(config.dataset_path); % Origional

%check mri file name
if ~isfield(config,'mrifile_name') || isempty(config.mrifile_name),         error('Please specify the mri file name as "config.mrifile_name"'); end
if ~exist(fullfile(subject_path, config.mrifile_name), 'file'),             fprintf(['*** Tips ***\nCheck .mri file name. Move .mri file into the folder "' subject_path '". Do NOT provide full path, only .mri file name.\n']);
                                                                            error(['File "' config.mrifile_name '" cannot be found inside folder "' subject_path '".']); end
mrifile_name = config.mrifile_name;

end %function
