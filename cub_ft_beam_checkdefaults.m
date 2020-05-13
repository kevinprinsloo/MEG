function [ config ] = cub_ft_beam_checkdefaults( config )
%[ config ] = cub_ft_beam_checkdefaults( config )
%   LM 2-Feb-2016


%% definitions

subject_path    = config.subject_path;
mrifile_name    = config.mrifile_name;
beam_method     = config.method.beam;
dip_fitting     = config.method.dfit;
forw_string     = config.method.forw_str;
diff_metric     = config.method.diff;


%% template file

%default
if ~isfield(config,'template_filename') || isempty(config.template_filename),   config.template_filename = 'standard_sourcemodel3d5mm.mat'; end% 'standard_sourcemodel3d4mm.mat'; end

%check if the specified template exists
if ~exist(fullfile(config.fieldtrip_path, '/template/sourcemodel/', config.template_filename),'file'),
    error('Please specify the template file as "config.template_filename". The full path is required, unless the file is in the fieldtrip "/template/sourcemodel/" directory'); end


%% segmented mri, head model, leadfields

param_array = {'seg' 'hdm' 'ldf'};

for iP = 1:length(param_array)
    
    param = param_array{iP};
    
    %fix filename
    if strcmp(param,'hdm') || strcmp(param,'ldf'),                                          forw_string = ['_' config.method.forw_str];
    else                                                                                    forw_string = '';
    end
    
    
    %make the parameter field exist
    if ~isfield(config, param), config.(param) = []; end
    
    %parameter LOAD option NOT specified
    if ~isfield(config.(param),'load') || isempty(config.(param).load),
        
        %parameter LOAD file NOT specified - set default
        if ~isfield(config.(param),'load_file') || isempty(config.(param).load_file),       config.(param).load_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' param forw_string '.mat']));
                                                                                            config.(param).load = 1; %try to load
        end
        
        %default parameter LOAD file NOT found
        if ~exist(config.(param).load_file, 'file'),                                        config.(param).load = 0; %don't load
            if ~isfield(config.(param),'save') || isempty(config.(param).save),             config.(param).save = 1; end %save unless otherwise specified
        %default parameter LOAD file FOUND
        else                                                                                config.(param).load = 1; %load
            if ~isfield(config.(param),'save') || isempty(config.(param).save),             config.(param).save = 0; end %don't save unless otherwise specified
        end
                                                                                            fprintf(['\nSetting config.' param '.load = ' num2str(config.(param).load) ', by default.'])
                                                                                            fprintf(['\nSetting config.' param '.save = ' num2str(config.(param).save) ', by default.'])
    
    %parameter LOAD option = YES
    elseif config.(param).load == 1
        
        %parameter load file NOT specified - set default
        if ~isfield(config.(param),'load_file') || isempty(config.(param).load_file),     	config.(param).load_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' param forw_string '.mat'])); end
        
        %specified/default parameter file NOT found
        if ~exist(config.(param).load_file, 'file'),                                        error(['Please specify the full path to the "' param '" file as "config.' param '.load_file".']); end
        
        %default saving option
        if ~isfield(config.(param),'save') || isempty(config.(param).save),                 config.(param).save = 0; end %don't save
                                                                                            fprintf(['\nSetting config.' param '.save = ' num2str(config.(param).save) ', by default.'])
    
    %parameter LOAD option = NO
    elseif config.(param).load == 0
        
        %default saving option
        if ~isfield(config.(param),'save') || isempty(config.(param).save),                 config.(param).save = 1; end %do and save
                                                                                            fprintf(['\nSetting config.' param '.save = ' num2str(config.(param).save) ', by default.'])
    
    end
    
    %parameter saving option
    if ~isfield(config.(param),'save_file') || isempty(config.(param).save_file),           config.(param).save_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' param forw_string '.mat'])); end
    
    %parameter string to append
    if ~isfield(config.(param),'save_str') || isempty(config.(param).save_str),             config.(param).save_str = ''; end
    
    %re-define save name
    if ~isempty(config.(param).save_str),                                                   config.(param).save_file = strrep(config.(param).save_file, '.mat', ['_' config.(param).save_str '.mat']); end
    
    %get plot save name and clear unused load/save file names
    if config.(param).load == 0,                                                            config.(param).load_file = ''; end
    if config.(param).save == 0,
                                                                                            config.(param).plot_file = strrep(config.(param).save_file, '.mat', '.png');
                                                                                            config.(param).save_file = '';
    else
                                                                                            config.(param).plot_file = strrep(config.(param).save_file, '.mat', '.png');
    end
        
    
end %for iP
% fprintf('\n')

%fix string
if ~isempty(forw_string), forw_string = config.method.forw_str; end

%leadfiled normalisation
if ~isfield(config.ldf,'norm') || isempty(config.ldf.norm),                                	config.ldf.norm = 'no'; end %default = do NOT normalise leadfields

if      ~ischar(config.ldf.norm)
    if      config.ldf.norm == 0,                                                           config.ldf.norm = 'no';
    elseif  config.ldf.norm == 1,                                                           config.ldf.norm = 'yes';
    else                                                                                   	error('Please specify a valid normalisation option for "config.ldf.norm"');
    end
elseif  ischar(config.ldf.norm)
    if      strcmp(config.ldf.norm,'no'),
    elseif  strcmp(config.ldf.norm,'yes'),
    else                                                                                   	error('Please specify a valid normalisation option for "config.ldf.norm"');
    end
end


%% weights options

%make wts field exist
if ~isfield(config,'wts'),                                                  config.wts = []; end

%default config.wts options
if isempty(config.wts),
                                                                            config.wts.load = 0;
                                                                            config.wts.save = 0;
end

%default config.wts options
if ~isfield(config.wts,'load') || isempty(config.wts.load),                 config.wts.load = 0; end
if ~isfield(config.wts,'save') || isempty(config.wts.save),                 config.wts.save = 0; end

%check load options are valid
if      config.wts.load == 0,
elseif  config.wts.load == 1
else                                                                        error('Please specify a valid option for "config.wts.load".');
end

%load file options
if ~isfield(config.wts,'load_file') || isempty(config.wts.load_file),
    %default weights file, NOT loading
    if      config.wts.load == 0,                                           config.wts.load_file = '';
    %default weights fiel, to LOAD
    elseif  config.wts.load == 1,                                           config.wts.load_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' beam_method '_' dip_fitting '_' forw_string '_wts.mat']));
        %default file does not exist
        if  ~exist(config.wts.load_file,'file'),                            error('Please specify the full path to the weights file as "config.wts.load_file".'); end
    end
else
    %specified weights file, NOT loading
    if      config.wts.load == 0,                                           error('Please set "config.wts.load" to zero, or specify "config.wts.load_file" as empty.');
    %specified weights file, to LOAD
    elseif  config.wts.load == 1,
        %specified file does not exist
        if 	~exist(config.wts.load_file,'file'),                            error('Please specify the full path to the weights file as "config.wts.load_file".'); end
    end
end

%check save options are valid
if      config.wts.save == 0,
elseif  config.wts.save == 1
else                                                                        error('Please specify a valid option for "config.wts.save".');
end

fprintf(['\nSetting config.wts.load = ' num2str(config.wts.load) '.'])
fprintf(['\nSetting config.wts.save = ' num2str(config.wts.save) '.'])
fprintf('\n')

%load the weights config variable. We will try to set the covariance 
%frequencies to the same as the frequencies used to estimate the weights!!!
if config.wts.load == 1,
    try load(config.wts.load_file, 'wts_config'); end
    
    %covariance
    if ~exist('wts_config','var') || ~isfield(wts_config,'cov') || ~isfield(wts_config.cov,'bpfreq') || isempty(wts_config.cov.bpfreq),
                                                                            fprintf('\nCannot find the covariance band-pass frequencies that were used to estimate the beamformer filters. User should check!\n');
    else
                                                                            config.wts.bpfreq = wts_config.cov.bpfreq;
                                                                            fprintf('\nWeights covariance frequencies found: [%.1f %.1f] Hz.\n', wts_config.cov.bpfreq(1), wts_config.cov.bpfreq(2))
    end
    
    if ~exist('wts_config','var') || ~isfield(wts_config,'cov') || ~isfield(wts_config.cov,'bpfreq') || isempty(wts_config.cov.bpfreq),
                                                                            fprintf('\nCannot find the covariance band-pass frequencies that were used to estimate the beamformer filters. User should check!\n');
    else
                                                                            config.wts.bpfreq = wts_config.cov.bpfreq;
                                                                            fprintf('\nWeights covariance frequencies found: [%.1f %.1f] Hz.\n', wts_config.cov.bpfreq(1), wts_config.cov.bpfreq(2))
    end
    
    %CSD
    if ~exist('wts_config','var') || ~isfield(wts_config,'csd') || ~isfield(wts_config.csd,'freqarr') || isempty(wts_config.csd.freqarr),
                                                                            fprintf('\nCannot find the CSD frequencies that were used to estimate the beamformer filters. User should check!\n');
    else
                                                                            config.wts.freqarr = wts_config.csd.freqarr;
                                                                            fprintf('\nWeights CSD frequencies found: from %.1f to %.1f Hz.\n', wts_config.csd.freqarr(1), wts_config.csd.freqarr(end))
    end
    
    if ~exist('wts_config','var') || ~isfield(wts_config,'csd') || ~isfield(wts_config.csd,'freqarr') || isempty(wts_config.csd.freqarr),
                                                                            fprintf('\nCannot find the CSD frequencies that were used to estimate the beamformer filters. User should check!\n');
    else
                                                                            config.wts.freqarr = wts_config.csd.freqarr;
                                                                            fprintf('\nWeights CSD frequencies found: from %.1f to %.1f Hz.\n', wts_config.csd.freqarr(1), wts_config.csd.freqarr(end))
    end
end

%save file options
if ~isfield(config.wts,'save_file') || isempty(config.wts.save_file),
    %default weights file, NOT saving
    if      config.wts.save == 0,                                           config.wts.save_file = '';
    %default weights fiel, to SAVE
    elseif  config.wts.save == 1,                                           config.wts.save_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' beam_method '_' dip_fitting '_' forw_string '_wts.mat']));
    end
else
    %specified weights file, NOT saving
    if      config.wts.save == 0,                                           error('Please set "config.wts.save" to zero, or specify "config.wts.save_file" as empty.'); end
end

%parameter string to append
if ~isfield(config.wts,'save_str') || isempty(config.wts.save_str),         config.wts.save_str = ''; end

%re-define save name
if ~isempty(config.wts.save_str) && ~isempty(config.wts.save_file),         config.wts.save_file = strrep(config.wts.save_file, '.mat', ['_' config.wts.save_str '.mat']); end


%% source options

%distinguish between 'statistics' and 'percentage'
if      strcmp(diff_metric,'statistics'),                                   diff_string = 'stat';
elseif  strcmp(diff_metric,'percentage'),                                   diff_string = 'perc';
else                                                                        error(['Unknown diff_metric "' diff_metric '"']);
end

if ~isfield(config,'src'),                                                  config.src = []; end
if ~isfield(config.src,'save') || isempty(config.src.save),                 config.src.save = 0;
                                                                            config.src.save_file = '';
else
    if config.src.save == 0
        if ~isfield(config.src,'save_file') || isempty(config.src.save_file),
                                                                            config.src.save_file = '';
        end
    elseif config.src.save == 1
        if ~isfield(config.src,'save_file') || isempty(config.src.save_file),
                                                                            config.src.save_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' beam_method '_' dip_fitting '_' forw_string '_' diff_string '-src.mat'])); end
    else
                                                                            error('Please specify a valid option for "config.src.save".');
    end
end
if ~isfield(config.src,'save_str') || isempty(config.src.save_str),         config.src.save_str = ''; end

%the following is now a bit useless, so...
if ~isempty(config.src.save_file) && ~isempty(config.src.save_str),         config.src.save_file = strrep(config.src.save_file, '.mat', ['_' config.src.save_str '.mat']); end
if      ~isempty(config.src.save_str),                                      config.src.plot_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' beam_method '_' dip_fitting '_' forw_string '_' diff_string '-src_' config.src.save_str '.png']));
elseif  isempty(config.src.save_str),                                       config.src.plot_file = fullfile(subject_path, strrep(mrifile_name, '.mri', ['_' beam_method '_' dip_fitting '_' forw_string '_' diff_string '-src.png']));
end


%% plot options

if ~isfield(config,'plot'), config.plot = []; end

param_array = {'seg' 'hdm' 'src'};

for iP = 1:length(param_array)
    
    param = param_array{iP};
    
    if ~isfield(config.plot, param) || isempty(config.plot.(param)),        config.plot.(param) = 0; end
    if      config.plot.(param) == 0,
    elseif  config.plot.(param) == 1,
    else                                                                    error(['Please specify a valid option for "config.plot.' param '".']);
    end
    
end

if ~isfield(config.plot, 'save') || isempty(config.plot.save),              config.plot.save = 1; end %default to 1, as this is only called if something is plotted first.
if      config.plot.save == 0,
elseif  config.plot.save == 1,
else                                                                        error(['Please specify a valid option for "config.plot.save".']);
end


end %function
