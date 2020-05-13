function [ config ] = cub_ft_beam_checkdata( config )
%[ config ] = cub_ft_beam_checkdata( config )
%   LM 3-Feb-2016


%% data options

%data
if ~isfield(config,'data') || isempty(config.data),                         error('Please specify the config.data field.'); end
if ~isfield(config.data,'markertype') || isempty(config.data.markertype),   config.data.markertype = 'UPPT002'; end %default is parallel port 2
if ischar(config.data.markertype) && strcmp(config.data.markertype,'UPPT002'),
    if ~isfield(config.data,'markervalue') || isempty(config.data.markervalue),
                                                                            error('Please specify the config.data.markervalue field.');
    end
end
if ~isfield(config.data,'timewin') || isempty(config.data.timewin),         error('Please specify the config.data.timewin field.');
elseif length(config.data.timewin) ~= 2,                                    error('Please check the config.data.timewin field. Specify a time-window in seconds, relative to the marker, for the whole data epoch.');
end
%baseline
if ~isfield(config,'bsln') || isempty(config.bsln),                         error('Please specify the config.bsln field.'); end
if ~isfield(config.bsln,'timewin') || isempty(config.bsln.timewin),         error('Please specify the config.bsln.timewin field.');
elseif length(config.bsln.timewin) ~= 2,                                    error('Please check the config.bsln.timewin field. Specify a time-window in seconds, relative to the marker, for the baseline epoch.');
end
%active
if ~isfield(config,'actv') || isempty(config.actv),                         error('Please specify the config.actv field.'); end
if ~isfield(config.actv,'timewin') || isempty(config.actv.timewin),         error('Please specify the config.actv.timewin field.');
elseif length(config.actv.timewin) ~= 2,                                    error('Please check the config.actv.timewin field. Specify a time-window in seconds, relative to the marker, for the active epoch.');
end


%% covariance and CSD

%time-domain beamformers
if strcmp(config.method.beam,'sam') || strcmp(config.method.beam,'lcmv')
    
    if ~isfield(config,'cov') || isempty(config.cov),                           error('Please specify the config.cov field.'); end
    if ~isfield(config.cov,'timewin') || isempty(config.cov.timewin),           error('Please specify the config.cov.timewin field.');
    elseif length(config.cov.timewin) ~= 2,                                     error('Please check the config.cov.timewin field. Specify a time-window in seconds, relative to the marker, for the covariance.');
    end
    if ~isfield(config.cov,'bpfreq') || isempty(config.cov.bpfreq),
        %check if pre-computed filters are used, and if so try to find the freqs
        if config.wts.load == 1 && isfield(config.wts,'bpfreq'),
            config.cov.bpfreq = config.wts.bpfreq;
            fprintf('\nUsing the same covariance frequencies as those used to estimate the weights\n');
        else
            error('Please specify the config.cov.bpfreq field.');
        end
    else
        %check if pre-computed filters are used, and force the freqs
        if config.wts.load == 1 && isfield(config.wts,'bpfreq'),
            config.cov.bpfreq = config.wts.bpfreq;
            fprintf('\nForcing the covariance frequencies to be the same as those used to estimate the weights\n');
        end
    end
    if length(config.cov.bpfreq) ~= 2,                                          error('Please check the config.cov.bpfreq field. Specify the band-pass frequencies in Hz, for the covariance.');
    end
    
%freq-domain beamformers
elseif strcmp(config.method.beam,'dics')
    
    if ~isfield(config,'csd') || isempty(config.csd),                           error('Please specify the config.csd field.'); end
    if ~isfield(config.csd,'timewin') || isempty(config.csd.timewin),           error('Please specify the config.csd.timewin field.');
    elseif length(config.csd.timewin) ~= 2,                                     error('Please check the config.csd.timewin field. Specify a time-window in seconds, relative to the marker, for the CSD.');
    end
    if ~isfield(config.csd,'freqarr') || isempty(config.csd.freqarr),
        %check if pre-computed filters are used, and if so try to find the freqs
        if config.wts.load == 1 && isfield(config.wts,'freqarr'),
            config.csd.freqarr = config.wts.freqarr;
            fprintf('\nUsing the same CSD frequencies as those used to estimate the weights\n');
        else
            error('Please specify the config.csd.freqarr field.');
        end
    else
        %check if pre-computed filters are used, and force the freqs
        if config.wts.load == 1 && isfield(config.wts,'freqarr'),
            config.csd.freqarr = config.wts.freqarr;
            fprintf('\nForcing the CSD frequencies to be the same as those used to estimate the weights\n');
        end
    end
end


%% BAD trials

if ~isfield(config.data,'excludebad') || isempty(config.data.excludebad),
                                                                            config.data.excludebad = 0;
                                                                            fprintf('\n***** BAD trials will NOT be excluded. *****\n');
else
    if      config.data.excludebad == 0,
                                                                            fprintf('\n***** BAD trials will NOT be excluded. *****\n');
    elseif  config.data.excludebad == 1,
                                                                            fprintf('\n***** BAD trials will be EXCLUDED. *****\n');
    else
                                                                            error('Please specify a valid option for config.data.excludebad.');
    end
end


%% always SAVE data, by default

%FIXME
config.data.save = 1;
config.data.save_file = fullfile(config.subject_path, strrep(config.mrifile_name, '.mri', ['_' config.method.beam '_data.mat']));


end

