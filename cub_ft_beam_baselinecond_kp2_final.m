function [ src_diff, src_peak, data_virt, config ] = cub_ft_beam_baselinecond_kp2_final( config )
%[ src_diff, src_peak, data_virt, config ] = cub_ft_beam_baselinecond( config )
%   
%   Input:
%   "config",           configuration structure with various fields. 
%                       
%                       The mandatory config fields are: 
%                       - 'subject_path' (full path to the analysis folder)
%                       - 'dataset_path' (full path to the analysis dataset)
%                       - 'mrifile_name' (name of the coregistered .mri file 
%                           NOTE: this will also be used as prefix 
%                           to the analysis output filenames)
%                       - 'data', with struct subfields 'markertype',
%                           'markervalue', 'timewin', 'excludebad', ...
%                       - 'bsln', with struct subfields 'timewin', ...
%                       - 'actv', with struct subfields 'timewin', ...
%                       - 'cov', with struct subfields 'timewin', 
%                           'bpfreq', ...
%                       - 'csd', with struct subfields 'timewin',
%                           'freqarr', 'freqsmo', ...
%                       
%                       Other optional config fields are:
%                       - 'method', struct with subfields...
%                       - 'src', struct with subfields...
%                       - 'ldf', struct with subfields...
%                       - 'plot', struct with subfields...
%                       - 'template_filename', string...
%                       - 'seg', struct with subfields...
%                       - 'hdm', struct with subfields...
%   
%   Output:
%   "src_diff",         ft-structure of the source solution, computer as: 
%                       - paired-t values (uncorrected), when
%                           config.method.diff = 'statistics' [default]
%                       - percentage change from baseline, when
%                           config.method.diff = 'percentage'
%   "src_peak",         structure with fields: 
%                       'max', either t-value or percentage value of the 
%                       peak voxel (across all voxels inside the brain).
%                       'idx', index to the peak voxel relative to the
%                       source model.
%                       'mni', mni coordinates (in mm) of the peak voxel.
%   "data_virt",        ft-like structure of the virtual sensor
%                       reconstruction for the peak voxel location. The
%                       whole trials are reconstructed by projecting the
%                       sensor level data through the common filter. The 
%                       data is epoched according to config.data.timewin,
%                       around the marker specified by config.data.markertype
%                       and config.data.markervalue.
%   "config",           the input config, with extra fields set by default, 
%                       is returned for bookkeeping purposes.

% First implementation by LM, 15-Feb-2016

% To-do: produce an error when the loaded template filename does not match
% the input template filename


%% config

[config, fieldtrip_path, subject_path, dataset_path, mrifile_name] = cub_ft_beam_checkconfig(config);

%display settings
fprintf(['\nUsing fieldtrip version: ' fieldtrip_path '\n'])
fprintf(['Running analysis in subject''s folder: ' subject_path '\n'])
fprintf(['Analysing dataset: ' dataset_path '\n'])
fprintf(['Using mri file: ' mrifile_name '\n'])

%% config.method

[config, beam_method, dip_fitting, fixedori, forw_method, diff_metric] = cub_ft_beam_checkmethod(config);

%print out some stuff
fprintf('\n*** Using: %s "%s" with %s ***\n', beam_method, dip_fitting, forw_method)

%% more config options

%check and set defaults
config = cub_ft_beam_checkdefaults(config);

%template
template_filename = config.template_filename;
%segmented mri
seg_load        = 0;
seg_loadfile    = config.seg.load_file;
seg_save        = config.seg.save;
seg_savefile    = config.seg.save_file;
%seg_savefile    = [seg_savefile(1:end-5) '_seg.mat']
seg_plotfile    = config.seg.plot_file;
%headmodel
hdm_load        = 0;
hdm_loadfile    = config.hdm.load_file;
hdm_save        = config.hdm.save;
hdm_savefile    = config.hdm.save_file;
hdm_plotfile    = config.hdm.plot_file;
%leadfield
ldf_load        = 0;
ldf_loadfile    = config.ldf.load_file;
ldf_save        = config.ldf.save;
ldf_savefile    = config.ldf.save_file;
ldf_norm        = config.ldf.norm;
%weights
wts_load        = config.wts.load;
wts_loadfile    = config.wts.load_file;
wts_save        = config.wts.save;
wts_savefile    = config.wts.save_file;
%source
src_save        = config.src.save;
src_savefile    = config.src.save_file;
src_plotfile   = config.src.plot_file;
%src_plotfile = config.src.sub;
%plot
plot_seg    = config.plot.seg;
plot_hdm    = 1;
plot_src    = config.plot.src;
plot_save   = config.plot.save;

%% data and covariance options

config = cub_ft_beam_checkdata(config);

%print out some stuff
fprintf('\n')
fprintf('The data will be epoched from %.1f to %.1f s, relative to marker(s): %s.\n', config.data.timewin(1), config.data.timewin(2), num2str(config.data.markervalue))
fprintf('The beamformer will contrast active [%.1f %.1f] vs. baseline [%.1f %.1f].\n', config.actv.timewin(1), config.actv.timewin(2), config.bsln.timewin(1), config.bsln.timewin(2))

if strcmp(beam_method,'sam') || strcmp(beam_method,'lcmv')
    fprintf('The covariance will be calculated from %.1f to %.1f s, for frequencies between %d and %d Hz.\n', config.cov.timewin(1), config.cov.timewin(2), config.cov.bpfreq(1), config.cov.bpfreq(2))
elseif strcmp(beam_method,'dics')
    fprintf('The CSD will be calculated from %.1f to %.1f s, for frequencies: %s Hz.\nSmoothing +/- %s Hz.\n', config.csd.timewin(1), config.csd.timewin(2), num2str(config.csd.freqarr), num2str(config.csd.freqsmo))
end

%data will be saved by default
dat_save = config.data.save;
dat_savefile = config.data.save_file;

%definitions
markertype      = config.data.markertype;
markervalue     = config.data.markervalue;
markerpre       = config.data.timewin(1)*-1;
markerpos       = config.data.timewin(2);
bsln_timewin    = config.bsln.timewin;
actv_timewin    = config.actv.timewin;

%covariance and CSD
if strcmp(beam_method,'sam') || strcmp(beam_method,'lcmv')
    cov_timewin     = config.cov.timewin;
    cov_bpfreq      = config.cov.bpfreq;
elseif strcmp(beam_method,'dics')
    csd_timewin     = config.csd.timewin;
    csd_freqarr     = config.csd.freqarr;
    csd_freqsmo     = config.csd.freqsmo;
end


%% mri stuff

%read in mri file
load(fullfile(subject_path, mrifile_name),'mri'); %% problem
mri = ft_convert_units(mri, 'cm');

%do mri segmentation
if seg_load == 1
    fprintf('\nLoading existing segmented mri from file %s\n', seg_loadfile)
    mri_seg = load(seg_loadfile);
else
    cfg = [];
    cfg.output = {'brain'; 'skull'; 'scalp'};
    mri_seg = ft_volumesegment(cfg, mri);
    %add anatomical information to the segmentation
    mri_seg.transform = mri.transform;
    mri_seg.anatomy   = mri.anatomy;
end

%check segmentation
if plot_seg == 1
    cfg              = [];
    cfg.funparameter = 'brain';% 'gray';% 
    ft_sourceplot(cfg,mri_seg);
    drawnow;
    
    %quick save
    if plot_save == 1
        saveas(gcf, strrep(seg_plotfile, '.mat', '.png'))
        % close
    end
end

%save segmented mri, if requested
if seg_save == 1
    fprintf('\nSaving segmented mri file as %s\n', seg_savefile)
    save(seg_savefile, '-struct', 'mri_seg')
end


%% head model

%get grad structure from header file
hdr = ft_read_header(dataset_path);
grad = hdr.grad;
grad = ft_convert_units(grad, 'cm'); %make sure units match...

%load headmodel, if requested
if hdm_load == 1
    
    fprintf('\nLoading existing volume conduction model (headmodel) from file %s\n', hdm_loadfile)
    hdm = load(hdm_loadfile);
    
else %prepare headmodel
    
    if strcmp(forw_method,'singleshell')
        
        %singleshell in native space
        cfg = [];
        cfg.method = forw_method;% 'singleshell';
        hdm = ft_prepare_headmodel(cfg, mri_seg);
        hdm = ft_convert_units(hdm, 'cm');
        
    elseif strcmp(forw_method,'localspheres')
        
        %localspheres in native space
        cfg = [];
        cfg.method = forw_method;% 'localspheres';
        cfg.grad = grad;
        cfg.feedback = 'no';% 'yes';
        hdm = ft_prepare_headmodel(cfg, mri_seg);
        hdm = ft_convert_units(hdm, 'cm');
        
    else
        error(['Sorry, can''t deal with the ' forw_method ' forward model, yet...'])
    end
    
end

%make space
clear mri_seg

%save headmodel
if hdm_save == 1
    fprintf('\nSaving volume conduction model (headmodel) as %s\n', hdm_savefile)
    save(hdm_savefile, '-struct', 'hdm')
end


%% source model and leadfield

%first load the template sourcemodel
fprintf('\nLoading template sourcemodel from file %s\n', template_filename)
load(template_filename)
if ~exist('template_sourcemodel','var')
    template_sourcemodel = sourcemodel; %sourcemodel is Fieldtrip's varname
    clear sourcemodel
end

%then load sourcemodel and leadfiled if requested, or compute them otherwise
if ldf_load == 1
    
    fprintf('\nLoading pre-computed leadfield from file %s\n', ldf_loadfile)
    load(ldf_loadfile, 'leadfield', 'sourcemodel', 'template_filename');
    fprintf('NOTE: "sourcemodel" was normalised to template using file %s.\nUser should check!\n', template_filename)
    
else
    
    %convert to MNI template space
    cfg                 = [];
    cfg.grid.warpmni    = 'yes';
    cfg.grid.template   = template_sourcemodel;
    cfg.grid.nonlinear  = 'yes';
    cfg.mri             = mri;
    cfg.inwardshift     = -1.5;
    
    sourcemodel = ft_prepare_sourcemodel(cfg);
    
    %check headmodel
    if plot_hdm == 1
        figure; hold on;
        ft_plot_vol(hdm);%, 'edgecolor', 'none')
        %alpha 0.4
        ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
        ft_plot_sens(grad);
        hold off;
        drawnow;
        
        %quick save
        if plot_save == 1
            saveas(gcf, strrep(hdm_plotfile, '.mat', '.png'))
            % close
        end
    end
    
    %compute the leadfield matrices
    %note: we now define the grid resolution based on the pre-computed
    %sourcemodel, which is the irregular grid warped onto the template...
    cfg                 = [];
    cfg.channel         = {'MEG'};% chan;
    cfg.grid            = sourcemodel;
    cfg.vol             = hdm;
    cfg.grad            = grad;
    cfg.normalize       = ldf_norm; %'yes' would remove depth bias (Q in eq. 27 of van Veen et al, 1997)
    
    %leadfield = ft_prepare_leadfield(cfg);
    
    %save sourcemodel and leadfield variables, along with the template used
    if ldf_save == 1
        fprintf('\nSaving sourcemodel and leadfield as %s\n', ldf_savefile)
        save(ldf_savefile, 'sourcemodel', 'leadfield', 'template_filename')
    end
end

%make space - we won't need this anymore, as now go into template space
clear mri %could clear also "sourcemodel"?


%% data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in epoched data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.dataset = dataset_path;
cfg.trialdef.eventtype  = markertype;
cfg.trialdef.eventvalue = markervalue;
cfg.trialdef.prestim    = markerpre;
cfg.trialdef.poststim   = markerpos;

cfg_data = ft_definetrial(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess and band-pass filter the data (covariance frequencies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time-domain preprocessing
if strcmp(beam_method,'sam') || strcmp(beam_method,'lcmv')
    
    %FIXME: lower frequencies may require padding before band-pass filtering ??
    
    %preprocessing options
    cfg_data.channel    = {'all'};
    cfg_data.precision  = 'single';% 'double';% 
    cfg_data.detrend     = 'yes'; %Glasgow replaced demean with detrend
    cfg_data.bpfilter   = 'yes';
    cfg_data.bpfreq     = cov_bpfreq;
    cfg_data.bpfilttype = 'but';
    cfg_data.bpfiltdir  = 'twopass';
    
    data_prep_all = ft_preprocessing(cfg_data);
    
%freq-domain preprocessing
elseif strcmp(beam_method,'dics')
    
    %preprocessing options
    cfg_data.channel    = {'all'};
    cfg_data.precision  = 'single';% 'double';% 
    cfg_data.detrend     = 'yes'; %Glasgow replaced demean with detrend
    cfg_data.hpfilter   = 'yes';
    cfg_data.hpfreq     = 1; %high-pass at 1 Hz shouldn't hurt 
    
    data_prep_all = ft_preprocessing(cfg_data);
    
end

%downsample to XXX Hz - speed things up!
cfg = [];
cfg.resamplefs  = 300;
cfg.detrend     = 'no';
data_prep_res = ft_resampledata(cfg, data_prep_all);


cfg=[];
data_prep_res=ft_denoise_pca(cfg,data_prep_res);

cfg=[];
cfg.channel={'MEG', '-A102' , '-A107', '-A146', '-A40'}; %'-A102',
data_prep = ft_preprocessing(cfg,data_prep_res);
%Use this if using ft_rejectvisual
%data_prep_res = ft_preprocessing(cfg,data_prep_res);

%cfg=[];
%cfg.method='summary';
%data_prep=ft_rejectvisual(cfg,data_prep_res);
clear data_prep_all data_prep_res cfg_data class_*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% work out BAD trials and exclude them (if requested)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %read .cls ClassFile
% class_file = fullfile(dataset_path, 'ClassFile.cls');
% if ~exist(class_file,'file'), error(['ClassFile.cls does cannot be found inside dataset folder ' dataset_path]); end
% class_info = readClassFile(class_file);
% if max(size(class_info))~=1 || ~strcmp(class_info.Name,'BAD'), error(['Cannot corretly read BAD trials. Please check file ' class_file]); end
% 
% %get BAD trials index
% trials_bad = class_info.trial;
% trials_all = 1:size(data_prep_res.trialinfo,1);
% 
% %exclude BAD trials
% cfg = [];
% cfg.trials = trials_all(~ismember(trials_all,trials_bad));
% 
% if ~isempty(cfg.trials) && config.data.excludebad == 1
%     fprintf('\nBAD trials will be removed (%d out of %d).\n', length(trials_bad), length(trials_all))
%     data_prep = ft_preprocessing(cfg, data_prep_res);
% else
%     data_prep = data_prep_res;
% end
% 
% %make space, avoid mistakes
% clear data_prep_all data_prep_res cfg_data class_*
% 

%% beamformer weights

%NOTE: I've removed the option to load the source solution, as this was 
%making things too complicated... and it didn't turn out to be very useful
%LM 4-Feb-2016

if wts_load == 1
    
    fprintf('\nLoading beamformer weights from file %s', wts_loadfile)
    fprintf('\nSkipping computation of common source solution')
    fprintf('\nUser should check "wts_cfg" and "wts_config" variables...\n')
    load(wts_loadfile, 'wts');
    
    if src_save == 1
        fprintf('\n***** Warning: *****')
        fprintf('\nSource common data will NOT be saved.\n')
    end
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute global covariance, i.e. using the whole time-window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %time-domain: covariance
    if strcmp(beam_method,'sam') || strcmp(beam_method,'lcmv')
        
        %covariance matrix
        cfg = [];
        %%%cfg.channel = {'MEG'}; % Orig script
        cfg.channel={'MEG', '-A102' , '-A107', '-A146', '-A40'}; %'-A102',
        cfg.removemean = 'no';
        cfg.covariance = 'yes';
        cfg.covariancewindow = cov_timewin;% 'all';
        
        tlck = ft_timelockanalysis(cfg, data_prep);
        data_beam = tlck;
        
    %freq-domain: cross-spectral density
    elseif strcmp(beam_method,'dics')
        
        %csd matrix
        cfg = [];
        %%%cfg.channel     = {'MEG'}; % orig script
        cfg.channel={'MEG', '-A102' , '-A107', '-A146', '-A40'}; %'-A102',
        cfg.method      = 'mtmfft';
        cfg.output      = 'powandcsd';
        cfg.taper       = 'dpss';
        cfg.toi         = csd_timewin;
        cfg.foi         = csd_freqarr;
        cfg.tapsmofrq   = csd_freqsmo;
        
        freq = ft_freqanalysis(cfg, data_prep);
        data_beam = freq;
        clear freq
        
        %NOTE: it would be more efficient to first compute frequency 
        %separately for the two conditions, and then combine for the common
        %filter... but we keep it this way instead... for consistency with 
        %the time-domain beamformer, the covariance matrix is computed on a 
        %time range that can exceed the individual time ranges of bsln and
        %actv... this way, we can then use the weights to reconstruct the
        %time-series also *outside* the bsln+actv window!
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute common filters, i.e. using the whole time-window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %common spatial filter
    cfg = [];
    cfg.method      = beam_method;
    cfg.grid        = sourcemodel;
    cfg.headmodel   = hdm;
    cfg.grad        = grad;
    cfg.channel=data_beam.label;
    leadfield = ft_prepare_leadfield(cfg);
    leadfield = rmfield(leadfield, 'cfg'); % Kev added from JG
    cfg.grid         = leadfield;

    %beamformer options
    cfg.(beam_method).fixedori      = fixedori;
    cfg.(beam_method).keepfilter    = 'yes';
    cfg.(beam_method).projectnoise  = 'yes';
    cfg.(beam_method).lambda        = '5%'; %regularization parameter
    
    %DICS-only options
    if strcmp(beam_method,'dics')
    cfg.frequency                   = mean(csd_freqarr);
    cfg.(beam_method).realfilter    = 'yes';
    end
    
    src = ft_sourceanalysis(cfg, data_beam);
    
    %replace the source coordinates with the template ones
    src.pos = template_sourcemodel.pos;
    src.dim = template_sourcemodel.dim;
    
    %save source
%     if src_save == 1
%         fprintf('\nSaving source common data as %s\n', src_savefile)
%         save(src_savefile, '-v7.3', 'src')
%     end
%     
    %store the common spatial filter (aka beamformer weights)
    if strcmp(beam_method, 'sam')
        wts = src.avg.filter';
    else
        wts = src.avg.filter;
    end
    wts_cfg = cfg;
    wts_config = config;
    
    %take a look at the norm of the weights
%     insidx = find(src.inside);
%     normfilt = nan(length(insidx),1);
%     for i=1:3294
%         normfilt(i) = norm(cell2mat(src.avg.filter(insidx(i))));
%     end
%     figure, plot(normfilt)
%     close
    
%     %save weights
%     if wts_save == 1
%         fprintf('\nSaving beamformer weights as %s\n', wts_savefile)
%         save(wts_savefile, '-v7.3', 'wts', 'wts_cfg', 'wts_config')
%     end
    
end %if wts_load


%% plot, just for fun, as there's no contrast here yet!

% templatefile = fullfile(fieldtrip_path, '/template/anatomy/single_subj_T1.nii');
% template_mri = ft_read_mri(templatefile);
% template_mri = ft_volumereslice([], template_mri);
% 
% cfg              = [];
% cfg.parameter    = 'avg.pow';
% cfg.downsample=4;
% cfg.interpmethod = 'nearest';
% 
% src_plot = ft_sourceinterpolate(cfg, src, template_mri);
% src_plot.coordsys = 'mni';
% src_plot = ft_convert_units(src_plot, 'mm'); %make sure it's mm
% 
% %plot orthogonally
% cfg                     = [];
% cfg.method              = 'slice';
% cfg.funparameter        = 'pow';
% cfg.locationcoordinates = 'head';
% cfg.funcolormap         = 'jet';
% %cfg.atlas               = fullfile(fieldtrip_path, '/template/atlas/aal/ROI_MNI_V4.nii');
% cfg.downsample          = 4;
% 
%  ft_sourceplot(cfg,src_plot);
% 

%% baseline and active windows

%we will then beam these using the common filter computed above
cfg = [];
cfg.toilim = bsln_timewin;
data_prep_bsln = ft_redefinetrial(cfg, data_prep);
cfg = [];
cfg.toilim = actv_timewin;
data_prep_actv = ft_redefinetrial(cfg, data_prep);

%time-domain: covariance
if strcmp(beam_method,'sam') || strcmp(beam_method,'lcmv')
    
    %timelock analysis of bandpass-filtered baseline and active windows
    cfg             = [];
    cfg.channel     = {'MEG'};
    cfg.keeptrials  = 'yes';
    cfg.covariance  = 'yes'; %this is crucial to actually get the power estimates
    tlck_bsln = ft_timelockanalysis(cfg, data_prep_bsln);
    tlck_actv = ft_timelockanalysis(cfg, data_prep_actv);
    
    data_beam_bsln = tlck_bsln;
    data_beam_actv = tlck_actv;
    
%freq-domain: cross-spectral density
elseif strcmp(beam_method,'dics')
    
    %frequency analysis of baseline and active windows
    cfg = [];
    cfg.channel     = {'MEG'};
    cfg.keeptrials  = 'yes';
    cfg.method      = 'mtmfft';
    cfg.output      = 'powandcsd';
    cfg.taper       = 'dpss';
    cfg.foi         = csd_freqarr;
    cfg.tapsmofrq   = csd_freqsmo;
    
    freq_bsln = ft_freqanalysis(cfg, data_prep_bsln);
    freq_actv = ft_freqanalysis(cfg, data_prep_actv);
    
    data_beam_bsln = freq_bsln;
    data_beam_actv = freq_actv;
    clear freq_bsln freq_actv
    
end

%%%%%%%%%%%%%
% SAVE data %
%%%%%%%%%%%%%

%save baseline and active data 'to be beamed' (either tlck or freq)
% if dat_save == 1
%     fprintf('\nSaving baseline and active ''data_beam'' into %s\n', dat_savefile)
%     save(dat_savefile, '-v7.3', 'data_beam_bsln', 'data_beam_actv')
% end


%% beamforming

%NOTE: I've removed the option to load the source solution, see above...
%LM 4-Feb-2016

%get source-power estimates of bsln and actv separately
cfg = [];
cfg.method      = beam_method;
cfg.grid        = leadfield;
cfg.grid.filter = wts;
cfg.headmodel   = hdm;
cfg.grad        = grad;

if strcmp(diff_metric,'statistics')
    cfg.rawtrial    = 'yes';
    cfg.keeptrials  = 'no';
else
    %only remember trialinfo... otherwise use in conjunction with 'rawtrial'...
    cfg.keeptrials  = 'yes';
end

%beamformer options
cfg.(beam_method).fixedori      = fixedori;
cfg.(beam_method).keepfilter    = 'yes';
cfg.(beam_method).projectnoise  = 'yes';
cfg.(beam_method).lambda        = '5%'; %regularization parameter

%DICS-only options
if strcmp(beam_method,'dics')
cfg.frequency                   = mean(csd_freqarr);
cfg.(beam_method).realfilter    = 'yes';
end

src_bsln = ft_sourceanalysis(cfg, data_beam_bsln);
src_actv = ft_sourceanalysis(cfg, data_beam_actv);

%update the pos field...
src_bsln.pos = template_sourcemodel.pos;
src_bsln.dim = template_sourcemodel.dim;
src_actv.pos = template_sourcemodel.pos;
src_actv.dim = template_sourcemodel.dim;

%save
% if src_save == 1
%     
%     if wts_load == 1
%         %do not append, overwrite existing file
%         fprintf('\nSource common data was NOT saved into %s', src_savefile)
%         fprintf('\nSaving source baseline data as %s\n', src_savefile)
%         save(src_savefile, '-v7.3', 'src_bsln')
%     else
%         %append to existing file
%         fprintf('\nSaving source baseline data as %s\n', src_savefile)
%         save(src_savefile, '-v7.3', '-append', 'src_bsln')
%     end
%     
%     %append in both scenarios
%     fprintf('\nSaving source active data as %s\n', src_savefile)
%     save(src_savefile, '-v7.3', '-append', 'src_actv')
%     
% end %if src_save
% 

%% compute difference in power as...

%...percentage change from baseline
if      strcmp(diff_metric,'percentage')
    src_diff = cub_ft_beam_perc([], src_actv, src_bsln, template_sourcemodel);
%...paired-t statistics
elseif  strcmp(diff_metric,'statistics')
    src_diff = cub_ft_beam_stat([], src_actv, src_bsln, template_sourcemodel);
end

%make space
clear src_bsln src_actv


%% peak coordinate

%FIXME: it would be nice to automatically force the peak to be occipital...

%get peak voxel
if      strcmp(diff_metric,'percentage')
    [vox_max, vox_idx] = max(src_diff.pow);
elseif  strcmp(diff_metric,'statistics')
    [vox_max, vox_idx] = max(src_diff.stat);
end
if length(vox_idx)>1, fprintf('\n***** WARNING: *****\nSelecting only one of multiple peak voxels. User should check!\n'); vox_idx = vox_idx(1); end

%get mni coordinates
mni_coord_cm = src_diff.pos(vox_idx,:);
mni_coord_mm = mni_coord_cm * 10; %convert to mm

%rename (function output)
src_peak = struct;
src_peak.max = vox_max;
src_peak.idx = vox_idx;
src_peak.mni = mni_coord_mm;


%% plot beamformer contrast

if plot_src == 1
    
    %read in template mri
    templatefile = fullfile(fieldtrip_path, '/template/anatomy/single_subj_T1.nii');
    template_mri = ft_read_mri(templatefile);
    
    %distinguish between perc and stat
    if      strcmp(diff_metric,'percentage')
        parameter = 'pow';
    elseif  strcmp(diff_metric,'statistics')
        parameter = {'stat' 'prob' 'mask'};
    end
    
    %interpolate source with template mri
    cfg              = [];
    cfg.parameter    = parameter;
    cfg.interpmethod = 'nearest';
    src_intr  = ft_sourceinterpolate(cfg, src_diff, template_mri);
    clear parameter
    
    %add coordinates and check units
    src_intr.coordsys = 'mni';
    src_intr = ft_convert_units(src_intr, 'mm'); %make sure it's mm
    
    %distinguish between perc and stat
    if      strcmp(diff_metric,'percentage')
        funparameter        = 'pow';
        maskparameter       = [];
    elseif  strcmp(diff_metric,'statistics')
        funparameter        = 'stat';
        maskparameter       = 'mask';
    end
    
    %plot orthogonally
    cfg                     = [];
    cfg.method              = 'ortho';
    cfg.funparameter        = funparameter;
    cfg.maskparameter       = maskparameter;
    cfg.location            = mni_coord_mm; %[x y z], head coordinates
    cfg.locationcoordinates = 'head'; %coordinate system used in cfg.location, 'head' or 'voxel' (head coordinates as mm or cm)
    cfg.funcolormap         = 'jet';
%     cfg.funcolorlim         = [-100 100]; %FIXME
    cfg.atlas               = fullfile(fieldtrip_path, '/template/atlas/aal/ROI_MNI_V4.nii');
    
    ft_sourceplot(cfg, src_intr);
    drawnow;
    
    %quick save
    if plot_save == 1
       saveas(gcf, strrep(src_plotfile, '.mat','.png')) % replace '.mat' with '.png'
        % close
    end
    
    %make space
    clear src_intr
end


%% virtual sensors

%read epoched data again, but no filtering this time...
cfg = [];
cfg.dataset = dataset_path;
cfg.trialdef.eventtype  = markertype;
cfg.trialdef.eventvalue = markervalue;
cfg.trialdef.prestim    = markerpre;
cfg.trialdef.poststim   = markerpos;
cfg_data = ft_definetrial(cfg);

%...not even high-pass this time...
%cfg_data.channel    = {'MEG'};
%cfg_data.channel={'MEG', '-A102' , '-A107', '-A146', '-A40'}; %'-A102',

cfg_data.precision  = 'single';% 'double';% 
cfg_data.detrend     = 'yes';% cfg_data.demean='yes'; %LM - replaced demean with detrend (same as above)
%cfg_data.hpfilter   = 'yes';
%cfg_data.hpfreq     = 3;
data_all = ft_preprocessing(cfg_data);

cfg=[];
data_all=ft_denoise_pca(cfg,data_all);

cfg=[];
cfg.channel={'MEG', '-A102' , '-A107', '-A146', '-A40'}; %'-A102',
data_all = ft_preprocessing(cfg,data_all);


%we no longer downsample to 1kHz - we get Aston to match by interpolation
if data_all.fsample == 1200
    data_res = data_all;
else
    cfg = [];
    cfg.resamplefs  = 1200;
    cfg.detrend     = 'no';
    data_res = ft_resampledata(cfg, data_all);
end
clear data_all

data_raw = data_res; %LM - rename here, just for consistency..

% %exclude BAD trials
% cfg = [];
% cfg.trials = trials_all(~ismember(trials_all,trials_bad));
% 
% if ~isempty(cfg.trials) && config.data.excludebad == 1
%     fprintf('\nBAD trials will be removed (%d out of %d).\n', length(trials_bad), length(trials_all))
%     data_raw = ft_preprocessing(cfg, data_res);
% else
%     data_raw = data_res;
% end

%%%%%%%%%%%%%
% SAVE data %
%%%%%%%%%%%%%

% %save raw data (only demeaned)
% if dat_save == 1
%     fprintf('\nSaving raw data into %s\n', dat_savefile)
%     save(dat_savefile, '-v7.3', '-append', 'data_raw')
% end


%convert beamformer weights
vox_flt = cell2mat(wts(vox_idx));

%simulate ft struct
virt_label = ['[' num2str(mni_coord_cm(1), '%.1f') ' ' num2str(mni_coord_cm(2), '%.1f') ' ' num2str(mni_coord_cm(3), '%.1f') ']'];
data_virt = struct;
data_virt.label = {['mni ' virt_label]};
data_virt.time = data_raw.time;
data_virt.fsample = data_raw.fsample;

[a1,a2]=match_str(src.cfg.channel,data_raw.label);
for iTrial=1:length(data_raw.trial)
  data_virt.trial{iTrial} = vox_flt * data_raw.trial{iTrial}(a2,:);
end

%make space
clear cfg_data data_all data_res vox_flt virt_label


end %function

% 
% cfg=[];
% cfg.method='mtmconvol';
% %cfg.channel={'A72'};
% cfg.output='pow';
% cfg.taper='hanning';
% %cfg.keeptrials='yes';
% cfg.taper='dpss';
% cfg.tapsmofrq=4;
% cfg.foi=[30:100];
% cfg.toi=[-1:0.02:2];
% cfg.t_ftimwin=ones(size(cfg.foi))*0.5;
% freq=ft_freqanalysis(cfg,data_virt);
% 
% cfg=[];
% cfg.baseline=[-1 -0.25];
% cfg.baselinetype='relchange';
% fb=ft_freqbaseline(cfg,freq);
% 
% 


% 
%         