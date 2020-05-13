

%----------------------
%  Source anlyses
%----------------------


% addpath '/analyse/Project0179/Scripts/'
% addpath '/analyse/Project0179/Data/MRI/'
% addpath '/share/apps/fieldtrip_new/'
% ft_defaults
% 
% %addpath('/.../full-path-to-folder-containing-the-analysis-functions')
% %megdir = dir('MEG');
% 
% beam_method_arr = {'dics' 'lcmv'};% 'sam'
% 
% cd '/analyse/Project0179/Data/MRI/'
% 
% goodsub=[1001:1012 1014 1016:1023];
% 
% for isubj=1:length(goodsub);
%     isub=goodsub(isubj);
%     cd '/analyse/Project0179/Data/MRI/'
%     
%     % starting here
%     dd=dir([num2str(isub),'.mri']);
%     mri=ft_read_mri([num2str(isub), '.mri/' dd(5).name]);
%     
%     %cd to directory of MEG data
%     cd (['/analyse/Project0179/Data/MEG/', num2str(isub), '/0179a04KP/']);
%     dir('*');
%     filelist=dir('*');
%     cd(filelist(1).name);temp_folder1=dir('*');
%     cd(temp_folder1(3).name);temp_folder2=dir('*');
%     cd(temp_folder2(3).name);temp_folder4=dir('*');
%     hs=ft_read_headshape('hs_file');
%     hs.pnt=hs.pos;
%     
%     %upper is left
%     cfg=[];
%     cfg.method='interactive';
%     cfg.coordsys='4d';
%     mri2=ft_volumerealign(cfg,mri);
%     
%     cfg.method='headshape';
%     cfg.headshape=hs;
%     mri=ft_volumerealign(cfg,mri2);
%     cd '/analyse/Project0179/Subject/'
%     save([num2str(isub),'MRI.mat'],'mri')
% end




clear all

goodsub=[1001:1005 1007 1008:1011 1014:1020 1022:1023];

for isubj=1:length(goodsub);
    isub=goodsub(isubj);
    for icond=1:6;
    
    subject_path=(['/analyse/Project0179/Data/Subject_cond_',num2str(icond), '/S', num2str(isub), '/']);
    %subject_path=(['k:/Data/Subject_cond_',num2str(icond), '/S', num2str(isub), '/']);
    
    % Change here for different conditions
    cd (['/analyse/Project0179/Data/MEG_CMPR/S', num2str(isub), '/', num2str(icond), '/']); % chnage/1 for diff cond
    %cd (['K:/Data/MEG_CMPR/S', num2str(isub), '/', num2str(icond), '/']); % chnage/1 for diff cond

    dir('*');
    path_now = pwd;
    path_now1 = pwd;
    path_now2 = ([path_now1,'/']);
    dataset_path = path_now2;
    
    get_dir = ([path_now, '/c,rfDC']);
    dataset_path2 = (['' get_dir'']);
    
    mrifile_name     = ([ num2str(isub),'MRI.mat']);  
    
    %data
    markertype          = 'TRIGGER'; 
    markervalue         = 4; 
    timewin             = [-2.0 2.0];
    bsln.timewin        = [-1.0 0.0];
    actv.timewin        = [0.0 1.0];
    %sam and lcmv options
    cov.timewin  = [0 0.2];
    cov.bpfreq   = [1 40];
    
    %definitions
    markerpre           = timewin(1)*-1;
    markerpos           = timewin(2);
    bsln_timewin        = bsln.timewin;
    actv_timewin        = actv.timewin;
    
    %covariance and CSD
    cov_timewin         = cov.timewin;
    cov_bpfreq          = cov.bpfreq;
    
    %load mri
    load(fullfile(subject_path, mrifile_name),'mri'); %% problem
    mri = ft_convert_units(mri, 'cm');
    
    %compute segmentation
    cfg = [];
    cfg.output = {'brain'; 'skull'; 'scalp'};
    mri_seg = ft_volumesegment(cfg, mri);
    %add anatomical information to the segmentation
    mri_seg.transform = mri.transform;
    mri_seg.anatomy   = mri.anatomy;
    
%     %check segmentation
%     if plot_seg == 1
%         cfg              = [];
%         cfg.funparameter = 'brain';% 'gray';%
%         ft_sourceplot(cfg,mri_seg);
%         drawnow;
%         
%         %quick save
%         if plot_save == 1
%             saveas(gcf, strrep(seg_plotfile, '.mat', '.png'))
%             % close
%         end
%     end
%     
%     %save segmented mri, if requested
%     if seg_save == 1
%         fprintf('\nSaving segmented mri file as %s\n', seg_savefile)
%         save(seg_savefile, '-struct', 'mri_seg')
%     end
    
    %-----------------
    %% head model
    %----------------
    
    %get grad structure from header file
    hdr = ft_read_header(dataset_path2);
    grad = hdr.grad;
    grad = ft_convert_units(grad, 'cm'); %make sure units match...
    
    %singleshell in native space
    cfg = [];
    cfg.method = 'singleshell';
    hdm = ft_prepare_headmodel(cfg, mri_seg);
    hdm = ft_convert_units(hdm, 'cm');
    
    %make space
    clear mri_seg
    
    %save headmodel
    %fprintf('\nSaving volume conduction model (headmodel) as %s\n', hdm_savefile)
    save(['/analyse/Project0179/Data/Project1/MEG/Struct/Source/HDM_' num2str(isub),'_icond', num2str(icond), '.mat'], 'hdm','-v7.3');
    %save(hdm_savefile, '-struct', 'hdm')
    
    
    %-------------------------------------
    %% source model and leadfield
    %------------------------------------
    
    %% load template source model
    load /analyse/Project0179/Fieldtrip_versions/fieldtrip-20160116/template/sourcemodel/standard_sourcemodel3d8mm.mat
    %load M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502\template/sourcemodel/standard_sourcemodel3d8mm.mat
    
    template_sourcemodel = sourcemodel; %sourcemodel is Fieldtrip's varname
    clear sourcemodel
    
    % Inverse-warp the subject specific grid to the atlas based template grid
    % now compute %convert to MNI template space
    cfg                 = [];
    cfg.grid.warpmni    = 'yes';
    cfg.grid.template   = template_sourcemodel;
    cfg.grid.nonlinear  = 'yes';
    cfg.mri             = mri;
    cfg.inwardshift     = -1.5;
    
    sourcemodel = ft_prepare_sourcemodel(cfg);
    
%        %check headmodel
%     if plot_hdm == 1
%         figure; hold on;
%         ft_plot_vol(hdm);%, 'edgecolor', 'none')
%         %alpha 0.4
%         ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
%         ft_plot_sens(grad);
%         hold off;
%         drawnow;
%         
%         %quick save
%         if plot_save == 1
%             saveas(gcf, strrep(hdm_plotfile, '.mat', '.png'))
%             % close
%         end
%     end
    
    %----------------------------------
    %compute the leadfield matrices
    %----------------------------------
    
    %%%%% ** MIGHT NOT NEED THIS SECTION REALLY ** %%%%%%%
    
    %note: we now define the grid resolution based on the pre-computed
    %sourcemodel, which is the irregular grid warped onto the template...
    cfg                 = [];
    cfg.channel         = {'MEG'}; %REMOVE BAD CHAN??
    cfg.grid            = sourcemodel; %  individual head model from the previous step, the sensor array and the sourcemodel.
    cfg.vol             = hdm;
    cfg.grad            = grad;
    cfg.normalize       = 'yes'; %'yes' would remove depth bias (Q in eq. 27 of van Veen et al, 1997)
    
    leadfield = ft_prepare_leadfield(cfg); % **** this was commented out.. need it??
    
    %save sourcemodel and leadfield variables, along with the template used
    %fprintf('\nSaving sourcemodel and leadfield as %s\n', ldf_savefile)
    save(['/analyse/Project0179/Data/Project1/MEG/Struct/Source/ldf_' num2str(isub),'_icond', num2str(icond), '.mat'], 'sourcemodel', 'leadfield', '-v7.3');

        %save(ldf_savefile, 'sourcemodel', 'leadfield', 'template_filename')
  
    
    %----------------------------
    % read in epoched data
    %---------------------------
    
    cfg = [];
    cfg.dataset = dataset_path2;
    cfg.trialdef.eventtype  = markertype;
    cfg.trialdef.eventvalue = markervalue;
    cfg.trialdef.prestim    = markerpre;
    cfg.trialdef.poststim   = markerpos;
    
    cfg_data = ft_definetrial(cfg);
    
    %----------------------------------------------------------------------
    % preprocess and band-pass filter the data (covariance frequencies)
    %----------------------------------------------------------------------
    
    %time-domain preprocessing
    
    %FIXME: lower frequencies may require padding before band-pass filtering ??
    
    %preprocessing options
    cfg_data.channel    = {'all'};
    cfg_data.precision  = 'single';% 'double';% 
    cfg_data.detrend    = 'yes'; %Glasgow replaced demean with detrend
    cfg_data.bpfilter   = 'yes';
    cfg_data.bpfreq     = cov_bpfreq;
    %cfg.hilbert         = 'complex';
    cfg_data.bpfilttype = 'but';
    cfg_data.bpfiltdir  = 'twopass';
    
    data_prep_all = ft_preprocessing(cfg_data);
    
    %downsample to -.-.- Hz - speed things up!
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
    
    %-----------------------------------------------------------------
    % compute global covariance, i.e. using the whole time-window
    %-----------------------------------------------------------------
    
    %covariance matrix
    cfg = [];
    cfg.channel={'MEG', '-A102' , '-A107', '-A146', '-A40'}; %'-A102',
    cfg.removemean = 'no';
    cfg.covariance = 'yes';
    cfg.covariancewindow = cov_timewin;% 'all';
    
    tlck = ft_timelockanalysis(cfg, data_prep);
    data_beam = tlck;
        
    %cfg=[];
    %cfg.hilbert='complex';
    %databp=ft_preprocessing(cfg,tmp);
    %cfg=[];
    %cfg.toilim=[-0.5 0.5];
    %databp=ft_redefinetrial(cfg,databp);
        
    %--------------------------------------------------------------
    % compute common filters, i.e. using the whole time-window
    %--------------------------------------------------------------
    
    %common spatial filter
    cfg = [];
    cfg.method      = 'lcmv';
    cfg.grid        = sourcemodel; % computed above from temp
    cfg.headmodel   = hdm;
    cfg.grad        = grad;
    cfg.normalize   = 'yes';
    cfg.channel     = data_beam.label;
    leadfield       = ft_prepare_leadfield(cfg);
    leadfield       = rmfield(leadfield, 'cfg'); % Kev added from JG
    cfg.grid        = leadfield;
        
    %beamformer options
    cfg.lcmv.fixedori      = 'yes';
    %cfg.lcmv.keepfilter    = 'yes';
    cfg.lcmv.projectnoise  = 'yes';
    cfg.lcmv.lambda        = '5%'; %regularization parameter
    
    src = ft_sourceanalysis(cfg, data_beam);
    
    %replace the source coordinates with the template ones
    src.pos = template_sourcemodel.pos;
    src.dim = template_sourcemodel.dim;
    src=rmfield(src,'cfg');
    src=rmfield(src,'mom');
 
    %save source
    %if src_save == 1
    %fprintf('\nSaving source common data as %s\n', src_savefile)
    save(['/analyse/Project0179/Data/Project1/MEG/Struct/Source/LCMV_' num2str(isub),'_icond', num2str(icond), '.mat'], 'src','-v7.3');
    %save(src_savefile, '-v7.3', 'src')
    %end
    
    end
end
    
    
    
    
    
    
    
    
    
    %store the common spatial filter (aka beamformer weights)
    wts = src.avg.filter;
    
    %-------------------------------------------------------------
    %% plot, just for fun, as there's no contrast here yet!
    %-------------------------------------------------------------    
    
%     templatefile =  'M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502/template/anatomy/single_subj_T1.nii';
%     template_mri = ft_read_mri(templatefile);
%     template_mri = ft_volumereslice([], template_mri);
%     
%     cfg              = [];
%     cfg.parameter    = 'avg.pow';
%     %cfg.downsample=4;
%     cfg.interpmethod = 'nearest';
%     
%     src_plot = ft_sourceinterpolate(cfg, src, template_mri);
%     src_plot.coordsys = 'mni';
%     src_plot = ft_convert_units(src_plot, 'mm'); %make sure it's mm
%     
%     %plot orthogonally
%     cfg                     = [];
%     cfg.method              = 'slice';
%     cfg.funparameter        = 'pow';
%     cfg.locationcoordinates = 'head';
%     %cfg.funcolormap         = 'jet';
%     cfg.atlas               = 'M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502/template/atlas/aal/ROI_MNI_V4.nii';
%     %cfg.downsample           = 4;
%     
%     ft_sourceplot(cfg,src_plot);
    
    %----------------------------------
    %% baseline and active windows
    %----------------------------------
    
    %we will then beam these using the common filter computed above
    %cfg=[];
    %cfg.hilbert='complex';
    %data_prep=ft_preprocessing(cfg,data_prep1);
    cfg = [];
    cfg.toilim = bsln_timewin;
    data_prep_bsln = ft_redefinetrial(cfg, data_prep);
    
    %cfg=[];
    %cfg.hilbert='complex';
    %data_prep=ft_preprocessing(cfg,data_prep2);    
    cfg = [];
    cfg.toilim = actv_timewin;
    data_prep_actv = ft_redefinetrial(cfg, data_prep);
    
    %------------------------------
    % time-domain: covariance
    %-----------------------------
    
    % KEV - added
    %timelock analysis of bandpass-filtered baseline and active windows
    cfg             = [];
    cfg.channel     = {'MEG','-A102' , '-A107', '-A146', '-A40'}; %'-A102',
    cfg.keeptrials  = 'yes';
    cfg.covariance  = 'yes'; %this is crucial to actually get the power estimates
    cfg.covariancewindow = cov_timewin;% KEV ADDED NEW
    avg_kp = ft_timelockanalysis(cfg,data_prep); % KEV - added
       
    
    %timelock analysis of bandpass-filtered baseline and active windows
    cfg             = [];
    cfg.channel     = {'MEG','-A102' , '-A107', '-A146', '-A40'}; %'-A102',
    cfg.keeptrials  = 'yes';
    cfg.covariance  = 'yes'; %this is crucial to actually get the power estimates
    tlck_bsln = ft_timelockanalysis(cfg, data_prep_bsln);
    tlck_actv = ft_timelockanalysis(cfg, data_prep_actv);
    
    data_beam_bsln = tlck_bsln;
    data_beam_actv = tlck_actv;
    
    %--------------
    % SAVE data 
    %--------------
    
    %save baseline and active data 'to be beamed' (either tlck or freq)
    % if dat_save == 1
    %     fprintf('\nSaving baseline and active ''data_beam'' into %s\n', dat_savefile)
    %     save(dat_savefile, '-v7.3', 'data_beam_bsln', 'data_beam_actv')
    % end
    
    %--------------------
    %% Beamforming                          **** ERROR ****
    %--------------------
    
    %get source-power estimates of bsln and actv separately
    cfg             = [];
    cfg.method      = 'lcmv';
    cfg.grid        = leadfield; % grid we renamed it to leadfield
    cfg.grid.filter = wts;       % ** src.avg.filter; ** taken from whole cov timwin
    cfg.headmodel   = hdm;       % headmodel = vol ?
    %cfg.vol        = hdm;       % headmodel = vol ?
    cfg.grad        = grad;
    
    cfg.channel     = data_beam.label;
    
    cfg.rawtrial    = 'yes';  % avg.pow ?? lost
    %cfg.keeptrials  = 'yes'; % need this as well?
    
    %beamformer options
    cfg.lcmv.fixedori      = 'yes';
    cfg.lcmv.keepfilter    = 'yes';
    cfg.lcmv.projectnoise  = 'yes';
    cfg.lcmv.lambda        = '5%'; %regularization parameter
    
    src_bsln = ft_sourceanalysis(cfg, data_beam_bsln);
    src_actv = ft_sourceanalysis(cfg, data_beam_actv);
    
    %update the pos field...
    src_bsln.pos = template_sourcemodel.pos;
    src_bsln.dim = template_sourcemodel.dim;
    src_actv.pos = template_sourcemodel.pos;
    src_actv.dim = template_sourcemodel.dim;
    
    
%     % Kev -added: FIXME: 
%     cfg = [];
%     cfg.parameter = 'avg.pow';
%     cfg.operation = '((x1-x2)./x2)*100';
%     S1bl=ft_math(cfg,src_actv,src_bsln);
    
    %-------------------------------------------------------------
    % plot, No contrast here yet!
    %-------------------------------------------------------------    
    
    %templatefile =  'M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502/template/anatomy/single_subj_T1.nii';
    templatefile =   '/analyse/Project0179/Fieldtrip_versions/fieldtrip-20160116/template/anatomy/single_subj_T1.nii';
    template_mri = ft_read_mri(templatefile);
    template_mri = ft_volumereslice([], template_mri);
    
    cfg              = [];
    cfg.parameter    = 'avg.pow';
    %cfg.downsample  = 4;
    cfg.interpmethod = 'nearest';
    
    src_plot = ft_sourceinterpolate(cfg, src_actv, template_mri);
    src_plot.coordsys = 'mni';
    src_plot = ft_convert_units(src_plot, 'mm'); %make sure it's mm
    
    %plot orthogonally
    cfg                     = [];
    cfg.method              = 'slice';
    cfg.funparameter        = 'pow';
    cfg.locationcoordinates = 'head';
    %cfg.funcolormap        = 'jet';
    %cfg.atlas               = 'M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502/template/atlas/aal/ROI_MNI_V4.nii';
    %cfg.downsample         = 4;
    
    ft_sourceplot(cfg,src_plot);
    
    %----------------------------------------
    %% compute difference in power as...
    %----------------------------------------
    
    cfg = [];
    cfg.parameter   = 'avg.pow';
    cfg.method      = 'analytic';
    cfg.statistic   = 'depsamplesT';
    cfg.tail        = 0;
    cfg.alpha       = 0.05;% 0.01;%
    
    nTrls = length(src_actv.trial);
    design = zeros(2,2*nTrls);
    for i = 1:nTrls
        design(1,i)         = i;
        design(1,nTrls+i)   = i;
    end
    design(2,1:nTrls)           = 1;
    design(2,nTrls+1:2*nTrls)   = 2;
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
  
    src_diff = ft_sourcestatistics(cfg, src_actv1, src_bsln1); %, template_sourcemodel);
    
    %make mask double to have smooth masking later
    src_diff.mask = double(src_diff.mask);
    
    %update the pos field...
    src_diff.pos = template_sourcemodel.pos;
    src_diff.dim = template_sourcemodel.dim;
    
    
% cfg = [];
% cfg.parameter        = 'avg.pow';
% cfg.dim              = src_actv.dim;
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 500;
% 
% ntrials = numel(src_actv.time);
% design  = zeros(2,2*ntrials);
% design(1,1:ntrials) = 1;
% design(1,ntrials+1:2*ntrials) = 2;
% design(2,1:ntrials) = [1:ntrials];
% design(2,ntrials+1:2*ntrials) = [1:ntrials];
%  
% cfg.design   = design;
% cfg.ivar     = 1;
% cfg.uvar     = 2;
    

    
    %----------------------
    %% peak coordinate
    %----------------------
    
    [vox_max, vox_idx] = max(src_diff.stat);
    if length(vox_idx)>1, fprintf('\n***** WARNING: *****\nSelecting only one of multiple peak voxels. User should check!\n'); vox_idx = vox_idx(1); end
    
    %get mni coordinates
    mni_coord_cm = src_diff.pos(vox_idx,:);
    mni_coord_mm = mni_coord_cm * 10; %convert to mm
    
    %rename (function output)
    src_peak = struct;
    src_peak.max = vox_max;
    src_peak.idx = vox_idx;
    src_peak.mni = mni_coord_mm;
    
    %---------------------------------
    %% plot beamformer contrast
    %----------------------------------
    
    %read in template mri
    %templatefile = 'M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502/template/anatomy/single_subj_T1.nii');
    templatefile =   '/analyse/Project0179/Fieldtrip_versions/fieldtrip-20160116/template/anatomy/single_subj_T1.nii';
    template_mri = ft_read_mri(templatefile);
    
    parameter = {'stat' 'prob' 'mask'};
    
    %interpolate source with template mri
    cfg              = [];
    cfg.parameter    = parameter;
    cfg.interpmethod = 'nearest';
    src_intr         = ft_sourceinterpolate(cfg, src_diff, template_mri);
    clear parameter
    
    %add coordinates and check units
    src_intr.coordsys = 'mni';
    src_intr = ft_convert_units(src_intr, 'mm'); %make sure it's mm
    
    funparameter        = 'stat';
    maskparameter       = 'mask';
    
    
    %plot orthogonally
    cfg                     = [];
    cfg.method              = 'ortho';
    cfg.funparameter        = funparameter;
    cfg.maskparameter       = maskparameter;
    cfg.location            = mni_coord_mm; %[x y z], head coordinates
    cfg.locationcoordinates = 'head'; %coordinate system used in cfg.location, 'head' or 'voxel' (head coordinates as mm or cm)
    cfg.funcolormap         = 'jet';
    %cfg.funcolorlim         = [-100 100]; %FIXME
    cfg.atlas               = 'M:\Toolbox\fieldtrip-20160502\fieldtrip-20160502/template/atlas/aal/ROI_MNI_V4.nii');
    
    ft_sourceplot(cfg, src_intr);
    drawnow;
    
    %quick save
    if plot_save == 1
        saveas(gcf, strrep(src_plotfile, '.mat','.png')) % replace '.mat' with '.png'
        % close
    end
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


        
        
        
        
        

             

        

        

        














 















fg=[];
cfg.bpfilter='yes';
cfg.bpfiltord=3;

%cfg.hilbert='complex';
tmp=ft_preprocessing(cfg,data3);


cfg=[];
cfg.covariance='yes';
cfg.covariancewindow=[-0.5 0.5];
avg=ft_timelockanalysis(cfg,tmp);

cfg=[];
cfg.hilbert='complex';
databp=ft_preprocessing(cfg,tmp);
cfg=[];
cfg.toilim=[-0.5 0.5];
databp=ft_redefinetrial(cfg,databp);

cd ../Source
%outname=['MNE_filter' num2str(is) '-cond' num2str(3)];
%load(outname)


load(['Vol_' num2str(isub)])
load(['Grid_' num2str(isub)])

cfg=[];
cfg.method= 'lcmv';
cfg.lcmv.lambda='5%';
cfg.lcmv.keepfilter='yes';
cfg.lcmv.reducerank='yes';
cfg.lcmv.normalize='yes';
cfg.grid= indgrid;
cfg.vol= vol;
cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
cfg.grad=avg.grad;
cfg.reducerank=2;
cfg.normalize='yes';
if (f1==1),
    grid=ft_prepare_leadfield(cfg,avg);
end;
cfg.grid=grid;
source = ft_sourceanalysis(cfg, avg);
filt=cat(1,source.avg.filter{source.inside})*1e12;


































































