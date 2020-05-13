
clear

addpath '/analyse/Project0179/Scripts/'
addpath '/analyse/Project0179/Data/MRI/'
%addpath('/.../full-path-to-folder-containing-the-analysis-functions')

%megdir = dir('MEG');

beam_method_arr = {'dics' 'lcmv'};% 'sam'

cd '/analyse/Project0179/Data/MRI/'

goodsub=[1001 1014];

for isubj=1:length(goodsub);
    isub=goodsub(isubj);
    cd '/analyse/Project0179/Data/MRI/'
    
    % starting here
    dd=dir([num2str(isub),'.mri']);
    mri=ft_read_mri([num2str(isub), '.mri/' dd(5).name]);
    
    %cd to directory of MEG data
    cd (['/analyse/Project0179/Data/MEG/', num2str(isub), '/0179a04KP/']);
    dir('*');
    filelist=dir('*');
    cd(filelist(1).name);temp_folder1=dir('*');
    cd(temp_folder1(3).name);temp_folder2=dir('*');
    cd(temp_folder2(3).name);temp_folder4=dir('*');
    hs=ft_read_headshape('hs_file');
    hs.pnt=hs.pos;
    
    %upper is left
    cfg=[];
    cfg.method='interactive';
    cfg.coordsys='4d';
    mri2=ft_volumerealign(cfg,mri);
    
    cfg.method='headshape';
    cfg.headshape=hs;
    mri=ft_volumerealign(cfg,mri2);
    cd (['/analyse/Project0179/Subject/S', num2str(isub)]);
    save([num2str(isub),'MRI.mat'],'mri')
end

clear all

beam_method_arr = {'dics' 'lcmv'};% 'sam'

% All subjects working now
goodsub=[1002]; %% [1001:1023];

%loop over beamformers
for iBeam = 1:1
    
    %make sure the config structured is cleared!
    config = struct;
    
    %mandatory fields
    config.fieldtrip_path   = '/share/apps/fieldtrip_new/';  %'/.../full-path-to-fieldtrip-20160116';
    
    for isubj=1:length(goodsub);
        isub=goodsub(isubj);
        
        config.subject_path=(['/analyse/Project0179/Subject/S', num2str(isub), '/']);    %'/.../full-path-to-folder-X001';
        cd (['/analyse/Project0179/Data/MEG/', num2str(isub), '/0179a04KP/']); 
        dir('*');
        filelist=dir('*');
        cd(filelist(1).name);temp_folder1=dir('*');
        cd(temp_folder1(3).name);temp_folder2=dir('*');
        cd(temp_folder2(3).name);temp_folder4=dir('*');    
        path_now = pwd;
        get_dir = ([path_now, '/']); %c,rfDC']);
        addpath (['' get_dir''])
        config.dataset_path = (['' get_dir'']);
        
        get_dir = ([path_now, '/c,rfDC']);
        addpath (['' get_dir''])
        config.dataset_path2 = (['' get_dir''])
        
        config.mrifile_name     = ([ num2str(isub),'MRI.mat']);    %'X001.mri'; %name of coregistered mri file - needs to be inside the subject's folder
        
        %config.subject_path     = '/analyse/Project0179/Subject/1001/';    %'/.../full-path-to-folder-X001';
        %config.dataset_path     = '/analyse/Project0179/Data/MEG/1001/0179a04KP/16-01-25@1405/1/';    %'/.../full-path-to-dataset';
        %config.mrifile_name     = '1001MRI.mat ';    %'X001.mri'; %name of coregistered mri file - needs to be inside the subject's folder
        
        
        
        %data
        config.data.markertype  = 'TRIGGER'; %'UPPT002';
        config.data.markervalue = 10; %whatever value corresponds to the Grating Onset
        config.data.excludebad  = 1;
        config.data.timewin     = [-2.0 2.0];
        config.bsln.timewin     = [-1.2 0.0];
        config.actv.timewin     = [ 0.3 1.5];
        
        %sam and lcmv options
        config.cov.timewin      = [-1.5 1.5];
        config.cov.bpfreq       = [35 75];
        
        %dics options
        config.csd.timewin      = [-1.5 1.5];
        config.csd.freqarr      = mean(config.cov.bpfreq); %i.e. 55 Hz
        config.csd.freqsmo      = max(config.cov.bpfreq)-mean(config.cov.bpfreq); %i.e. 20 Hz
        
        %save weights
        config.wts.save         = 1;
        
        %save source
        config.src.save         = 0; %make sure src is not saved, as stats take up too much disk space
        
        %leadfield normalisation
        config.ldf.norm         = 'yes'; %although it doesn't really matter, i guess...
        
        %plot options
        config.plot.src         = 1; %plot the source figure
        config.plot.save        = 1; %make sure it's saved
        
        %forward model
        config.method.forw = 'singleshell';
        
        %define beamformer
        config.method.beam = beam_method_arr{iBeam};
        
        %change directory
        %cd (['/analyse/Project0179/Subject'])
        cd(config.subject_path) % Orig script
        
        %run the analysis
        [src_diff, src_peak, data_virt, config] = cub_ft_beam_baselinecond(config);
        
        %save the output
        str_append = ['_' config.method.beam '_' config.method.dfit '_' config.method.forw_str '_' config.method.diff_str];
        file_savename = fullfile(config.subject_path, strrep(config.mrifile_name, '.mri', [str_append '-out.mat']));
        save(file_savename, '-v7.3', 'src_diff', 'src_peak', 'data_virt', 'config')
        
        close
    end
    
end
