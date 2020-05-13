function [ config, beam_method, dip_fitting, fixedori, forw_method, diff_metric ] = cub_ft_beam_checkmethod( config )
%[ config, beam_method, dip_fitting, fixedori, forw_method, diff_method ] = cub_ft_beam_checkmethod( config )
%   LM 2-Feb-2016


%missing config.method
if ~isfield(config,'method') || isempty(config.method),
                                                                            config.method.beam = 'sam'; %default beamformer
                                                                            config.method.dfit = 'spinning'; %default dipole fitting
                                                                            config.method.forw = 'localspheres'; %default forward model
end


%missing config.method.beam
if isfield(config,'method') && ~isempty(config.method),
    
    %NOTE: dip_fitting and forw_method are dealt with below...
    if  ~isfield(config.method,'beam') || isempty(config.method.beam),      config.method.beam = 'sam'; end %default beamformer
    
end
beam_method = config.method.beam;
beam_method_options = {'sam' 'lcmv' 'dics'};

%forbidden beamformer method
if      ~any(cell2mat(strfind(beam_method_options, beam_method))),          error('Please specify a valid beamformer method as config.method.beam'); end

%SAM defaults
if      strcmp(beam_method, 'sam'),
    if  ~isfield(config.method,'dfit') || isempty(config.method.dfit),      config.method.dfit = 'spinning'; end %default SAM dipole fitting
    dip_fitting = config.method.dfit;
    dip_fitting_options = {'spinning' 'robert' 'gareth'};
    
    %forbidden SAM dipole fitting method
    if  ~any(cell2mat(strfind(dip_fitting_options, dip_fitting))),          error('Please specify a valid SAM dipole fitting method as config.method.dfit');
    
    %rename SAM dipole fitting method
    else                                                                    fixedori = dip_fitting;
    end
    
%LCMV defaults
elseif  strcmp(beam_method, 'lcmv'),
    if  ~isfield(config.method,'dfit') || isempty(config.method.dfit),      config.method.dfit = 'fixedori'; end %default LCMV dipole fitting
    dip_fitting = config.method.dfit;
    dip_fitting_options = {'fixedori' 'free'};
    
    %forbidden LCMV dipole fitting method
    if  ~any(cell2mat(strfind(dip_fitting_options, dip_fitting))),          error('Please specify a valid LCMV dipole fitting method as config.method.dfit'); end
    
    %rename LCMV dipole fitting method
    if      strcmp(dip_fitting,'fixedori'),                                 fixedori = 'yes';
    elseif  strcmp(dip_fitting,'free'),                                     fixedori = 'no';
    end
    
%DICS defaults
elseif  strcmp(beam_method, 'dics'),
    if  ~isfield(config.method,'dfit') || isempty(config.method.dfit),      config.method.dfit = 'fixedori'; end %default DICS dipole fitting
    dip_fitting = config.method.dfit;
    dip_fitting_options = {'fixedori' 'free'};
    
    %forbidden DICS dipole fitting method
    if  ~any(cell2mat(strfind(dip_fitting_options, dip_fitting))),          error('Please specify a valid DICS dipole fitting method as config.method.dfit'); end
    
    %rename DICS dipole fitting method
    if      strcmp(dip_fitting,'fixedori'),                                 fixedori = 'yes';
    elseif  strcmp(dip_fitting,'free'),                                     fixedori = 'no';
    end
    
end %beam_method options


%missing config.method.forw
if  ~isfield(config.method,'forw') || isempty(config.method.forw),
    %SAM
    if      strcmp(beam_method, 'sam'),
        config.method.forw = 'localspheres'; %default SAM forward model
    %LCMV
    elseif  strcmp(beam_method, 'lcmv'),
        config.method.forw = 'singleshell'; %default LCMV forward model
    %DICS
    elseif  strcmp(beam_method, 'dics'),
        config.method.forw = 'singleshell'; %default DICS forward model
    end
    fprintf(['\nSetting ' config.method.forw ' as default forward model for the ' beam_method 'beamformer.\n'])
end
forw_method = config.method.forw;

%define label
if strcmp(forw_method,'singleshell'),
    forw_string = 'sinshel';
elseif strcmp(forw_method,'localspheres'),
    forw_string = 'locsphe';
elseif strcmp(forw_method,'singlesphere'),
    forw_string = 'sinsphe';
else
    forw_string = forw_method;
end
config.method.forw_str = forw_string;


%define the metric to use
if  ~isfield(config.method,'diff') || isempty(config.method.diff)
    config.method.diff = 'statistics';
else
    if      strcmp(config.method.diff, 'statistics'),
    elseif  strcmp(config.method.diff, 'percentage'),
    else
        error('Please specify a valid method to compute the source power difference between baseline and active');
    end
end
if      strcmp(config.method.diff, 'statistics'),
    fprintf(['\nComputing difference in power as t ' config.method.diff '.\n'])
elseif  strcmp(config.method.diff, 'percentage'),
    fprintf(['\nComputing difference in power as ' config.method.diff 'change.\n'])
end
diff_metric = config.method.diff;

%string
if      strcmp(config.method.diff, 'statistics'),
    config.method.diff_str = 'stat';
elseif  strcmp(config.method.diff, 'percentage'),
    config.method.diff_str = 'perc';
end

end %function

