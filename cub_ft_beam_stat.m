function [ src_diff ] = cub_ft_beam_stat( cfg, src_actv, src_bsln, template_sourcemodel )
%[ src_diff ] = cub_ft_beam_stat( cfg, src_actv, src_bsln, template_sourcemodel )
%   Detailed explanation goes here


if isempty(cfg)
    
    cfg = [];
    cfg.parameter   = 'pow';
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

end

src_diff = ft_sourcestatistics(cfg, src_actv, src_bsln);

%make mask double to have smooth masking later
src_diff.mask = double(src_diff.mask);

%update the pos field...
src_diff.pos = template_sourcemodel.pos;
src_diff.dim = template_sourcemodel.dim;

end

