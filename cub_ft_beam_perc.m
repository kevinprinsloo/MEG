function [ src_diff ] = cub_ft_beam_perc( cfg, src_actv, src_bsln, template_sourcemodel )
%[ src_diff ] = cub_ft_beam_perc( cfg, src_actv, src_bsln, template_sourcemodel )
%   Detailed explanation goes here


if isempty(cfg)
    
    %source percentage change
    cfg = [];
    cfg.parameter = 'avg.pow';
    cfg.operation = '((x1-x2)./x2)*100';% '(x1-x2)./(x1+x2)';%
    
end

src_diff = ft_math(cfg, src_actv, src_bsln);

%update the pos field...
src_diff.pos = template_sourcemodel.pos;
src_diff.dim = template_sourcemodel.dim;

end

