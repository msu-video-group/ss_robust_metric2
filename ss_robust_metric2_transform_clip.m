function [ sm_transformed ] = ss_robust_metric2_transform_clip( sm, opt, cp_img )
% Applies transformation computed by psnr_robust_metric2 metric
%
% Inputs:
% sm - a video object with source saliency maps
% opt - transformations computed by psnr_robust_metric2 function
% cp_img - center prior image used in psnr_robust_metric2 to compute opt structure
%
% Outputs: a video object with adjusted saliency maps

if exist('cp_img', 'var')
    opt.cp_img = cp_img;
end
opt.cp_img = xrgb2gray(im2double( opt.cp_img ));

if isfield(opt, 'do_norm') && opt.do_norm
    warning('ss_robust_metric2_transform_clip. Using "do_norm" field. Results maybe incosistent!');
end

sm_transformed = FrameTransformer( sm, @(fr) tans(fr, opt) );

end

function [ sm_out ] = tans( sm, opt )

nbins = numel(opt.map);

sm = xrgb2gray(im2double( sm ));
if isfield(opt, 'do_norm') && opt.do_norm
    sm = normalize(sm);
end

sm_bins = discretize(sm, nbins);

sm_out = reshape( opt.map(sm_bins(:) + 1), size(sm) );
sm_out = sm_out + opt.beta * opt.cp_img;
end

function [ sm_out ] = normalize(sm_in)
    tol = 0.01;
    sm_out = imadjust(sm_in, stretchlim(sm_in, tol)); 
end

function [ out_bins ] = discretize(im, nbins)
	out_bins = round((nbins - 1) * im2double(im));
    out_bins = int32(out_bins);
	assert(min(out_bins(:)) >= 0 && max(out_bins(:)) < nbins);
end

function [ out ] = xrgb2gray(in)
    if size(in, 3) ~= 1
        out = rgb2gray(in);
    else
        out = in;
    end
end
