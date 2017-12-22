function [ rmse, out ] = ss_robust_metric2( sm_reader, gt_reader, cp_img, varargin )
% Inputs:
% sm_reader - a video object with saliency maps that should be adjusted
% gt_reader - a video object with ground truth saliency maps
% cp_img - an image with precomputed center prior image
%
% Optional inputs:
% skip - sampling parameter, only every skip-th frame will be used in cost function
% nbins - set size of map vector (for example, you could use nbins=65536 if you use 16-bit images)
% skip_saccade - if is true then skips frames with saccades
%
% Outputs:
% rmse - root of mean square error between GT and adjusted saliency maps
% out - a structure describing computed transform that fits source saliency maps to ground truth saliency maps with the following fields:
%   beta - a coefficient to blend source saliency map with center prior image
%   map - 256-element vector describing contrast transform of source saliency maps, so that GT = map[SM(i, j)] + beta CP(i, j)
%   psnr - PSNR measure between GT and adjusted saliency maps


parser = inputParser;
addParamValue(parser, 'skip',       1,      @isnumeric);
addParamValue(parser, 'normalize',  false,  @islogical);
addParamValue(parser, 'nbins',      256,    @isnumeric);
addParamValue(parser, 'fair_mse',   false,  @islogical);
addParamValue(parser, 'writer',     []                );
addParamValue(parser, 'skip_saccade',false, @islogical);
parse(parser, varargin{:});
in = parser.Results;

if exist('cp_img', 'var')
    if ischar(cp_img)
        cp_img = imread(cp_img);
    end
    cp_img = xrgb2gray(im2double(cp_img));
else
    %TODO: estimate CP automatically
end

assert(in.nbins >= 2);
assert(sm_reader.Height == gt_reader.Height && sm_reader.Width == gt_reader.Width);
assert(sm_reader.Height == size(cp_img, 1) && sm_reader.Width == size(cp_img, 2));

num_frames = min(sm_reader.NumberOfFrames, gt_reader.NumberOfFrames);

dg = zeros(in.nbins, 1); %diagonal of H expect last elem
bc = zeros(in.nbins, 1); %column for center prior blending coef
sum_gt_per_bin  = zeros(in.nbins, 1);
sum_cp_sqr      = 0;
sum_gt_weighted = 0;
sum_gt_sqr      = 0;
frames_processed = 0;
log_frames = 1 + [1:2:9, 10:10:100]/100*(num_frames-1);

for i = 1:in.skip:num_frames
   frame_gt = xrgb2gray( im2double(gt_reader.read(i)) );
   frame_sm = xrgb2gray( im2double(sm_reader.read(i)) );
   
   if in.skip_saccade && mean(frame_sm(:)) < eps
       continue
   end
   
   cur_weighted_gt = sum(sum(frame_gt .* cp_img));
   sum_gt_weighted = sum_gt_weighted + cur_weighted_gt;
   
   sum_gt_sqr = sum_gt_sqr + sum(frame_gt(:) .^ 2);
   
   if in.normalize
       frame_sm = normalize(frame_sm);
   end
   frame_sm_bins = discretize(frame_sm, in.nbins);
   
   dg = dg + hist(frame_sm_bins(:), in.nbins)';
   bc = bc + get_sum_under_slices(cp_img, frame_sm_bins, in.nbins);
   sum_gt_per_bin = sum_gt_per_bin + get_sum_under_slices(frame_gt, frame_sm_bins, in.nbins);
   
   frames_processed = frames_processed + 1;
   
   if ~isempty( find(i-1 < log_frames & log_frames <= i, 1) )
    	fprintf('ss_robust_metric2: processed %d/%d frames\n', round((i-1)/in.skip) + 1, round(num_frames/in.skip));
   end
end

sum_cp_sqr = frames_processed * sum(cp_img(:) .^ 2);
H = spdiags([sum_cp_sqr; dg], 0, in.nbins+1, in.nbins+1);
H(2:end, 1) = bc;
H(1, 2:end) = bc;

f = -[sum_gt_weighted; sum_gt_per_bin];

A_le = spdiags(repmat([1 -1], [in.nbins, 1]), [1 2], in.nbins-1, in.nbins+1);
b = zeros(size(A_le, 1), 1);

lb = zeros(in.nbins+1, 1);
ub = ones(in.nbins+1, 1);

optimizer = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'Display', 'off');
solution = quadprog(H, f, A_le, b, [], [], lb, ub, [], optimizer);
%solution = H \ -f;

out.beta    = solution(1);
out.map     = solution(2:end);
out.H       = H;
out.f       = f;
out.c0      = sum_gt_sqr;

num_pixels = frames_processed * numel(frame_gt);
out.mse_predicted   = (solution' * H * solution + 2 * f' * solution + sum_gt_sqr) / num_pixels;
out.num_pixels      = num_pixels;
out.num_frames      = frames_processed;
out.do_norm         = in.normalize;
out.psnr            = -10 * log(out.mse_predicted) / log(10);

fprintf('ss_robust_metric2: total psnr = %f\n', out.psnr);

if in.fair_mse || ~isempty(in.writer)
    out.mse_real = 0;
    
    for i = 1:in.skip:num_frames
    	frame_gt = xrgb2gray( im2double(gt_reader.read(i)) );
        frame_sm = xrgb2gray( im2double(sm_reader.read(i)) );
        if in.skip_saccade && mean(frame_sm(:)) < eps
            continue
        end
        if in.normalize
            frame_sm = normalize(frame_sm);
        end
        frame_sm_bins = discretize(frame_sm, in.nbins);

        frame_sm_tuned = reshape( out.map(frame_sm_bins(:) + 1), size(frame_sm) );
        frame_sm_tuned = frame_sm_tuned + out.beta * cp_img;

        out.mse_real = out.mse_real + sum( (frame_gt(:) - frame_sm_tuned(:)) .^ 2);

        if ~isempty(in.writer)
           in.writer.writeVideo(frame_sm_tuned, i);
        end
    end
    
    out.mse_real = out.mse_real / num_pixels;
end

if in.fair_mse    
	rmse = sqrt(out.mse_real);
else
    rmse = sqrt(out.mse_predicted);
end

end

function [ out ] = xrgb2gray(in)
    if size(in, 3) ~= 1
        out = rgb2gray(in);
    else
        out = in;
    end
end

function [ out_bins ] = discretize(im, nbins)
	out_bins = round((nbins - 1) * im2double(im));
    %out_bins = int32(out_bins);
	assert(min(out_bins(:)) >= 0 && max(out_bins(:)) < nbins);
end

function [ sm_out ] = normalize(sm_in)
    tol = 0.01;
    sm_out = imadjust(sm_in, stretchlim(sm_in, tol)); 
end

function [ sum_slice ] = get_sum_under_slices( img, bin_map, nbins )

s = sparse(bin_map(:)+1, 1:numel(img), img(:), nbins, numel(img));
sum_slice = full(sum(s , 2));

end