% The algorithm runs in three steps:
% (1) PSF (E) estimation and 
% (2) Image (R) non-blind deconvolution
% (3) mobilenetv2 network
function [U H report] = ProposedAlgorithm(inputImage,size_PSF,ground_truth)

% load parameters
parameters;

if ~exist('size_PSF','var') || isempty(size_PSF)
    size_PSF = [3 3]; 
end

% ground-truth PSF (for comparison)
if ~exist('ground_truth','var')
    ground_truth  = [];
end

% normalize images
[inputImage, norm_m, norm_v] = ImagesNormalize(inputImage);

% PSF ESTIMATION
disp('Estimating PSFs...');

% prepare for multiscale process
L = MSlevels; % number of multiscale levels
if (L<1)
    L = 1;
end
sr = [ maxROIsize(1)./(2.^(L-1:-1:0).'), maxROIsize(2)./(2.^(L-1:-1:0).')]; % ROI sizes at each scale level

% precalculate ROI (on which PSF is calculated) for each scale level
ROI = cell(1,L);
hstarP = cell(1,L);
ROI{L} = RGBChannel(inputImage,sr(L,:));
if ~isempty(ground_truth)
    hstarP{L} = ground_truth;
end
for i = L-1:-1:1
    ROI{i} = imresize(ROI{i+1},0.5);
    if ~isempty(ground_truth)
        hstarP{i} = imresize(hstarP{i+1},0.5);
    end
end

% initial PSF size and set them to delta functions
size_PSF = ceil(size_PSF/2^(L-1));
cen = floor((size_PSF+1)/2);
hi = zeros(size_PSF);
hi(cen(1),cen(2)) = 1; % init PSF as delta impulse

% main estimation loop
report.ms = cell(1,L);
for i = 1:L
	if(PAR.verbose > 0) disp(['size_PSF: ',num2str(size(hi))]); end
    
	hi = hi/sum(hi(:));
	[h u report.ms{i}] = ERSteps(ROI{i},hi,PAR,hstarP{i}); 
    hi = imresize(h,2,'lanczos3'); % upsample for new scale
end


H = h;
H(H<0) = 0;
H = h/sum(H(:));

disp('PSF estimation done.'); 

% NON-BLIND DECONVOLUTION
U = LagrangianApproach(inputImage,H,PAR);

disp('Nonblind deconvolution done.'); 
U = U*norm_v + norm_m;

report.par = PAR;
end

function R = RGBChannel(G,win)
		
	isize = size(G);
	gsize = isize(1:2);
	if size(G,3) > 1 % RGB image
		cind = 2; %green channel
	else
		cind = 1;
	end
	
	if any(gsize < win)
		win = gsize;
	end
	margin = floor((gsize-win)/2);
	R = G(margin(1)+1:margin(1)+win(1),margin(2)+1:margin(2)+win(2),cind);
end

function [I m v] = ImagesNormalize(G)
lb = min(G(:));
ub = max(G(:));

v = ub-lb;
m = lb;
I = (double(G)-m)/v;
end