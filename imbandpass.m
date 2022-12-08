function image_out = imbandpass(image, low_cutoff, high_cutoff, varargin)
% 2D image bandpass filter with options and stripe supression,
% defaults replicating imageJ's FFT bandpass: gaussian filter, mirror padding
% author jklebes 2022

%%%%% parse inputs
p = inputParser;
%image
addRequired(p,'image',@(x) isnumeric(x)&&ismatrix(x) );
%parse image first to validate further inputs againsts its dimensions
parse(p,image);
dims=size(image);
maxdim = max(dims);
%low_cutoff: less than image size, 0->None
addRequired(p,'low_cutoff', @(x) isempty(x)||(isnumeric(x)&&0<=x&&x<=maxdim));
%high_cutoff: more than low cutoff, less than image size
addRequired(p,'high_cutoff', @(x) isempty(x)||(isnumeric(x) &&(isempty(low_cutoff)||low_cutoff<x)&&x<=maxdim));
%stripe supression: Horizontal, Vertical, or None (default)
addParameter(p,'stripes', 'None', @(x) any(validatestring(x,{'Horizontal', 'Vertical', 'None'})));
%stripe Filter option
addParameter(p,'stripeFilter', 'gaussian', @(x) any(validatestring(x,{'gaussian', 'hard'})));
%stripe supression width: less than image dimension, 0->None
addParameter(p,'stripeWidth', 3, @(x) isnumeric(x))
%filter type, default gaussian
addParameter(p,'filter', 'gaussian', @(x) any(validatestring(x,{'gaussian','butterworth', 'hard'})))
addParameter(p,'butterworthN', 1)
%FT padding type, default 'symmetric'
% accept options of matlab's padarray and likely synonyms (case insensitive)
padOptions={'symmetric', 'mirror','None', 'circular','replicate','zeros', 'periodic'};
addParameter(p,'padOption', 'symmetric', @(x) isnumeric(x)||any(validatestring(x,padOptions)));

%run the parser
parse(p,image, low_cutoff, high_cutoff,varargin{:})
%additional argument consequences
if isempty(low_cutoff)
    low_cutoff=0;
    %in the Gaussian filter, this automatically
    %adds no filter, as e^0=1.  In the other options the value is
    %symbolic
end
if isempty(high_cutoff)
    %would like to set to Inf, but a very small central disk
    %of the Fourier space mask needs to be removed for 
    %it to function properly
    high_cutoff=maxdim*100;
    %and e^(-infty)=0 hopefuly
end
if p.Results.stripeWidth<=0
    stripes='None';
else
    stripes=p.Results.stripes;
end

%save image size
image_size=size(image);
%pad and Fourier transform
fourier = fft_padded(image, p.Results.padOption);
fourier_shifted = fftshift(fourier);
%size of mask to construct
%same (usually square) size as image
masksize_x = size(fourier_shifted,1);
masksize_y = size(fourier_shifted,2);
center_coord_x = floor(masksize_x/2)+1;
center_coord_y = floor(masksize_y/2)+1;
switch p.Results.filter
    case 'gaussian'
        %construct Gaussian filter
        %calculate factor in exponent of Gaussian from filterLarge / filterSmall
        high_cutoff_ratio = high_cutoff/maxdim;
        low_cutoff_ratio = low_cutoff/maxdim;
        mask=gaussianMask(masksize_x, masksize_y, low_cutoff_ratio, high_cutoff_ratio);
    case 'butterworth'
        %construct Fourier space butterworth bandpass filter
        n=p.Results.butterworthN;
        mask=butterworthMask(masksize_x, masksize_y, low_cutoff, high_cutoff,n);
    case 'hard'
        %construct hard cutoff Fourier space mask
        mask = zeros([masksize_x, masksize_y]);
        for col = 1: masksize_y
            for row = 1:masksize_y
                distance = sqrt((col-center_coord_y)^2+(row-center_coord_x)^2);
                if distance > low_cutoff && (distance <= high_cutoff ||high_cutoff==0)
                    mask(row, col)=1;
                end
            end
        end
end
%potentially add stripe to Fourier space mask
switch p.Results.stripeFilter
    case "hard"
        %add a stripe of 0s to the mask
        switch stripes
            case 'Horizontal'
                %Fourier space stripe goes other way than in real space!
                width=p.Results.stripeWidth;
                for col = floor(center_coord_x-width/2):ceil(center_coord_x+width/2)
                    mask(:,col)=0;
                end
            case 'Vertical'
                width=p.Results.stripeWidth;
                for row = floor(center_coord_y-width/2):ceil(center_coord_y+width/2)
                    mask(row,:)=0;
                end
        end
    case "gaussian"
        %add a gaussian stripe filter to the mask
        switch stripes
            case 'Horizontal'
                %Fourier space stripe goes other way than in real space!
                xs= -(masksize_x/2):masksize_x/2-1;
                xs = repmat(xs, [masksize_y 1]);
                stripe_cutoff_ratio = p.Results.stripeWidth/(2*image_size(1));
                stripeMask=exp(-stripe_cutoff_ratio^2*xs.^2);
                mask=mask.*stripeMask;
            case 'Vertical'
                ys= -(masksize_y/2):masksize_y/2-1;
                ys = repmat(ys', [1 masksize_x]);
                stripe_cutoff_ratio = p.Results.stripeWidth/(2*image_size(2));
                stripeMask=exp(-stripe_cutoff_ratio^2*ys.^2);
                mask=mask.*stripeMask;
                %otherwise no stripe filtering action
        end

end
%apply mask
fourier_masked = fourier_shifted .* mask;

%center
fourier_masked = fftshift(fourier_masked);
%transform back and cut away padding
padded_image_out = real(ifft2(fourier_masked));
left_border = ceil((masksize_x-image_size(1))/2);
top_border = ceil((masksize_y-image_size(2))/2);
image_out= padded_image_out(left_border+1:left_border+image_size(1), ...
    top_border+1:top_border+image_size(2));
% results may have range shifted into the negatives, shift back into
% 0 to 255 range
image_out = image_out-min(image_out, [], 'all');
end

function mask = gaussianMask(masksize_x, masksize_y, low_cutoff_ratio, high_cutoff_ratio)
%array of coordinate values from center
%imageJ version centers on (N/2+1, N/2+1),
%but we want to match the centering of fftshift
%on (N/2, N/2)
%arrays of x, y coordinates
xs= -(masksize_x/2):masksize_x/2-1;
xs = repmat(xs, [masksize_y 1]);
ys= -(masksize_y/2):masksize_y/2-1;
ys = repmat(ys', [1 masksize_x]);
%array of distances
dist=sqrt(xs.^2+ys.^2);
%exponentials
mask_low = exp(-low_cutoff_ratio^2*dist.^2);
mask_high = exp(-high_cutoff_ratio^2*dist.^2);
mask = mask_low.*(1-mask_high);
end

function mask= butterworthMask(masksize_x, masksize_y, low_cutoff, high_cutoff,n)
%array of coordinate values from center
xs= -(masksize_x/2):masksize_x/2-1;
xs = repmat(xs, [masksize_y 1]);
ys= -(masksize_y/2):masksize_y/2-1;
ys = repmat(ys', [1 masksize_x]);
%array of distances
dist=sqrt(xs.^2+ys.^2);
n2=n*2; %exponent
%butterworth equations on low, high lengthscale side
high_filter=(1./(1 + (dist/high_cutoff).^(n2)));
%(no high cutoff: hopefully 1/Inf -> 0)
if low_cutoff>0 %in case of no low cutoff, avoid divide by 0
    low_filter=1./(1 + (dist/low_cutoff).^(n2));
else
    low_filter = zeros(masksize_x, masksize_y);
end
%combine
mask = high_filter.*(1-low_filter);
end

function image_out = fft_padded(image, varargin)
%image_out adds padding and carries out fft
% the possibly rectangular image is extended to a square of size N*N,
% where N=2^i is >= 1.5 the larger dimension (except in case "none")
% "symetric" is mirror padding, can be used to replicate imageJ FFt bandpass filter
% https://imagej.nih.gov/ij/plugins/fft-filter.html
if nargin >1
    padding_opt = varargin{1};
    if ~isnumeric(padding_opt) && padding_opt=="mirror"
        padding_opt = "symmetric";
    end
    if ~isnumeric(padding_opt) && padding_opt=="periodic"
        padding_opt = "replicate";
    end
    if ~isnumeric(padding_opt) && padding_opt=="zeros"
        padding_opt=0;
    end
    if isnumeric(padding_opt) | ismember(padding_opt,["circular","replicate","symmetric"])
        dims = size(image);
        maxdim = max(dims);
        i = ceil(log2(maxdim*1.5)); %TODO check
        N=2^i;
        x_pad=N-dims(1);
        y_pad=N-dims(2);
        image=padarray(image,[ceil(x_pad/2), ceil(y_pad/2)], padding_opt);
        if mod(x_pad,2)==1
            image=image(1:end-1,:);
        end
        if mod(y_pad,2)==1
            image=image(:,1:end-1);
        end
    end
end
%else if argument is nonexistent, empty, or invalid option:
%image stays unpadded
image_out = uint8(fft2(image));
end