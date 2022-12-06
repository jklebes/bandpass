function image_out = bandpass(image, low_cutoff, high_cutoff, varargin)
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
addRequired(p,'low_cutoff', @(x) isnumeric(x)&&0<=x&&x<=maxdim);
%high_cutoff: more than low cutoff, less than image size
addRequired(p,'high_cutoff', @(x) isnumeric(x) &&low_cutoff<x&&x<=maxdim);
%stripe supression: Horizontal, Vertical, or None (default)
addParameter(p,'stripeOption', 'None', @(x) any(validatestring(x,{'Horizontal', 'Vertical', 'None'})));
%stripe supression width: less than image dimension, 0->None
addParameter(p,'stripeWidth', 0, @(x) isnumeric(x))
%stripefilter, current hard, goal Gaussian
addParameter(p,'Stripefilter', 'hard', @(x) any(validatestring(x,{'gaussian','butterworth', 'hard'})))
%filter type, default Butterworth
addParameter(p,'filter', 'gaussian', @(x) any(validatestring(x,{'gaussian','butterworth', 'hard'})))
addParameter(p,'butterworthN', 1)
%FT padding type, default 'symmetric'
%accept options of matlab's padarray and likely synonyms (case insensitive)
padOptions={'symmetric', 'mirror','None', 'circular','replicate','zeros', 'periodic'}; 
addParameter(p,'padOption', 'symmetric', @(x) isnumeric(x)||any(validatestring(x,padOptions)));

%run the parser
parse(p,image, low_cutoff, high_cutoff,varargin{:})
%additional argument consequences
if low_cutoff==0
    low_cutoff=[];
end
if p.Results.stripeWidth<=0
    stripeOption='None';
else 
    stripeOption=p.Results.stripeOption;
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
if p.Results.filter=='gaussian'
    %construct Gaussian filter
    mask=gaussianMaskOriginal(masksize_x, masksize_y, low_cutoff, high_cutoff);
elseif p.Results.filter=='butterworth'
    %construct Fourier space butterworth bandpass filter
    n=p.Results.butterworthN;
    mask=butterworthMask(masksize_x, masksize_y, low_cutoff, high_cutoff,n);
elseif p.Results.filter=='hard'
    %construct hard cutoff Fourier space mask
    mask = zeros([masksize_x, masksize_y]);
    for col = 1: masksize_y
        for row = 1:masksize_y
            distance = sqrt((col-center_coord_y)^2+(row-center_coord_x)^2);
            if distance <= high_cutoff && distance > low_cutoff
                mask(row, col)=1;
            end
        end
    end
end
%potentially add stripe to Fourier space mask
switch stripeOption
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
%apply mask
fourier_masked = fourier_shifted .* mask;
%center
fourier_masked_shift = fftshift(fourier_masked);
%transform back and cut away padding
padded_image_out = real(ifft2(fourier_masked_shift));
left_border = ceil((masksize_x-image_size(1))/2);
top_border = ceil((masksize_y-image_size(2))/2);
image_out= padded_image_out(left_border+1:left_border+image_size(1), ...
    top_border+1:top_border+image_size(2));
end

function mask = gaussianMaskOriginal(masksize_x, masksize_y, low_cutoff, high_cutoff)
%%direct transcirption of imageJ version filter construction (leaving out stripes)
% https://imagej.nih.gov/ij/plugins/download/FFT_Filter.java by Joachim Walter
%used for testing equivalence-
%the loops can be made more efficient by array operations
mask = zeros([masksize_x,masksize_y]);

%calculate factor in exponent of Gaussian from filterLarge / filterSmall
imsize=512; %TODO input original imsize/ input masksize as fraction of
scaleLarge = (high_cutoff/imsize)^2;
scaleSmall = (low_cutoff/imsize)^2;

%loop over rows
counter=1;
for j=1:masksize_x/2-1 %position away from center - start counting at 1
    row = j * masksize_x; %index for caulculating memory location also
    backrow = (masksize_x-j)*masksize_x; %keep as if 0-indexed
    rowFactLarge = exp(-(j*j) * scaleLarge);
    rowFactSmall = exp(-(j*j) * scaleSmall);


    % loop over columns
    for col=1:masksize_x/2-1
        backcol = masksize_x-(col+1); %but matlab array index is 2
        colFactLarge = exp(- (col*col) * scaleLarge);
        colFactSmall = exp(- (col*col) * scaleSmall);
        factor = (1 - rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
        factor = counter;
        mask(col+1+row) = factor; %even though array is 2D, in matlab assigning
                                %to 1D index works - but counts columnwise
        mask(col+1+backrow) = factor;
        mask(backcol+row) = factor;
        mask(backcol+backrow) = factor;
        counter=counter+1;
    end
end

%process meeting points (masksize_x/2,0) , (0,masksize_x/2), and (masksize_x/2,masksize_x/2)
rowmid = masksize_x * (masksize_x/2);
rowFactLarge = exp(- (masksize_x/2)*(masksize_x/2) * scaleLarge);
rowFactSmall = exp(- (masksize_x/2)*(masksize_x/2) * scaleSmall);

mask(masksize_x/2) = 500+masksize_x/2%(1 - rowFactLarge) * rowFactSmall; % (masksize_x/2,0)
mask(rowmid) = 500+rowmid%(1 - rowFactLarge) * rowFactSmall; % (0,masksize_x/2)
mask(masksize_x/2 + rowmid) = 500+masksize_x/2 + rowmid%(1 - rowFactLarge*rowFactLarge) * rowFactSmall*rowFactSmall; % (masksize_x/2,masksize_x/2)

%loop along row 0 and masksize_x/2
rowFactLarge = exp(- (masksize_x/2)*(masksize_x/2) * scaleLarge);
rowFactSmall = exp(- (masksize_x/2)*(masksize_x/2) * scaleSmall);
for col=1:masksize_x/2-1
    backcol = masksize_x-(col+1);
    colFactLarge = exp(- (col*col) * scaleLarge);
    colFactSmall = exp(- (col*col) * scaleSmall);
    mask(col+1)= 200%(1 - colFactLarge) * colFactSmall;
    mask(backcol)= 200%(1 - colFactLarge) * colFactSmall;
    mask(col+1+rowmid)= 200%(1 - colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
    mask(backcol+rowmid) = 200%(1 - colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
end

% loop along column 0 and masksize_x/2
colFactLarge = exp(- (masksize_x/2)*(masksize_x/2) * scaleLarge);
colFactSmall = exp(- (masksize_x/2)*(masksize_x/2) * scaleSmall);
for j=1:masksize_x/2-1
    row = j * masksize_x;
    backrow = (masksize_x-j)*masksize_x;
    rowFactLarge = exp(- (j*j) * scaleLarge);
    rowFactSmall = exp(- (j*j) * scaleSmall);
    mask(row) = 300%(1 - rowFactLarge) * rowFactSmall;
    mask(backrow) = 300%(1 - rowFactLarge) * rowFactSmall;
    mask(row+masksize_x/2)= 300%(1 - rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
    mask(backrow+masksize_x/2) = 300%(1 - rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
end
mask = swapQuadrants(mask);
imshow(mask,[]);
end

function array_out=swapQuadrants(array)
dims = size(array);
xsize=dims(1);
ysize=dims(2);
upperleft=array(1:xsize/2, 1:ysize/2);
upperright=array(1:xsize/2, ysize/2+1:end);
lowerleft=array(xsize/2+1:end,1:ysize/2);
lowerright=array(xsize/2+1:end, ysize/2+1:end);
array_out=[lowerright, lowerleft; upperright, upperleft];
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
    if ~isempty(high_cutoff)
        high_filter=(1./(1 + (dist/high_cutoff).^(n2)));
    else 
        high_filter= ones(masksize_x, masksize_y);
    end
    if ~isempty(low_cutoff)
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
image_out = fft2(image);
end