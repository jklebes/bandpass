function image_out = bandpass(image, low_cutoff, high_cutoff, varargin)
% 2D image bandpass filter with options, 
% defaults replicating imageJ's FFT bandpass: butterworth, mirror padding
% author jklebes 2022

%%%%% parse inputs
p = inputParser;
%image
addRequired(p,'image',@isnumeric); %would like to validate array is 2D
%low_cutoff: less than image size, 0->None
addRequired(p,'low_cutoff', @isnumeric)
if low_cutoff==0
    low_cutoff=[];
end
%high_cutoff: more than low cutoff, less than image size
addRequired(p,'high_cutoff', @isnumeric)
%stripe supression: Horizontal, Vertical, or None (default)
addParameter(p,'stripeOption', "None", @(x) any(validatestring(x,["Horizontal", "Vertical", "None"])))
%stripe supression width: less than image dimension, 0->None
addParameter(p,'stripeWidth', 0, @isnumeric)
%filter type, default Butterworth
addParameter(p,'filter', "butterworth", @(x) any(validatestring(x,["butterworth", "hard"])))
addParameter(p,'butterworthN', 2)
%FT padding type, default "symmetric"
padOptions=["symmetric", "mirror","None"];
addParameter(p,'padOption', "symmetric", @(x) isnumeric(x)||any(validatestring(x,padOptions)));

%run the parser
parse(p,image, low_cutoff, high_cutoff,varargin{:})
%additional argument consequences
if p.Results.stripeWidth<=0
    stripeOption="None";
else 
    stripeOption=p.Results.stripeOption;
end

%save image size
image_size=size(image);
%pad and Fourier transform
fourier = fft_padded(image, p.Results.padOption);
fourier_shifted = fftshift(fourier);
if p.Results.filter=="butterworth"
    %construct Fourier space butterworth bandpass filter
    %same (usually square) size as image
    masksize_x = size(fourier_shifted,2); %switch dimensions in case no padding
    masksize_y = size(fourier_shifted,1);
    center_coord_x = floor(masksize_x/2)+1;
    center_coord_y = floor(masksize_y/2)+1;
    n=p.Results.butterworthN;
    mask=butterworthMask(masksize_x, masksize_y, low_cutoff, high_cutoff,n);
elseif p.Results.filter=="hard"
    %construct hard cutoff Fourier space mask
    mask = zeros(size(fourier_shifted));
    masksize_x = size(mask,2);
    masksize_y = size(mask,1);
    center_coord_x = floor(masksize_x/2)+1;
    center_coord_y = floor(masksize_y/2)+1;
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
if stripeOption=="Horizontal"
    %Fourier space stripe goes other way than in real space!
    for col = center_coord_x-width/2:center_coord_x+width/2
        mask(:,col)=0;
    end
elseif stripeOption=="Vertical"
    for row = center_coord_y-width/2:center_coord_y+width/2
        mask(row,:)=0;
    end
end
%apply mask
fourier_masked = fourier_shifted .* mask;
%center
fourier_masked_shift = fftshift(fourier_masked);
%transform back and cut away padding
padded_image_out = real(ifft2(fourier_masked_shift));
left_border = (masksize_x-image_size(1))/2;
top_border = (masksize_y-image_size(2))/2;
image_out= padded_image_out(left_border+1:left_border+image_size(1), ...
    top_border+1:top_border+image_size(2));
end

function mask= butterworthMask(masksize_x, masksize_y, low_cutoff, high_cutoff,n)
    %array of coordinate values from center
    xs= -(masksize_x/2):masksize_x/2-1;
    xs = repmat(xs, [masksize_y 1]);
    ys= -(masksize_y/2):masksize_y/2-1;
    ys = repmat(ys', [1 masksize_x]);
    %array of distances
    dist=sqrt(xs.^2+ys.^2);
    n2=n*2;
    %trivially scaled arrays dist/low_cutoff, dist/high_cutoff
    %butterworth equations on low, high side
    high_filter=(1./(1 + (dist/high_cutoff).^(n2)));
    low_filter=1./(1 + (dist/low_cutoff).^(n2));
    %combine
    mask = high_filter.*(1-low_filter);
end
