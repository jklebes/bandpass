function image_out = butterworth(image, low_cutoff, high_cutoff, n, varargin)
    %padding option mirror ("symmetric") by default to correspond to 
    %ImageJ FFT bandpass
    if nargin>4
        pad_option=varargin(1);
    else
        pad_option = "symmetric";
    end
    %save image size
    image_size=size(image);
    %pad and Fourier transform
    fourier = fft_padded(image, pad_option);
    %shift to center
    fourier_shifted = fftshift(fourier);
    %construct Fourier space butterworth bandpass filter
    mask = zeros(size(fourier_shifted)); %same (square) size as image
    masksize_x = size(mask,2); %switch dimensions in case no padding
    masksize_y = size(mask,1);
    center_coord_x = floor(masksize_x/2)+1;
    center_coord_y = floor(masksize_y/2)+1;
    for col = 1:masksize_x
        for row =1:masksize_y
            dist = sqrt((col-center_coord_y)^2+(row-center_coord_x)^2);
            dist_from_low_cutoff= dist/low_cutoff; %catch possible divide by 0!
            dist_from_high_cutoff = dist/high_cutoff;
            % Butterworth filter - equation
            mask(row,col)= (1/(1 + (dist_from_high_cutoff)^(2*n))).*(1-1/(1 + (dist_from_low_cutoff)^(2*n)));
        end
    end
    %apply filter to Fourier space image
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