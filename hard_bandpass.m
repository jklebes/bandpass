%author jk

function image_out = hard_bandpass(image,low_cutoff, high_cutoff)
%UNTITLED2 hard low and high frequency filter
%   Detailed explanation goes here
% problem: only works on perfectly square images
fourier = fft2(image);
fourier_shifted = fftshift(fourier);
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
fourier_masked = fourier_shifted .* mask;
fourier_masked_shift = fftshift(fourier_masked);
image_out = real(ifft2(fourier_masked_shift));
end
