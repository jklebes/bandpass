%author jk

function image_out = stripe_filter(image, width, direction)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fourier = fft2(image);
fourier_shifted = fftshift(fourier);
mask = zeros(size(fourier_shifted));
masksize_x = size(mask,2);
masksize_y = size(mask,1);
center_coord_x = floor(masksize_x/2)+1;
center_coord_y = floor(masksize_y/2)+1;
for col = 1: masksize_y
    for row = 1:masksize_x
        if abs((center_coord_x-col) /2) > width
            mask(row, col)=1;
        end
    end
end
fourier_masked = fourier_shifted .* mask;
fourier_masked_shift = fftshift(fourier_masked);
image_out = real(ifft2(fourier_masked_shift));
end