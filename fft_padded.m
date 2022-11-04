function image_out = fft_padded(image, varargin)
%image_out adds padding and carries out fft
% the possibly rectangular image is extended to a square of size N*N,
% where N is >= 1.5 the larger dimension (except in case "none")
% "symetric" is mirror padding, can be used to replicate imageJ FFt bandpass filter
% https://imagej.nih.gov/ij/plugins/fft-filter.html
if nargin >1
    padding_opt = varargin{1};
    if ~isnumeric(padding_opt) && padding_opt=="mirror"
        padding_opt = "symmetric";
    end
    if ~isnumeric(padding_opt) && padding_opt=="zeros"
        padding_opt=0;
    end
    if isnumeric(padding_opt) | ismember(padding_opt,["circular","replicate","symmetric"])
        dims = size(image);
        maxdim = max(dims);
        i = ceil(log2(maxdim*1.5)); %TODO check
        N=i^2;
        image=padarray(image,[N/2, N/2], padding_opt);
    end
end
%else if argument is nonexistent, empty, or invalid option:
%image stays unpadded
image_out = fft2(image);
end