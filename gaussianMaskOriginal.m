function mask = gaussianMaskOriginal(masksize_x, masksize_y, low_cutoff_ratio, high_cutoff_ratio)
%%direct transcirption of imageJ version filter construction (leaving out stripes)
% https://imagej.nih.gov/ij/plugins/download/FFT_Filter.java by Joachim Walter
%used for testing equivalence-
%the loops can be made more efficient by array operations
mask = zeros([masksize_x,masksize_y]);

%loop over rows
counter=1;
for j=1:masksize_x/2-1 %position away from center - start counting at 1
    row = j * masksize_x; %index for caulculating memory location also
    backrow = (masksize_x-j)*masksize_x; %keep as if 0-indexed
    rowFactLarge = exp(-(j*j) * scaleLarge);
    rowFactSmall = exp(-(j*j) * scaleSmall);


    % loop over columns
    for col=1:masksize_x/2-1
        backcol = masksize_x+1-col; %but matlab array index is 2
        colFactLarge = exp(- (col*col) * scaleLarge);
        colFactSmall = exp(- (col*col) * scaleSmall);
        factor = (1 - rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
        %factor = counter;
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

mask(masksize_x/2+1) = (1 - rowFactLarge) * rowFactSmall; % (masksize_x/2,0)
mask(rowmid+1) = (1 - rowFactLarge) * rowFactSmall; % (0,masksize_x/2)
mask(masksize_x/2 + rowmid+1) = (1 - rowFactLarge*rowFactLarge) * rowFactSmall*rowFactSmall; % (masksize_x/2,masksize_x/2)

%loop along row 0 and masksize_x/2
rowFactLarge = exp(- (masksize_x/2)*(masksize_x/2) * scaleLarge);
rowFactSmall = exp(- (masksize_x/2)*(masksize_x/2) * scaleSmall);
for col=1:masksize_x/2-1
    backcol = masksize_x-(col);
    colFactLarge = exp(- (col*col) * scaleLarge);
    colFactSmall = exp(- (col*col) * scaleSmall);
    mask(col+1)= (1 - colFactLarge) * colFactSmall;
    mask(backcol+1)= (1 - colFactLarge) * colFactSmall;
    mask(col+rowmid+1)= (1 - colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
    mask(backcol+rowmid+1) = (1 - colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
end

% loop along column 0 and masksize_x/2
colFactLarge = exp(- (masksize_x/2)*(masksize_x/2) * scaleLarge);
colFactSmall = exp(- (masksize_x/2)*(masksize_x/2) * scaleSmall);
for j=1:masksize_x/2-1
    row = j * masksize_x;
    backrow = (masksize_x-j)*masksize_x;
    rowFactLarge = exp(- (j*j) * scaleLarge);
    rowFactSmall = exp(- (j*j) * scaleSmall);
    mask(row+1) = (1 - rowFactLarge) * rowFactSmall;
    mask(backrow+1) = (1 - rowFactLarge) * rowFactSmall;
    mask(row+masksize_x/2+1)= (1 - rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
    mask(backrow+masksize_x/2+1) = (1 - rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
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