function image = read_image_file(path)
l=split(path, '.');
file_extension = l(end);
image_file_types={'png', 'jpg', 'jpeg', 'tif', 'tiff', 'mat'};
assert(any(validatestring(file_extension,image_file_types)));
%if .mat object, accept if array
if strcmp(file_extension,'mat')
end
%try tiff
if strcmp(file_extension,'tif')||strcmp(file_extension, 'tiff')
end
%everything else should be ok with imread
image=imread(path);
end