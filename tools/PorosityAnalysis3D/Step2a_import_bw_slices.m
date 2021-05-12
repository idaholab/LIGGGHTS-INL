clc;
clear;
close all
dbstop  if error

[fileName, pathName, ~] =...
    uigetfile(...
    {'*.jpg;*.tif;*.tiff;*.png;*.gif;*.bmp;','All Image Files'},...
    'Select stack images to be processed',...
    'MultiSelect','on'...
    );

%{
Determine the number of files to be processed; >= 0 and <= 2^16-1
    %}
    nFiles = uint16( size(fileName, 2) );
    
    %{
Convert file name from cell to string.
    %}
    fileName = char(fileName);
    
    %{
Get the pixel resolution of each image
    %}
    img = imread( [pathName fileName(1,:)] );
    imgSize = size(img);
    %{
Initiate a 3D matrix for the stack images
    %}
    if ( isa(img, 'double') )
        BW3D = zeros( imgSize(1), imgSize(2), nFiles );
    elseif ( isa(img, 'single') )
        BW3D = zeros( imgSize(1), imgSize(2), nFiles, 'single' );
    elseif ( isa(img, 'uint16') )
        BW3D = zeros( imgSize(1), imgSize(2), nFiles, 'uint16' );
    elseif ( isa(img, 'uint8') )
        BW3D = zeros( imgSize(1), imgSize(2), nFiles, 'uint8' );
    else
        BW3D = false( imgSize(1), imgSize(2), nFiles );
    end
    
    %{
Serial-reading of images to assign values in BW3D matrix
    %}
    for iFile=1:nFiles
        BW3D(:,:,iFile) = imread( [pathName fileName(iFile,:)] );
    end
    %
    
    
    BW = BW3D;
    %{ 
     Optional: binarize the image to logical format and crop the margin (if the image is not logical type)
    %}
    BW = imbinarize(BW3D);
    stats = regionprops3(BW,'Image','BoundingBox');
    BW=stats.Image{1,1};

%{
    save the 3D binary matrix
%}
    save ('Steel_packing.mat','BW');
%     save ('Loose_packing.mat','BW');
%     save ('Dense_packing.mat','BW');
    %{
  Visualize your packing by z-slices
    
   Thanks to Maysam Shahedi for the visualization function.
   (https://www.mathworks.com/matlabcentral/fileexchange/41334-imshow3d)
    %}
imshow3Dfull(BW)