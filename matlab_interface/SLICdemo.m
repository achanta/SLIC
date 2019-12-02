%======================================================================
%SLIC superpixels demo
%======================================================================
% Copyright (C) 2019, Radhakrishna Achanta,
% Ecole Polytechnique Federale de Lausanne, Switzerland.
% 
% Please also read the copyright notice in the file slicmex.c and slic.c
%
% This version of the code can handle grey, color, as well as
% images with any other number of channels.
%======================================================================
%Input parameters are:
%[1] Image with any number of channels
%[2] Number of required superpixels (optional, default is 200)
%[3] Compactness factor (optional, default is 10)
%[4] Color conversion flag, if true, converts RGB colors to CIELAB. This
%    flag only takes effect if the number of image channels is 3.
%
%Ouputs are:
%[1] labels (in raster scan order)
%[2] number of labels in the image (same as the number of returned
%superpixels
%
% NOTE:
% The number of returned superpixels may be slightly different from the
% input number of superpixels.
%======================================================================

%----------------------------------
% compile the mex file before usage
%----------------------------------
mex 'slicmex.c'; % mex file compiled

img = imread('bee.png');
% img = rgb2gray(img);
size(img);

tic
[labels, numlabels] = slicmex(img,500,20.0,true);%numlabels is the same as number of superpixels
toc
fprintf('number of superpixels created: %d\n',numlabels);
imagesc(labels); % view the labels

