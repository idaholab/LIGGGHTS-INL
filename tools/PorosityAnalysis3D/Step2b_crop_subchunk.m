%% This is an optional function to crop a subchunk

clc;
clear;
close all
dbstop  if error
% read the 3D binary matrix

% type = 'Loose';
% type = 'Dense';
type = 'Steel';

load([type,'_packing.mat']);
sz = size(BW);

%% user inputs to select the position and size of the subchun.
% e.g. a 100^3 subchunk
mix_x = 1; 
min_y = 1;
min_z = 1;
wid = 99;
len = 99;
hei = 99;

% crop and svae
BW = imcrop3(BW,[mix_x min_y min_z wid len hei]);
save (['Sub_',type,'_packing.mat'],'BW');
