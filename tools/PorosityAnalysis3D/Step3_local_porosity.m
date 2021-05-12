%% Local porosity calculation 
clear;clc;
close all;
addpath(genpath('BW_figures'));

% read the 3D binary matrix
% type = 'Loose';
% type = 'Dense';
% type = 'Steel';
type = 'Sub_Steel';
load([type,'_packing.mat']);

% select kernel size (L_e by voxel)
n=5; % must be an odd numer!

% prepare scanning window matrix and output matxix
scan_window = zeros(n,n,n);
sz=size(BW);
window_center = (n+1)/2;
output = zeros(sz(1)-n,sz(2)-n,sz(3)-n);

%% traverse the entire 3D binary matrix


for i = window_center:sz(1)-window_center+1
    for j = window_center:sz(2)-window_center+1
        for k = window_center:sz(3)-window_center+1
            scan_window = BW(i-(window_center-1):i+(window_center-1),j-(window_center-1):j+(window_center-1),k-(window_center-1):k+(window_center-1));
            output(i-(window_center-1),j-(window_center-1),k-(window_center-1))= 1- (sum(scan_window,'all')/n^3);
            
        end
    end
end

% save the results
save([type,'_porosity_',num2str(n),'.mat'],'output','i','-v7.3');


