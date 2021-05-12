clear;clc;
close all;
% read and prepare the data 
n=5;% kernel size (L_e by voxel)

% type = 'Loose';
% type = 'Dense';
% type = 'Steel';
type = 'Sub_Steel';
load([type,'_porosity_',num2str(n),'.mat']);
load([type,'_packing.mat']);

sz = size(output);
SE=strel('cube',n); %Erode strel
strel_sz=(n-1)/2;

output=padarray(output,[strel_sz strel_sz strel_sz ]);

% numOfske = 0;%skeleton voxel numher
% numOfvox = 0;%total voxel numher

%% Directional porosity calculation
%Directional porosity along z-axis
for i = 1: sz(1)
    
    x_layers = bwconvhull(squeeze(BW(i,:,:)));
    J = imerode(x_layers,SE);
    convex_ind=find(J==1);
    por_layers=output(i,:,:);
    convex_por_x(i)=sum(por_layers(convex_ind))/length(convex_ind);
end
for i = 1: sz(2)
    
    y_layers = bwconvhull(squeeze(BW(:,i,:)));
    J = imerode(y_layers,SE);
    convex_ind=find(J==1);
    por_layers=output(:,i,:);
    convex_por_y(i)=sum(por_layers(convex_ind))/length(convex_ind);
end

for i = 1: sz(3)
    
    z_layers = bwconvhull(BW(:,:,i));
    J = imerode(z_layers,SE);
    convex_ind=find(J==1);
    por_layers=output(:,:,i);
    convex_por_z(i)=sum(por_layers(convex_ind))/length(convex_ind);
end

%% Envelop porosity calculation
% bounding box 
Env_porosity =1 - sum(BW,'all')/(sz(1)*sz(2)*sz(3));

% Env_porosity_theo =1- ((4/3*pi*(50)^3*64)/(400^3)); % Theoretical env_prosity of loose packing

% convex_hull
stats = regionprops3(BW,'ConvexVolume');


Convol=stats.ConvexVolume;
Convex_porosity = sum(Convol)/400^3;

save([type '_' num2str(n) '_Dimensional_porosity.mat'],'convex_por_x','convex_por_y','convex_por_z','Env_porosity');