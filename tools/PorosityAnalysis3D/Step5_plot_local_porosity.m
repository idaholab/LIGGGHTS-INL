clear;clc;close all
addpath(genpath('BW_figures'));

% type = 'Loose';
% type = 'Dense';
% type = 'Steel';
type = 'Sub_Steel';
n=5; % Kernel size
kernels = num2str(n);  %% select kernel size
filename = [type '_porosity_' kernels,'.mat'];
load(filename);

%% Prepare variables


rows=sum(sum(output,1),2);
cols=sum(sum(output,1),3);
zs=sum(sum(output,2),3);
rows = squeeze(rows);
cols = squeeze(cols);
[z,ZI]=max(rows);
[y,YI]=max(cols);
[x,XI]=max(zs);

xslice = [XI];
yslice = [YI];
zslice = [ZI];


%% Plot the 3D slices
 f = figure


fontsize = 16;

% Select colormap
map = [flipud(jet)];
% map = [parula];
% map = [hsv];
% map = [cool];
% map = [spring];
% map = [summer];
% map = [autumn];
% map = [winter];
% map = [gray];


colormap(map)
hold on
h = slice(output,xslice,[],[]);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;

h = slice(output,[],yslice,[]);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;
%
h = slice(output,[],[],zslice);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;

axis image

box on
bx = gca;
bx.LineWidth = 2;
colorbar
caxis([0 1])

xlabel X(voxel), ylabel Y(voxel), zlabel Z(voxel)

ax = gca;
ax.FontSize = fontsize;
grid on
view([39, 29])
axis image


 
%% Plot 2D slices
%Prepare section position
sz = size(output);
n = 4; %select how many slices you need
xpos = linspace(1,sz(1),n);
xpos = floor(xpos);
ypos = linspace(1,sz(2),n);
ypos = floor(ypos);
zpos = linspace(1,sz(3),n);
zpos = floor(zpos);

%% Local porosity distribution in three directions

%%Xslices
f = figure
s=slice(output,xpos,[],[])
for i = 1:n
    s(i).EdgeColor = 'none'
end
colormap(map);
grid off 
axis image
axis off
c = colorbar('FontSize', fontsize);
c.Label.String = 'Porosity';



%%Yslices

f = figure
s=slice(output,[],ypos,[]);
for i = 1:n
    s(i).EdgeColor = 'none'
end
axis image
axis off
colormap(map);

c = colorbar('FontSize', fontsize);
c.Label.String = 'Porosity';



%%Zslices

f = figure;

s=slice(output,[],[],zpos);
for i = 1:n
    s(i).EdgeColor = 'none'
end
colormap(map);

c = colorbar('FontSize', 16);
c.Label.String = 'Porosity';
axis image
axis off
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

