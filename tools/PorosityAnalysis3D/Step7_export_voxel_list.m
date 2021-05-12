clear;clc;
close all;
% read the 3D binary matrix
type = 'Sub_Steel';
% type = 'Loose';

load([type,'_packing.mat']);
sz = size(BW);


BW=ones(19,19,19);
BW(1:6,1:19,1:6)=0;
    numVoxAll = numel(BW);
    numVoxObj = sum(sum(sum(BW)));
    CC = bwconncomp(BW,6);
    [NumVoxConnObj,Idx] = sort(cellfun(@numel,CC.PixelIdxList),'descend');
    VOX = zeros(numVoxObj, 5);
    
    [R,C,Z] = ind2sub( size(BW), CC.PixelIdxList{Idx(1)} );
numStart = 1;
numEnd = size(R,1);
resolution = 2; %% e.g. 2mm/voxel
for ii = numStart : numEnd
    iii = ii - numStart + 1;
    VOX(ii,1) = R(iii)*resolution;
    VOX(ii,2) = C(iii)*resolution;
    VOX(ii,3) = Z(iii)*resolution;
    VOX(ii,4) = resolution/2; % radius of atom
    VOX(ii,5) = 2; % single atom type in LIGGGHTS

end


tStart = cputime;
outFileName = "vox_list.lmp";
msg = strcat("Export file ", outFileName);
fprintf('%-60s', msg);
fileID = fopen(outFileName,'w');
for i = 1 : numEnd
    fprintf(fileID,'%g %g %g %g %g\n',VOX(i,1),VOX(i,2),VOX(i,3),VOX(i,4),VOX(i,5));
end
fclose(fileID);
tElapsed = cputime - tStart;
fprintf('  %.3f\n', tElapsed);

