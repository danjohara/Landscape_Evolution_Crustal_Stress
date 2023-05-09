function createModelMovie(modelFile,filePath,fileName,FR,endTime)
% Name: createModelMovie
% Author: Daniel O'Hara
% Date: 5/14/18

% Description: Creates a movie of the inputted model.

% Input:
%   modelFile: File (with folder path) of the landscape evolution model.
%   filePath: Folder path of the movie.
%   fileName: File name of the movie.
%   FR: Movie frame rate (20-40 seems to be good).
%   endTime: Model time that will be made a movie.

% Load model file information.
load(modelFile)
zz = [fv_all.zt];
zzz = zz;
zzz = zzz-zzz(1,:);
tt = [fv_all.ts];
maxZ = max(zz(:));
maxZZ = max(zzz(:));
minZZ = min(zzz(:));

% Create video writer and set up figure window.
vi = VideoWriter([filePath,'\',fileName,'.mp4']);
vi.FrameRate = FR;
open(vi);
ff = figure;
set(ff,'Position',get(0,'Screensize'))
rr = find(tt > endTime,1);
if isempty(rr)
    rr = length(fv_all.ts);
end

% Create movie.
tt = tt(1:rr);
for i = 1:rr
    plot(fv_all.v{1}.x./1000,zz(i,:)./1000,'k')
    hold on
    c = i;
    plot(fv_all.v{1}.x(fv_all.basCs{c}.singleBCNodesI{1})./1000,zz(i,fv_all.basCs{c}.singleBCNodesI{1})./1000,'LineWidth',2)
    for j = 2:length(fv_all.basCs{c}.singleBCNodesI)
       xC = fv_all.v{1}.x(fv_all.basCs{c}.singleBCNodesI{j})./1000;
       zC = zz(i,fv_all.basCs{c}.singleBCNodesI{j})./1000;
       plot(xC,zC,'LineWidth',2)
    end

    ti1 = find(~isnan(fv_all.zLps{i}),1);
    ti2 = find(~isnan(fv_all.zLps{i}),1,'last');
    ti3 = find(~isnan(fv_all.zLns{i}),1,'last');
    ti4 = find(~isnan(fv_all.zLns{i}),1);
    lX = [fv_all.v{1}.x(ti1:ti2),fliplr(fv_all.v{1}.x(ti4:ti3)),fv_all.v{1}.x(ti1)]./1000;
    lZ = [fv_all.zLps{i}(ti1:ti2),fliplr(fv_all.zLns{i}(ti4:ti3)),fv_all.zLps{i}(ti1)]./1000;
    plot(lX,lZ,'-k','linewidth',2)
   
    hold off
    xlabel('X (km)')
    ylabel('Z (km)')
    tmpTitle = sprintf('Topography - T = %.3e yr',tt(i));
    tmpTitle = strrep(tmpTitle,'e+',' x 10^');
    tmpTitle = [tmpTitle(1:end-5),'{',tmpTitle(end-3),'}',tmpTitle(end-2:end)];
    title(tmpTitle)
    ylim([0 maxZ./1000])

    frame = getframe(ff);
    writeVideo(vi,frame)
end

% Close writer.
close(vi)