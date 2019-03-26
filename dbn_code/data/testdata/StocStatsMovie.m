function StocStatsMovie(s_stats, n, fname)
%% script for ex-stats images of stochastic search (simulated anealing)
%
% show image of 
% s_stats.addededge{} list of edges added in sequence

if (nargin < 3)
    fname = 'searchstats_movie2.mp4';
end

SAVEINT = 150;

diffs = s_stats.lldiffs(2:end) - s_stats.lldiffs(1:end-1);
diffs = int16(diffs);

figh = figure();
set(figh,'DefaultImageCDataMapping','direct');
mainh = subplot_tight(3,1,[1,2],[0.1,0.06]);
title(mainh,'Adjacency Matrix, Change in Log Likelihood');
% get a color map with white in the middle, blue and red at opposite ends:
axismax = max(abs(diffs));
cmap = cbrewer('div', 'PiYG',1 + 2*double(axismax),'PCHIP');
credblue = cbrewer('qual','Set1',9);
cmap = [cmap;credblue(1,:);credblue(2,:)];
colormap(cmap);
cremove = size(cmap,1)-1; % red
cadd = size(cmap,1); % blue


% set it up so zero is in the middle:
%caxis([0,2*axismax]);
%caxis([0,2*axismax]+1);
maxdiff = int16(axismax);

% set up image, a matrix:
img = zeros(n,'int16') + maxdiff;
previmg = zeros(n,'int16') + maxdiff;
%img = zeros(n,'int16') + maxdiff+1;
%previmg = zeros(n,'int16') + maxdiff+1;

%% do an image
image(img,'Parent',mainh);
hCB = colorbar('peer',mainh,'YTick',1+[0,axismax/2,axismax,1.5*axismax,2*axismax],...
    'YTickLabel',{num2str(-1*axismax),num2str(-0.5*axismax),'0',num2str(0.5*axismax),num2str(axismax)});
sh1 = subplot_tight(3,1,3,[0.1,0.06]);
plot(sh1,1,10*s_stats.numedges(1),'b',1,s_stats.lldiffs(1),'r');
title(sh1,'Network Log Likelihood (red), # of Edges (blue)');

% MPEG-4 is maybe 10 times better than AVI:
aviobj = VideoWriter(fname,'MPEG-4');
aviobj.FrameRate = 60;
open(aviobj);
M = getframe(figh);
writeVideo(aviobj, M);
adjmat = zeros(n);
for i = 1:length(s_stats.addededge)
    img = previmg;
    % new edges "flash" when first appearing:
    if (~iscell(s_stats.addededge{i}))
        if (isempty(s_stats.addededge{i}))
            continue;
        end
        if (1 == adjmat(s_stats.addededge{i}(1),s_stats.addededge{i}(2)))
            % actually removing this edge
            img(s_stats.addededge{i}(1),s_stats.addededge{i}(2)) = cremove;
            adjmat(s_stats.addededge{i}(1),s_stats.addededge{i}(2)) = 0;
        else
            img(s_stats.addededge{i}(1),s_stats.addededge{i}(2)) = cadd;
            adjmat(s_stats.addededge{i}(1),s_stats.addededge{i}(2)) = 1;
        end        
    else
        for j = 1:length(s_stats.addededge{i})
            if (1 == adjmat(s_stats.addededge{i}{j}(1),s_stats.addededge{i}{j}(2)))
                % actually removing this edge
                img(s_stats.addededge{i}{j}(1),s_stats.addededge{i}{j}(2)) = cremove;
                adjmat(s_stats.addededge{i}{j}(1),s_stats.addededge{i}{j}(2)) = 0;
            else
                img(s_stats.addededge{i}{j}(1),s_stats.addededge{i}{j}(2)) = cadd;
                adjmat(s_stats.addededge{i}{j}(1),s_stats.addededge{i}{j}(2)) = 1;
            end
        end
    end
    % edges that were added in the last frame go back to a normal color:
    if (i > 1)
        for j = 1:length(s_stats.addededge{i-1})
            if (~iscell(s_stats.addededge{i-1}))
                img(s_stats.addededge{i-1}(1),s_stats.addededge{i-1}(2)) = ...
                    diffs(i-1) + maxdiff;
            else
                img(s_stats.addededge{i-1}{j}(1),s_stats.addededge{i-1}{j}(2)) = ...
                    diffs(i-1) + maxdiff;
            end
        end
    end
    image(img,'Parent',mainh);
    title(mainh,'Adjacency Matrix, Change in Log Likelihood');
    previmg = img;
    hCB = colorbar('peer',mainh,'YTick',1+[0,axismax/2,axismax,1.5*axismax,2*axismax],...
        'YTickLabel',{num2str(-1*axismax),num2str(-0.5*axismax),'0',num2str(0.5*axismax),num2str(axismax)});
    plot(sh1,1:i,10 * s_stats.numedges(1:i),'b',1:i,s_stats.lldiffs(1:i),'r');
    title(sh1,'Network Log Likelihood (red), # of Edges (blue)');

    M(rem(i-1,SAVEINT)+1) = getframe(figh);
    if (mod(i,SAVEINT) == 0)
        writeVideo(aviobj, M);
    elseif (i == length(s_stats.addededge))
        writeVideo(aviobj, M(1:rem(i-1,SAVEINT)+1));
    end
end
close(aviobj);

%%


