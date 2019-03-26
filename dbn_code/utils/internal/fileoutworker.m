function fileoutworker(fname, titlearray, data, middletitles)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

titlearray = char(titlearray);
[cs1,cs2] = size(titlearray);
c = repmat('\t',cs1,1);
titlearray = [titlearray c];

if (nargin > 3)
    middletitles = char(middletitles);
    [ts1,ts2] = size(middletitles);
    tabs = repmat('\t',ts1,1);
    middletitles = [middletitles tabs];
end



fid = fopen(fname,'wt');
if (nargin < 4)
    fprintf(fid, titlearray');
else
    fprintf(fid, titlearray(1,:)');
    fprintf(fid, middletitles');
    fprintf(fid, titlearray(2:end,:)');
end
fprintf(fid, '\n');
fclose(fid);
dlmwrite(fname,data,'-append','delimiter', '\t', 'newline','pc','precision', 2);