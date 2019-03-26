function keep = IDMatch(idlist,data)
% filter down the data matrix by data(:,1) matching an element of idlist
%
%

if (iscell(idlist))
    DOCELL = true;
else
    DOCELL = false;
end

if (DOCELL)
    keep = false(size(data));
    for i = 1:length(data)
        val = data{i};
        if (sum(strcmp(val,idlist)) > 0)
            keep(i) = true;
        end
    end
else
    keep = false(size(data(:,1)));
    for i = 1:length(data(:,1))
        val = data(i,1);
        if (sum(val == idlist) > 0)
            keep(i) = true;
        end
    end
end
