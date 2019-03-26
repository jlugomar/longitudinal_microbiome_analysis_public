function outfile = bdnmbwrite(infile, mb, varname)
% writes a .BDN file that BayesWare Discoverer will open
% consisting of the markov blanket (input MB)
% for the network represented in the .BDN file named in INFILE
% parameter VARNAME defaults to 'phenotype'
% 
% 
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


if (nargin < 3)
    varname = 'phenotype';
end
mb{length(mb)+1} = varname;

fid=fopen(infile);
outtext = '';
line = fgetl(fid);
while ischar(line);
    cmd = textscan(line, '%s %[^\n]');
    if (~strcmpi(cmd{1,1}, 'defnode'))
        outtext = strcat(outtext,line,'\n');
    else
        break;
    end
    line = fgetl(fid);
end

while ischar(line)
    lcells = textscan(line, '%s %s %[^\n]', 'BufSize', 50000);
    nodename = lcells{1,2};
    
    % only the defnodes are important:
    cmd = lcells{1,1};
    if (~strcmpi(cmd, 'defnode'))
        line = fgetl(fid);
        continue;
    end

    % find this node name in our MB list
    if (~findinmb(nodename, mb))
        line = fgetl(fid);
        continue;
    end

    % take parsed elements and reconstruct them:
    linetext = sprintf('%s%s%s', lcells{1,1}{1}, ' ', nodename{1});
    
    rest = lcells{1,3};
    lineremainder = textscan(rest{1}, '(%[^)]) (%[^)]) %[^\n]', 'BufSize', 50000);
    nplist = {};
    parents = lineremainder{1,2};
    if (~isempty(parents))
        % we found parents:
        linetext = sprintf('%s%s%s%s', linetext, ' (', lineremainder{1,1}{1}, ')');
        parentlist = textscan(parents{1}, '%s');
        % scan through parents to make sure each is in the MB
        nump = 0;
        for i = 1:length(parentlist{1})
            p = parentlist{1}{i,1};
            if (findinmb(p,mb))
                nump = nump + 1;
                nplist{nump} = p;
            end
        end
        % now assemble the parent string:
        if (nump > 0)
            pstring = '';
            for i = 1:length(nplist)
                if (i > 1)
                    pstring = sprintf('%s%s%s', pstring, ' ', nplist{i});
                else
                    pstring = sprintf('%s%s', pstring, nplist{i});
                end
            end
            linetext = sprintf('%s%s%s%s', linetext, ' (', pstring, ')');
        else
            % none of the parents listed are in the MB
            linetext = sprintf('%s%s', linetext, ' NIL');
        end
        linetext = sprintf('%s%s%s%s', linetext, ' ', lineremainder{1,3}{1}, '\n');        
    else
        % no parents
        linetext = sprintf('%s%s', line, '\n');
    end
    outtext = sprintf('%s%s', outtext, linetext);
    line = fgetl(fid);
end
fclose(fid);

fnamecells = textscan(infile, '%[^.]%[.]%s');
outfile = strcat(fnamecells{1,1}, '_MB.', fnamecells{1,3});

fid = fopen(outfile{1},'w');
fprintf(fid, outtext);
fclose(fid);


function bool = findinmb(name, mblist)

bool = false;
for i = 1:length(mblist)
    if (strcmpi(mblist{i}, name))
        bool = true;
        return;
    end
end





