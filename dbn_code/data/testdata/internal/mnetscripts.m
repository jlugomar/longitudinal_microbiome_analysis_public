% mixed network scripts of nodes -

cowell = repmat(struct('filename',{},'self',{},'parents',{},'children',{},...
    'values',{},'xy',{},'discrete',{},'index',{},'cindex',{},'pindex',{}),...
    1,6);

cfig4 = repmat(struct('filename',{},'self',{},'parents',{},'children',{},...
    'values',{},'xy',{},'discrete',{},'index',{},'cindex',{},'pindex',{}),...
    1,6);

% cfig1 = repmat(struct('filename',{},'self',{},'parents',{},'children',{},...
%     'values',{},'xy',{},'discrete',{},'index',{},'cindex',{},'pindex',{}),...
%     1,8);



cowell(1).self = 'A';
cowell(1).children{1} = 'X';
cowell(1).discrete = true;
cowell(1).values = {1, 0};

cowell(2).self = 'B';
cowell(2).children{1} = 'X';
cowell(2).discrete = true;
cowell(2).values = {1,2,3};

cowell(3).self = 'X';
cowell(3).parents{1} = 'A';
cowell(3).parents{2} = 'B';
cowell(3).children{1} = 'Y';
cowell(3).discrete = false;

cowell(4).self = 'C';
cowell(4).children{1} = 'Y';
cowell(4).children{2} = 'Z';
cowell(4).discrete = true;
cowell(4).values = {1,0};

cowell(5).self = 'Y';
cowell(5).parents{1} = 'X';
cowell(5).parents{2} = 'C';
cowell(5).children{1} = 'Z';
cowell(5).discrete = false;

cowell(6).self = 'Z';
cowell(6).parents{1} = 'Y';
cowell(6).parents{2} = 'C';
cowell(6).discrete = false;


cfig4(1).self = 'A';
cfig4(1).children{1} = 'B';
cfig4(1).children{2} = 'C';
cfig4(1).discrete = true;
cfig4(1).values = {0, 1};

cfig4(2).self = 'B';
cfig4(2).parents{1} = 'A';
cfig4(2).parents{2} = 'C';
cfig4(2).parents{3} = 'E';
cfig4(2).parents{4} = 'D';
cfig4(2).discrete = false;

cfig4(3).self = 'C';
cfig4(3).children{1} = 'B';
cfig4(3).children{2} = 'E';
cfig4(3).parents{1} = 'A';
cfig4(3).parents{2} = 'D';
cfig4(3).discrete = false;

cfig4(4).self = 'D';
cfig4(4).children{1} = 'C';
cfig4(4).children{2} = 'B';
cfig4(4).parents{1} = 'F';
cfig4(4).discrete = false;

cfig4(5).self = 'E';
cfig4(5).children{1} = 'B';
cfig4(5).parents{1} = 'C';
cfig4(5).parents{2} = 'F';
cfig4(5).discrete = false;

cfig4(6).self = 'F';
cfig4(6).children{1} = 'D';
cfig4(6).children{2} = 'E';
cfig4(6).discrete = false;


% cfig1(1).self = 'A';
% cfig1(1).parents{1} = 'B';
% cfig1(1).discrete = false;
% 
% cfig1(2).self = 'B';
% cfig1(2).parents{1} = 'E';
% cfig1(2).children{1} = 'A';
% cfig1(2).children{2} = 'D';
% cfig1(2).children{3} = 'F';
% cfig1(2).discrete = false;
% 
% cfig1(3).self = 'C';
% cfig1(3).parents{1} = 'E';
% cfig1(3).parents{2} = 'G';
% cfig1(3).children{1} = 'F';
% cfig1(3).discrete = false;
% 
% cfig1(4).self = 'D';
% cfig1(4).parents{1} = 'B';
% cfig1(4).parents{2} = 'E';
% cfig1(4).children{1} = 'F';
% cfig1(4).discrete = false;
% 
% cfig1(5).self = 'E';
% cfig1(5).children{1} = 'B';
% cfig1(5).children{2} = 'D';
% cfig1(5).children{3} = 'C';
% cfig1(5).discrete = false;
% 
% cfig1(6).self = 'F';
% cfig1(6).parents{1} = 'B';
% cfig1(6).parents{2} = 'C';
% cfig1(6).parents{3} = 'D';
% cfig1(6).children{1} = 'H';
% cfig1(6).discrete = false;
% 
% cfig1(7).self = 'G';
% cfig1(7).children{1} = 'C';
% cfig1(7).discrete = false;
% 
% cfig1(8).self = 'H';
% cfig1(8).parents{1} = 'F';
% cfig1(8).discrete = false;

