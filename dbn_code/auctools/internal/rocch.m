function [tp, fp] = rocch(t, y)
%
% ROCCH - generate a receiver operating characteristic convex hull 
%
%    [TP,FP] = ROCCH(T,Y) gives the true-positive rate (TP) and false positive
%    rate (FP), corresponding to the convex hull of the receiver operating
%    characteristic, where Y is a column vector giving the score assigned to
%    each pattern and T indicates the true class (a value above zero
%    represents the positive class and anything else represents the negative
%    class).  To plot the ROC convex hull,
%
%       PLOT(FP,TP);
%       XLABEL('FALSE POSITIVE RATE');
%       YLABEL('TRUE POSITIVE RATE');
%       TITLE('RECEIVER OPERATING CHARACTERISTIC CONVEX HULL (ROCCH)');
%
%    See [1,2] for further information.
%
%    [1] Fawcett, T., "ROC graphs : Notes and practical
%        considerations for researchers", Technical report, HP
%        Laboratories, MS 1143, 1501 Page Mill Road, Palo Alto
%        CA 94304, USA, April 2004.
%
%    [2] Provost, F. and Fawcett, T., "Robust classification for
%        imprecise environments", Machine Learning, vol. 42,
%        no. 3, pp. 203-231, 2001.
%
%    See also : ROCCH, AUROC

%
% File        : rocch.m
%
% Date        : Wednesday 10th November 2004 
%
% Author      : Dr Gavin C. Cawley
%
% Description : Generate the convex hull of an ROC curve for a two-class
%               classifier (see [1] and [2] for details).
%
% References  : [1] Fawcett, T., "ROC graphs : Notes and practical
%                   considerations for researchers", Technical report, HP
%                   Laboratories, MS 1143, 1501 Page Mill Road, Palo Alto
%                   CA 94304, USA, April 2004.
%
%               [2] Provost, F. and Fawcett, T., "Robust classification for
%                   imprecise environments", Machine Learning, vol. 42,
%                   no. 3, pp. 203-231, 2001.
%
% History     : 10/11/2004 - v1.00
%
% To do       : Add an option to specify the maximum number of points to
%               consider to minimise memory usage.
%
% Copyright   : (c) G. C. Cawley, November 2004.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%

% generate the ROC curve

[tp,fp] = roc(t,y);

tp = [tp ; -0.1];
fp = [fp ; 1.1];

% we are really interested in the convex hull

idx = unique(convhull(fp, tp),'legacy');
fp  = fp(idx(1:end-1));
tp  = tp(idx(1:end-1));

% bye bye...

