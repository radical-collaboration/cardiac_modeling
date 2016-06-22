function [E] = ucb(x)
global hyp;
global meanfunc;
global covfunc;
global likfunc;
global data_x;
global post;
global D;
global si;
x=x';
[m sd2] = testGp(hyp,meanfunc, covfunc, likfunc, full(data_x(1:D,:)), post, x);
%
% if (length(x)<5)
%     beta =si;
% else 
   beta=si* sqrt(2*log(D^2*pi^2 /(6*0.1)));
% else
 
E = -(m + beta*sqrt(sd2));
end