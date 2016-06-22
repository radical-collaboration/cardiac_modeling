function ei = EI(x)

x=x';
global data_x;
global hyp;
global post;
global D;
global xbest;
[mu sd2] = testGp(hyp,{@meanZero}, {@covSEiso}, {@likGauss}, full(data_x(1:D,:)), post, x);
sig = sqrt(sd2);
%ei = -log_exp_imp(-xbest, -(mu-0.0001), sig);
rho = 0.001;

diff = (mu - xbest - rho);
Z = diff./sigma;
npdf = normpdf(Z);
ncdf=normcdf(Z);
acq = -(ncdf.*diff + npdf.*sig);


%expected_improvement = (diff) .* normcdf(Z) + sigma .* normpdf(Z')';

