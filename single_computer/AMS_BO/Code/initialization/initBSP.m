function obs= initBSP(truep)
% given a parameter setting simulate ECG data with noise
% trup: true parameter setting
load('modelFiles.mat','simu');
[timesteps,state] = inverseModel(simu,truep); % the simulation model
[data,sigPow,noiPow_SD2] = genIniNoi(state,30,'','',0); % add 30db Gaussian noise

% downsample the simulated data in time because the real ECG data is
% sparser than the simulated data
tsparse = 0:(170/450*0.5):170;
timestps = length(tsparse);
index = zeros(1,timestps);
for i  = 1: timestps
    [dummy,index(i)] = min(abs(tsparse(i)-timesteps));
end
obs.bsp  = data(:,index);
obs.time = timesteps(index);
obs.noiseSigma = (repmat(noiPow_SD2,length(tsparse),1));
end