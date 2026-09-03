function [T_onset,T_peak,T_01] = tgaparam(data)
% data = St01;

char_yield = data(end-1,3)
data = cleandata2(data);

% onset temperature

[peaks,locs] = findpeaks(data(:,3),'MinPeakHeight',0.005);
[dda_onset,idx_onset] = max(data(1:locs(1),5));
T_onset = data(idx_onset,1)

%peak temperature

[da_peak,idx_peak] = max(data(:,3));
T_peak = data(idx_peak,1)

% Tempearature at 10% loss

idx01 = data(:,4)>0.01;
idx = find(idx01,1,'first');
T_01 = data(idx,1)

%temperature of end of cellulose decomposition
[dda_c,idx_c] = min(data(:,5));
T_c = data(idx_c,1)

DTG_max = da_peak


