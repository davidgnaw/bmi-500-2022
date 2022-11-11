% reset the workspace
clear
close all

% load spiral drawing data
d = read_trc("lue-spiral.trc");

% set plotting parameters
TL = [0 5];
nr = 2;
nc = 3;

% plot the left hand marker in x-y-z
marker_name = "L.Finger3.M3";
marker_xyz = d{:,find(names(d) == "L.Finger3.M3") + [0:2]};

t = d{:,"Time"};
t_inds = t>min(TL)&t<max(TL);
t_secs = rem(t(t_inds),1)==0;


% plot
figure
subplot(nr,nc,1)
hold on

% Filter out large, slow movements with a high-pass butterworth filter at 2
% Hz cutoff and filter out jitter with a low-pass butterworth filter at 20
% Hz cutoff. A 6th order filter is fine.


% sampling freq fs is the reciprocal of the difference between two points
fs = 1/mean(diff(t));

% cutoff frequencies for the filter
fc_hi = 2;
fc_lo = 20;
n_order = 6;

% high-pass filter
[b, a] = butter(n_order, fc_hi/(fs/2), "high");
filt1 = filtfilt(b, a, marker_xyz);

% low-pass filter
[y, x] = butter(n_order, fc_lo/(fs/2), "low");
filt2 = filtfilt(y, x, filt1);

% [b,a] = butter(n,Wn) returns the transfer function coefficients of an 
% nth-order lowpass digital Butterworth filter with normalized 
% cutoff frequency Wn [https://www.mathworks.com/help/signal/ref/butter.html]

% calculate the first PC

% find covariance matrix
C = cov(filt2);

% find eigen vector to find PC
[evector, evalue] = eig(C);

% note from above that largest eigenvalue is 50.1646 where eigenvector is
% 0.9732 and at evector(:,3) third column
e = diag(evalue);
[max_evalue, max_index] = max(e);

u = evector(:,max_index);
% filt2 (3600 x 3) * u (3 x 1) = (3600 x 1)
PC = filt2 * u; % new principle component

% calculate projection onto first PC

% transpose u to 3 x 3600 to multiple with PC 3600 x 1 to get 3600 x 3
proj = PC * u';

% plot projection
s1 = scatter3(filt2(:,1), filt2(:,2), filt2(:,3));
alpha(s1, 0.2)
hold on;
s2 = scatter3(proj(:,1), proj(:,2), proj(:,3));

% smooth with a savitsky-golay smoother
proj_smooth = smoothdata(proj,'sgolay');

% count zero crossings
zcd = dsp.ZeroCrossingDetector();
numZeroCross = cast(zcd(proj_smooth(t_inds)),"double");
tremorFrequency = (numZeroCross/2)/max(TL);

% get envelope from 25 sample moving average
env_width = 25;
env = movmax(proj_smooth(t_inds),env_width);

% use the median of the moving maximum as the estimator of the amplitude
amp = median(env);

ttl = round(tremorFrequency,1) + " Hz, " + round(2*amp,1) + " mm amplitude";

% plot
subplot(nr,nc,[5 6])
hold on
plot(t,proj,'k.')
plot(t,proj_smooth,'r')
h1 = refline(0,amp);
h2 = refline(0,-amp);
h1.Color = 0.5*[1 1 1];
h2.Color = 0.5*[1 1 1];
xlim(TL)
title(ttl)
ylabel("mm")
xlabel("seconds")
