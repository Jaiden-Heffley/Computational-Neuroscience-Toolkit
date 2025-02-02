% main.m - Example of how to run the spike sorting analysis

% Load or generate waveforms (Example: random data for testing)
waveforms = randn(100, 50);  % Replace with your actual waveform data

% Specify the number of clusters for k-means (k)
k = 3;

% Call spike sorting function
[X, idx, k] = spike_sorting(waveforms, k);

% Visualize the clusters
visualization(X, idx, k, ['k-means Clusters for k = ' num2str(k)]);