function [X, idx, k] = spike_sorting(waveforms, k)
    % spike_sorting - Perform PCA and k-means clustering on neural waveforms
    %
    %   Inputs:
    %   waveforms - A matrix of neural waveforms, where each row represents
    %               a waveform and each column represents a data point.
    %   k         - Number of clusters for k-means.
    %
    %   Outputs:
    %   X         - The projected waveforms in PCA space.
    %   idx       - The cluster indices from k-means.
    %   k         - The number of clusters.

    % Perform PCA
    [coef, score, latent] = pca(waveforms);
    
    % Calculate the cumulative explained variance
    cumulative_variance = cumsum(latent) / sum(latent);
    
    % Determine how many components to keep based on the variance threshold
    variance_threshold = 0.95;  % Default threshold to keep components that explain 95% variance
    num_components = find(cumulative_variance >= variance_threshold, 1);
    
    % Project the waveforms onto the selected principal components
    X = score(:, 1:num_components);
    
    % Perform k-means clustering
    [idx, ~] = kmeans(X, k);
end
