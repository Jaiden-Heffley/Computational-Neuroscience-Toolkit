% Initialize dataset
data = [5.5,3.1; 5.1,4.8; 6.4,3.6; 5.6,4.7; 6.7,3.7];

% Initial cluster centers (as per the problem)
initial_centers = [5.3, 3.5; 5.1, 4.2; 6.0, 3.9];

figure;
subplot(1,2,1)
scatter(data(:,1),data(:,2),100,'filled')
hold on
scatter(initial_centers(:,1),initial_centers(:,2),200, 'kx', 'LineWidth', 2)
% Add plot formatting
xlim([5, 7]);
ylim([3, 5]);
xticks(5:0.5:7);
yticks(3:0.5:5);
xticklabels({'', '5.5', '6.0', '6.5', ''});
yticklabels({'', '3.5', '4', '4.5', ''});
title('Original Data','FontSize',20);
hold off

set(findall(gcf,'type','axes'),'FontSize',17);

% Apply k-means and let it converge
[idx_converged, final_centers] = kmeans(data, 3, 'Start', initial_centers);

% Plot the points and their cluster assignments after convergence
subplot(1,2,2)
scatter(data(:,1).*(idx_converged == 1), data(:,2).*(idx_converged == 1), 100, 'r', 'filled', 'DisplayName', 'Cluster 1'); % Cluster 1 points
hold on;
scatter(data(:,1).*(idx_converged == 2), data(:,2).*(idx_converged == 2), 100, 'g', 'filled', 'DisplayName', 'Cluster 2'); % Cluster 2 points
scatter(data(:,1).*(idx_converged == 3), data(:,2).*(idx_converged == 3), 100, 'b', 'filled', 'DisplayName', 'Cluster 3'); % Cluster 3 points

% Plot final cluster centers
scatter(final_centers(:,1), final_centers(:,2), 200, 'kx', 'LineWidth', 2, 'DisplayName', 'Final Centers');

% Add plot formatting
xlim([5, 7]);
ylim([3, 5]);
xticks(5:0.5:7);
yticks(3:0.5:5);
xticklabels({'', '5.5', '6.0', '6.5', ''});
yticklabels({'', '3.5', '4', '4.5', ''});
legend('Location', 'best'); % Automatically handles all 'DisplayName' labels
title('K-Means Clustering - After Convergence','FontSize',20);
hold off;
set(findall(gcf,'type','axes'),'FontSize',17);

% Display the final cluster centers
fprintf('Final cluster center for Cluster 1: (%.2f, %.2f)\n', final_centers(1,1), final_centers(1,2));
fprintf('Final cluster center for Cluster 2: (%.2f, %.2f)\n', final_centers(2,1), final_centers(2,2));
fprintf('Final cluster center for Cluster 3: (%.2f, %.2f)\n', final_centers(3,1), final_centers(3,2));

%% Spike Sorting
homeDir = getenv('USERPROFILE');
load(fullfile(homeDir,'OneDrive - University of Pittsburgh\Documents\MATLAB\BIOENG 1586\SpikeSorting\waveforms.mat'))
waveforms = data.wf;
timestamps = data.stamps;
% Number of waveforms to plot
num_waveforms = 100;

% Randomly sample indices to spread waveforms across the session
total_waveforms = size(waveforms, 1); % Total number of waveforms
indices = round(linspace(1, total_waveforms, num_waveforms)); % Spread-out indices
sampled_waveforms = waveforms(indices, :);

% Plot the waveforms
figure;
plot(sampled_waveforms');
xlabel('Time (s)','FontSize',17);
ylabel('Amplitude','FontSize',17);
title('100 Sampled Waveforms from Neural Recording','FontSize',20);
grid on;
set(findall(gcf,'type','axes'),'FontSize',17);

% Brief narrative: Visual inspection for neuron estimation
fprintf('Q1: 100 waveforms plotted. Judging by the waveform shapes, I estimate ~2-4 distinct neurons (clusters) in the data.\n\n');

%% PCA
% Q1) I would say line d corresponds withe the first principle component
% Q2) I woulkd say line c corresponds with the sceond principle component
% Q3) I would say line d corresponds with the third prinsiple component
% Q4) Total variance = 2.2 + 1.7 + 1.4 + 0.8 + 0.4 + 0.2 + 0.15 + 0.02 +
% 0.001 = 6.871
% (2.2 + 1.7 + 1.4 + 0.8 + 0.4)/6.871 = 0.946
% Therefore, to retain at least 90% variance, you need 5 components

%% Perform PCA on Spike set
%perform PCA
[coef, score, latent] = pca(waveforms);

%find  projected waveforms
wf_PC1 = waveforms*coef(:,1);
wf_PC2 = waveforms*coef(:,2);

%choosing inital k's
X = [wf_PC1 wf_PC2];

% Trying different values of k
k_values = [2, 3, 4];  % Different values of k to test
colors = {'r', 'g', 'b', 'm'};  % Color list for clusters: red, green, blue, magenta
figure;  % Create a new figure for the subplots

for i = 1:length(k_values)
    k = k_values(i);
    
    % Perform k-means clustering
    [idx, cluster_centers] = kmeans(X, k);

    % Plot the k-means clusters in the corresponding subplot
    subplot(1, 3, i);  % Create subplot for the current k value
    hold on;
    for j = 1:k
        % Use the specified colors for the clusters
        scatter(X(idx == j, 1), X(idx == j, 2), 100, colors{j}, 'DisplayName', ['Cluster ' num2str(j)]);
    end
    hold off;
    title(['k-means Clusters for k = ' num2str(k)], 'FontSize', 12);
    xlabel('1st PC', 'FontSize', 17);
    ylabel('2nd PC', 'FontSize', 17);
    legend('Location', 'best');
    set(findall(gcf,'type','axes'),'FontSize',17);
end

%plot first two pcs in 2D
figure;
subplot(1,3,1)
scatter(wf_PC1,wf_PC2);
xlabel('1st PC','FontSize', 17);
ylabel('2nd PC','FontSize',17);
title('PCA: Scatter of First Two Components','FontSize', 10);
set(findall(gcf,'type','axes'),'FontSize',17);

fprintf(['Q2: It looks like there may be 2 or 3 neurons present in the recording, \n ' ...
    'but it is hard to tell still. \n I would defintely say there isnt 4. \n \n']);

% Plot the first two PCs with color mapping for the event index
subplot(1,3,2)
scatter3(wf_PC1, score(:,2), (1:length(score))', 10, (1:length(score))', 'filled');
colormap jet;
colorbar;
xlabel('1st PC','FontSize',17);
ylabel('2nd PC','FontSize',17);
zlabel('Event Index','FontSize',17);
title('PCA: 3D Scatter of First Two Components','FontSize', 10);
set(findall(gcf,'type','axes'),'FontSize',17);
%looks like there are 3 clusters, (e.g. cluter centers 
% cls = [-184.856, 155.988; 77.1142, 150.671; -36.4561, -86.0849]

k=3;
[idx, clsuter_centers] = kmeans(X,k);

% Plot the points and their cluster assignments after convergence
subplot(1,3,3)
scatter(X(:,1).*(idx == 1), X(:,2).*(idx == 1), 100, 'r', 'DisplayName', 'Cluster 1'); % Cluster 1 points
hold on;
scatter(X(:,1).*(idx == 2), X(:,2).*(idx == 2), 100, 'g', 'DisplayName', 'Cluster 2'); % Cluster 2 points
scatter(X(:,1).*(idx == 3), X(:,2).*(idx == 3), 100, 'b', 'DisplayName', 'Cluster 3'); % Cluster 3 points
hold off
title("K-means Clusters",'FontSize', 10)
xlabel('1st PC','FontSize',17);
ylabel('2nd PC','FontSize', 17);
set(findall(gcf,'type','axes'),'FontSize',17);

%Color coding 100 waveforms by cluster assignmetns
sampled_idx=idx(indices); %sampling indexes to be same as sampled wavforms

%plot the waveforms in each cluster
figure;
h1 = plot(sampled_waveforms(sampled_idx == 1, :)', 'LineWidth', 1.5, 'Color','r');
hold on
h2 = plot(sampled_waveforms(sampled_idx == 2, :)', 'LineWidth', 1.5, 'Color','g');
h3 = plot(sampled_waveforms(sampled_idx == 3, :)', 'LineWidth', 1.5, 'Color','b');
legend([h1(1), h2(1), h3(1)], {'Cluster 1', 'Cluster 2', 'Cluster 3'}, 'Location', 'best');
hold off
set(findall(gcf,'type','axes'),'FontSize',17);

ylabel('Amplitude', 'FontSize', 17)
xlabel('Time (s)', 'FontSize', 17)
title('Waveforms Colored by PCA Cluster')

fprintf(['Q4: It looks like it worked. The green spikes all have the similar \n ' ...
    'pattern of spiking early, the red spikes \n have the pattern of having lower spikes \n ' ...
    'and the blue spikes have later and \n more drawn out spikes. \n \n '])

%% K-means with initalzied cluster centers
cls = [0, 0; 0, 0; 0, 0];
k=3;
[idx_origin, clsuter_centers_origin] = kmeans(X,k, 'Start',cls);

% Plot the points and their cluster assignments after convergence
figure;
subplot(1,3,1)
scatter(X(:,1).*(idx_origin == 1), X(:,2).*(idx_origin == 1), 100, 'r', 'DisplayName', 'Cluster 1'); % Cluster 1 points
hold on;
scatter(X(:,1).*(idx_origin == 2), X(:,2).*(idx_origin == 2), 100, 'g', 'DisplayName', 'Cluster 2'); % Cluster 2 points
scatter(X(:,1).*(idx_origin == 3), X(:,2).*(idx_origin == 3), 100, 'b', 'DisplayName', 'Cluster 3'); % Cluster 3 points
hold off
title("K-Means Clusters with Initial Centers",'FontSize', 12)
ylabel('PC 1','FontSize', 17);
xlabel('PC 2','FontSize', 17);
set(findall(gcf,'type','axes'),'FontSize',17);


cls = [-184.856, 155.988; 77.1142, 150.671; -36.4561, -86.0849];
k=3;
[idx_est , clsuter_centers_est] = kmeans(X,k, 'Start',cls);

% Plot the points and their cluster assignments after convergence
subplot(1,3,2)
scatter(X(:,1).*(idx_est == 1), X(:,2).*(idx_est == 1), 100, 'r', 'DisplayName', 'Cluster 1'); % Cluster 1 points
hold on;
scatter(X(:,1).*(idx_est == 2), X(:,2).*(idx_est == 2), 100, 'g', 'DisplayName', 'Cluster 2'); % Cluster 2 points
scatter(X(:,1).*(idx_est == 3), X(:,2).*(idx_est == 3), 100, 'b', 'DisplayName', 'Cluster 3'); % Cluster 3 points
hold off
title("K-Means Clusters with Estimated Centers",'FontSize', 15);
ylabel('PC 1','FontSize', 17);
xlabel('PC 2','FontSize', 17);
set(findall(gcf,'type','axes'),'FontSize',17);


cls = [-380, 267.711; -386, 46; -363,-129];
k=3;
[idx_off, clsuter_centers_off] = kmeans(X,k, 'Start',cls);

% Plot the points and their cluster assignments after convergence
subplot(1,3,3)
scatter(X(:,1).*(idx_off == 1), X(:,2).*(idx_off == 1), 100, 'r', 'DisplayName', 'Cluster 1'); % Cluster 1 points
hold on;
scatter(X(:,1).*(idx_off == 2), X(:,2).*(idx_off == 2), 100, 'g', 'DisplayName', 'Cluster 2'); % Cluster 2 points
scatter(X(:,1).*(idx_off == 3), X(:,2).*(idx_off == 3), 100, 'b', 'DisplayName', 'Cluster 3'); % Cluster 3 points
hold off
title("K-Means Clusters with Offset Centers",'FontSize', 15)
ylabel('PC 1','FontSize', 17);
xlabel('PC 2','FontSize', 17);
set(findall(gcf,'type','axes'),'FontSize',17);

fprintf("Q5: It seems like \n \n ")


%% Unit vector of 3 PCS
first_3_pcs = coef(:,1:3);

%Reconstruct the original waveform for each of the first three PCs
reconstructed_first_pc = first_3_pcs(:,1);  % First PC
reconstructed_second_pc = first_3_pcs(:,2); % Second PC
reconstructed_third_pc = first_3_pcs(:,3);  % Third PC

% Plot the waveforms for the first three principal components
figure;
subplot(3,1,1);
plot(reconstructed_first_pc);
title('Waveform of First Principal Component','FontSize', 20);

subplot(3,1,2);
plot(reconstructed_second_pc);
ylabel('Amplitude','FontSize',17)
title('Waveform of Second Principal Component','FontSize', 20);

subplot(3,1,3);
plot(reconstructed_third_pc);
xlabel('Time (s)','FontSize', 17);
title('Waveform of Third Principal Component','FontSize',20);

set(findall(gcf,'type','axes'),'FontSize',17);

fprintf("Q6: It looks like one of the general spike shapes in the initial data. \n" + ...
    " It seems like PCA splits the data into similar shape catgegories \n " + ...
    "and labels the data based on their shape. \n \n")

%% Calculaing explained variance 

% Calculate the variance explained by each principal component
variance_explained = latent / sum(latent);  % Normalize by total variance
cumulative_variance = cumsum(variance_explained);  % Cumulative variance

% Plot the eigenspectrum (variance explained by each PC)
figure;
plot(1:length(cumulative_variance), cumulative_variance, '-o');
xlabel('Number of Dimensions (Principal Components)','FontSize', 17);
ylabel('Cumulative Variance Explained', 'FontSize',17);
title('Eigenspectrum of the Data','FontSize', 20);
grid on;
set(findall(gcf,'type','axes'),'FontSize',17);

num_dimensions = size(waveforms, 2);  % Number of columns (features)

num_components_100 = find(cumulative_variance >= 1, 1);  % First PC where cumulative variance reaches 100%

num_components_90 = find(cumulative_variance >= 0.9, 1);  % First PC where cumulative variance reaches 90%

fprintf("Q7: There are %d dimensions in the original data \n", num_dimensions)
fprintf("It takes %d dimensions to account for 100%% of the variance in the waveform data \n", num_components_100)
fprintf("It takes %d dimensions to account for 90%% of the variance in the waveform data \n", num_components_90)
fprintf("I'd estimate that the elbow is at about 11 dimensions \n")

%% Rasters with Color-Coded Cells
% Define the 10s time window 
start_time = 20;  % Starting time for the window
end_time = 30;   % Ending time for the window

% Filter timestamps within the time window
valid_idx = (timestamps >= start_time & timestamps <= end_time);
filtered_stamps = timestamps(valid_idx);
filtered_cells = idx(valid_idx);

% Get unique cell IDs and assign colors
unique_cells = unique(filtered_cells);
colors = lines(length(unique_cells)); % Generate distinct colors

% Plot raster
figure;
hold on;

% Loop through each cell and plot its spikes
for i = 1:length(unique_cells)
    cell_id = unique_cells(i);
    cell_timestamps = filtered_stamps(filtered_cells == cell_id);
    
    % Plot spikes for this cell on a separate line with a unique color
    for j = 1:length(cell_timestamps)
        line([cell_timestamps(j), cell_timestamps(j)], [i - 0.4, i + 0.4], ...
            'Color', colors(i, :), 'LineWidth', 1.5);
    end
end

% Customize plot
xlabel('Time (s)','FontSize', 17);
ylabel('Cell ID', 'FontSize', 17);
title(sprintf('Raster Plot for Time Window %d to %d Seconds', start_time, end_time),'FontSize', 20);
ylim([0.5, length(unique_cells) + 0.5]); % Adjust y-limits to fit all cells
yticks(1:length(unique_cells));
yticklabels(arrayfun(@(x) sprintf('Cell %d', x), unique_cells, 'UniformOutput', false));
grid on;
hold off;
set(findall(gcf,'type','axes'),'FontSize',17);


%% Interspike Interval Histogram (ISI)
num_neurons = max(idx);  % Number of clusters (neurons)
timestamps = 1:length(waveforms);  % Adjust this based on your timestamp data

figure;
for neuron = 1:num_neurons
    % Get the indices of the points that belong to this neuron (cluster)
    neuron_spikes_idx = timestamps(idx == neuron);
    
    % Calculate ISIs for this neuron
    ISI = diff(neuron_spikes_idx);  % Time differences between consecutive spikes
    
    %Plot the ISI histogram for this neuron
    subplot(1,3,neuron)
    histogram(ISI, 'EdgeColor', 'black', 'FaceColor', [0.8, 0.8, 0.8]);
    xlabel('Interspike Interval (s)','FontSize',17);
    ylabel('Frequency','FontSize',17);
    title(['ISI Histogram for Neuron ' num2str(neuron)],'FontSize', 20);
    set(findall(gcf,'type','axes'),'FontSize',17);
    xlim([0,40])
end
