function visualization(X, idx, k, title_text)
    % plot_clusters - Visualizes k-means clustering results
    %
    %   Inputs:
    %   X          - A matrix where each row represents a waveform in PCA space.
    %   idx        - Cluster indices assigned by k-means clustering.
    %   k          - Number of clusters.
    %   title_text - Title for the plot.

    % Define colors for each cluster
    colors = {'r', 'g', 'b', 'm'}; 
    
    % Create a new figure for the cluster plot
    figure;
    hold on;
    
    % Plot each cluster with a different color
    for j = 1:k
        scatter(X(idx == j, 1), X(idx == j, 2), 100, colors{j}, 'DisplayName', ['Cluster ' num2str(j)]);
    end
    
    % Finalize the plot
    hold off;
    xlabel('1st Principal Component', 'FontSize', 17);
    ylabel('2nd Principal Component', 'FontSize', 17);
    title(title_text, 'FontSize', 12);
    legend('Location', 'best');
    set(gca, 'FontSize', 17);
    
    % Dynamically adjust axis limits
    xlim([min(X(:, 1)), max(X(:, 1))]);
    ylim([min(X(:, 2)), max(X(:, 2))]);
end
