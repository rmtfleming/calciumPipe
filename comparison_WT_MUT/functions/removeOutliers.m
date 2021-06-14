function potential_outlier = removeOutliers(X, Y, plotFlag, thr)

stats = regstats(Y,X,'linear');

% if Cook's Distance > n/4 is a typical treshold that is used to suggest
% the presence of an outlier
potential_outlier = stats.cookd > thr/length(X);

if plotFlag
    
    % Display the index of potential outliers and graph the results
    X(potential_outlier);
    figure,
    h = scatter(X,Y, 'r.');
    hold on
    currentunits = get(gca,'Units');
    set(gca, 'Units', 'Points');
    axpos = get(gca,'Position');
    set(gca, 'Units', currentunits);
    s = 0.2;
    markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
    set(h, 'SizeData', markerWidth^2)

    h = scatter(X(potential_outlier),Y(potential_outlier), 'b.');
    markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
    set(h, 'SizeData', markerWidth^2)
    set(gca,'FontSize', 16,'FontWeight','bold')
    grid on
    
end