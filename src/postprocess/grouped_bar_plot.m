function grouped_bar_plot(categories, data, series_names, y_label, title_str, filepath)
%GROUPED_BAR_PLOT  Salva em filepath um gráfico de barras agrupadas (PNG, 150 dpi).
%   data tem uma linha por categoria e uma coluna por série.
    fig = figure('Visible', 'off', 'Position', [100 100 900 500]);
    bar(data);
    set(gca, 'XTick', 1:numel(categories), 'XTickLabel', categories);
    try
        xtickangle(30);
    catch
        % xtickangle pode não existir em versões mais antigas — ignora
    end
    ylabel(y_label);
    title(title_str, 'FontWeight', 'normal');
    legend(series_names, 'Location', 'best');
    grid on;

    print(fig, filepath, '-dpng', '-r150');
    close(fig);
end
