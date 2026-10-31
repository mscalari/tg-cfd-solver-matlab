function conv_heatmap_plot(ax, iters_matrix, diverged_mask, alpha_u_vals, alpha_p_vals, panel_title, max_iter)
%CONV_HEATMAP_PLOT  Desenha em ax o mapa de iterações no plano (alpha_p, alpha_u).
%   Pontos convergidos mostram o número de iterações; os demais, um X vermelho.
%
%   ENTRADAS
%     ax            : eixo de destino
%     iters_matrix  : iterações, com linhas em alpha_u e colunas em alpha_p
%     diverged_mask : true onde não houve convergência
%     alpha_u_vals, alpha_p_vals : valores dos eixos
%     panel_title   : título do painel
%     max_iter      : limite superior da escala de cores
    axes(ax); 
    imagesc(ax, alpha_p_vals, alpha_u_vals, iters_matrix);
    set(ax, 'YDir', 'normal');
    try
        colormap(ax, parula);
    catch
        colormap(ax, jet);
    end
    cb = colorbar(ax);
    ylabel(cb, 'iterações');
    caxis(ax, [0 max_iter]);
    xlabel(ax, '\alpha_p');
    ylabel(ax, '\alpha_u');
    title(ax, panel_title, 'FontWeight', 'normal');

    hold(ax, 'on');
    for iu = 1:numel(alpha_u_vals)
        for ip = 1:numel(alpha_p_vals)
            if diverged_mask(iu, ip)
                plot(ax, alpha_p_vals(ip), alpha_u_vals(iu), 'rx', ...
                    'MarkerSize', 10, 'LineWidth', 2);
            else
                text(ax, alpha_p_vals(ip), alpha_u_vals(iu), ...
                    sprintf('%d', iters_matrix(iu, ip)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'w');
            end
        end
    end
    hold(ax, 'off');
end
