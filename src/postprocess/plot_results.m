function plot_results(c, r, algorithm, esquema_convec, pasta)
%PLOT_RESULTS  Gera e salva as figuras de um caso resolvido.
%
%   ENTRADAS
%     c, r           : caso de build_case e resultado de run_case
%     algorithm      : nome do algoritmo, usado nos títulos
%     esquema_convec : nome do esquema, usado nos títulos
%     pasta          : pasta existente onde os PNGs são gravados
%
%   Todo escoamento gera convergência, linhas de corrente, campos u e v,
%   magnitude com vetores e pressão. Couette acrescenta o perfil analítico e
%   a cavidade com tampa, o benchmark de perfis u e v.
    Nx = c.Nx;  Ny = c.Ny;  xp = c.xp;  yp = c.yp;
    [u_c, v_c] = cell_center_velocity(c, r.u, r.v);
    p_plot = r.p;
    p_plot(c.is_wall) = NaN;
    [X, Y] = meshgrid(xp, yp);
    tag = sprintf('%s, Re=%.0f, %s, %s', c.escoamento, c.Re, algorithm, esquema_convec);

    fig = figure('Name', 'Convergencia');
    semilogy(r.res_hist, 'b-', 'LineWidth', 1.5);
    grid on; box on;
    xlabel('Iteracao', 'FontSize', 12);
    ylabel('Residuo de massa (norma Frobenius)', 'FontSize', 12);
    title(['Convergencia - ' tag], 'FontSize', 11, 'FontWeight', 'normal');
    salvar(fig, pasta, 'convergencia');

    fig = figure('Name', 'Linhas de corrente');
    streamslice(X, Y, u_c', v_c', 2);
    hold on;
    formatar_campo(c, 'k', 2, ['Linhas de corrente - ' tag]);
    salvar(fig, pasta, 'linhas_corrente');

    campo(c, u_c', 'Velocidade u [m/s]', tag, pasta, 'campo_u');
    campo(c, v_c', 'Velocidade v [m/s]', tag, pasta, 'campo_v');
    campo(c, p_plot', 'Pressao', tag, pasta, 'pressao');

    fig = figure('Name', 'Magnitude');
    contourf(xp, yp, sqrt(u_c'.^2 + v_c'.^2), 20, 'LineStyle', 'none');
    colorbar; colormap('jet'); hold on;
    ix = 1:max(1, round(Nx/24)):Nx;
    iy = 1:max(1, round(Ny/12)):Ny;
    quiver(X(iy,ix), Y(iy,ix), u_c(ix,iy)', v_c(ix,iy)', 0.6, 'w', 'LineWidth', 1);
    formatar_campo(c, 'k', 2, ['Magnitude [m/s] + vetores - ' tag]);
    salvar(fig, pasta, 'magnitude_vetores');

    if c.escoamento == "CouettePoiseuille"
        perfil_couette(c, r, u_c, esquema_convec, pasta);
    elseif c.escoamento == "LidDrivenCavity"
        perfis_cavidade(c, u_c, v_c, esquema_convec, pasta);
    end
end

function campo(c, f, nome, tag, pasta, arquivo)
%CAMPO  Figura de contornos preenchidos de um campo f (Ny,Nx) nos centros dos VCs.
    fig = figure('Name', nome);
    contourf(c.xp, c.yp, f, 20, 'LineStyle', 'none');
    colorbar; colormap('jet'); hold on;
    formatar_campo(c, 'w', 1.5, sprintf('%s - %s', nome, tag));
    salvar(fig, pasta, arquivo);
end

function formatar_campo(c, cor, largura, titulo)
%FORMATAR_CAMPO  Desenha o contorno sólido e formata eixos e título da figura atual.
    for k = 1:size(c.outline, 1)
        plot(c.outline(k,1:2), c.outline(k,3:4), '-', 'Color', cor, 'LineWidth', largura);
    end
    if c.fechado, axis equal tight; else, axis tight; end
    xlabel('x  [m]', 'FontSize', 12);
    ylabel('y  [m]', 'FontSize', 12);
    title(titulo, 'FontSize', 11, 'FontWeight', 'normal', 'Interpreter', 'none');
end

function perfil_couette(c, r, u_c, esquema_convec, pasta)
%PERFIL_COUETTE  Compara o perfil u(y) em x = Lx/2 com a analítica e imprime o erro.
    y_plot = [0, c.yp, c.Ly];
    u_num  = [0; r.u(round(c.Nx/2),:)'; c.U];
    u_ref  = [0; c.u_analitica(c.yp)'; c.U];

    u_ref_grid = repmat(c.u_analitica(c.yp), c.Nx, 1);
    erro_L2 = norm(u_c - u_ref_grid, 'fro') / norm(u_ref_grid, 'fro');
    fprintf(' Erro L2 relativo: %.4e\n', erro_L2);

    fig = figure('Name', 'Perfil u');
    hold on; grid on; box on;
    plot(u_ref, y_plot, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analitica');
    plot(u_num, y_plot, 'b-',  'LineWidth', 2,   'DisplayName', esquema_convec);
    xlabel('u  [m/s]', 'FontSize', 12);
    ylabel('y  [m]',   'FontSize', 12);
    title(sprintf('Perfil u - %s\nErro L2 = %.3e', esquema_convec, erro_L2), ...
        'FontSize', 11, 'FontWeight', 'normal');
    legend('Location', 'northwest', 'FontSize', 10);
    salvar(fig, pasta, 'perfil_u');
end

function perfis_cavidade(c, u_c, v_c, esquema_convec, pasta)
%PERFIS_CAVIDADE  Compara os perfis centrais da cavidade com o benchmark tabelado.
%   A tabela é a do benchmark Acenumerics usada no projeto, com Re = 1000; em
%   outro Re a curva e o erro servem só de referência visual.
    if c.Re ~= 1000
        fprintf(' Aviso: o benchmark tabelado e de Re=1000; este caso usa Re=%.0f.\n', c.Re);
    end
    [ace_y, ace_u, ace_x, ace_v] = benchmark_cavidade();
    imid = round(c.Nx/2);
    jmid = round(c.Ny/2);

    u_ref_int  = interp1(ace_y, ace_u, c.yp/c.Ly, 'linear', 'extrap');
    u_mid_norm = u_c(imid,:) / c.U;
    erro_ref   = norm(u_mid_norm - u_ref_int) / norm(u_ref_int);
    fprintf(' Erro L2 vs Acenumerics benchmark (perfil u): %.4e\n', erro_ref);

    fig = figure('Name', 'Perfil u');
    hold on; grid on; box on;
    plot([0; u_c(imid,:)'; c.U]/c.U, [0, c.yp, c.Ly]/c.Ly, 'b-', 'LineWidth', 2, ...
        'DisplayName', esquema_convec);
    plot(ace_u, ace_y, 'ko', 'MarkerSize', 7, 'DisplayName', 'Acenumerics benchmark');
    xlabel('u / U_{lid}', 'FontSize', 12);
    ylabel('y / L', 'FontSize', 12);
    title(sprintf('Perfil u em x=L/2 - Re=%.0f\nErro L2 = %.3e', c.Re, erro_ref), ...
        'FontSize', 11, 'FontWeight', 'normal');
    legend('Location', 'northwest', 'FontSize', 10);
    salvar(fig, pasta, 'perfil_u');

    fig = figure('Name', 'Perfil v');
    hold on; grid on; box on;
    plot([0, c.xp', c.Lx]/c.Lx, [0, v_c(:,jmid)', 0]/c.U, 'r-', 'LineWidth', 2, ...
        'DisplayName', esquema_convec);
    plot(ace_x, ace_v, 'ko', 'MarkerSize', 7, 'DisplayName', 'Acenumerics benchmark');
    xlabel('x / L', 'FontSize', 12);
    ylabel('v / U_{lid}', 'FontSize', 12);
    title(sprintf('Perfil v em y=L/2 - Re=%.0f', c.Re), 'FontSize', 11, 'FontWeight', 'normal');
    legend('Location', 'northeast', 'FontSize', 10);
    salvar(fig, pasta, 'perfil_v');
end

function [ace_y, ace_u, ace_x, ace_v] = benchmark_cavidade()
%BENCHMARK_CAVIDADE  Perfis de referência da cavidade com tampa, normalizados.
    ace_y = [0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0547, 0.0625, ...
             0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000, 0.6172, ...
             0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 0.9800, ...
             0.9900, 1.0000];
    ace_u = [0.0000000, -0.0397486, -0.0759628, -0.1090870, -0.1396601, ...
            -0.1812881, -0.2023300, -0.2228955, -0.3004561, -0.3885690, ...
            -0.2803696, -0.1081999, -0.0620561,  0.0570178,  0.1886747, ...
             0.3372212,  0.4723329,  0.5169277,  0.5808359,  0.6644227, ...
             0.7070189,  0.8489396,  1.0000000];
    ace_x = [0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0625, 0.0703, ...
             0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.3500, 0.5000, ...
             0.6800, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, ...
             0.9688, 0.9700, 0.9800, 0.9900, 0.9950, 1.0000];
    ace_v = [0.0000000,  0.0709203,  0.1306210,  0.1796062,  0.2189314, ...
             0.2807056,  0.2962703,  0.3099097,  0.3330442,  0.3769189, ...
             0.3339924,  0.3253592,  0.1882246,  0.0257995, -0.1730887, ...
            -0.3202137, -0.4264545, -0.5264392, -0.4103754, -0.3553213, ...
            -0.2936869, -0.2279225, -0.2178657, -0.1359277, -0.0617528, ...
            -0.0291290,  0.0000000];
end

function salvar(fig, pasta, nome)
%SALVAR  Grava a figura como PNG em 150 dpi.
    print(fig, fullfile(pasta, [nome '.png']), '-dpng', '-r150');
end
