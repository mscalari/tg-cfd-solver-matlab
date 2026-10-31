function run_comparison(estudos, fast, pasta)
%RUN_COMPARISON  Suíte comparativa SIMPLE x SIMPLEC x PISO.
%
%   ENTRADAS
%     estudos : vetor com os estudos a rodar, entre 1, 3, 4 e 5
%               1 - precisão: esquema x algoritmo, com erros L2/Linf e CPU
%               3 - sub-relaxação: mapa alpha_u x alpha_p x algoritmo
%               4 - sensibilidade: variações de mu, G e U x algoritmo
%               5 - robustez: 4 escoamentos x algoritmo
%     fast    : true usa malhas menores e tolerâncias mais frouxas; false usa
%               os critérios próximos aos do relatório original
%     pasta   : pasta de saída; recebe mat/ e figuras/
    invalidos = setdiff(estudos, [1 3 4 5]);
    if ~isempty(invalidos)
        error('run_comparison: estudo inexistente: %s. Validos: 1 3 4 5.', mat2str(invalidos));
    end

    cfg = configuracao(fast);
    cfg.mat_dir = fullfile(pasta, 'mat');
    cfg.fig_dir = fullfile(pasta, 'figuras');
    if ~exist(cfg.mat_dir, 'dir'), mkdir(cfg.mat_dir); end
    if ~exist(cfg.fig_dir, 'dir'), mkdir(cfg.fig_dir); end

    t_suite = tic;
    fprintf('=== Suite comparativa (fast=%d) ===\n', fast);
    if ismember(1, estudos), estudo1_precisao(cfg); end
    if ismember(3, estudos), estudo3_subrelaxacao(cfg); end
    if ismember(4, estudos), estudo4_parametrico(cfg); end
    if ismember(5, estudos), estudo5_robustez(cfg); end

    t_total = toc(t_suite);
    fprintf('\n=== Suite completa em %.1f s (%.1f min) ===\n', t_total, t_total/60);
    fprintf('Resultados em: %s\n', pasta);
end

function cfg = configuracao(fast)
%CONFIGURACAO  Malhas, tolerâncias e grades de cada estudo.
    cfg.algorithms   = ["SIMPLE", "SIMPLEC", "PISO"];
    cfg.algo_labels  = cellstr(cfg.algorithms);
    cfg.esquemas7    = {'Upwind 1a ordem', 'Upwind 2a ordem', 'Diferenças Centrais', ...
                        'Power Law', 'QUICK', 'MUSCL', 'HRIC'};
    cfg.esquema_fixo = 'Upwind 1a ordem';   % estudos 3, 4 e 5
    cfg.alpha_u_base = 0.8;
    cfg.cores        = lines(numel(cfg.algorithms));

    if fast
        cfg.N1 = struct('Nx', 40, 'Ny', 16, 'tol', 1e-7, 'max_iter', 3000);
        cfg.alpha_u_vals = [0.3 0.5 0.7 0.9];
        cfg.alpha_p_vals = [0.1 0.3 0.5 0.7 0.9];
        cfg.N3 = struct('Nx', 25, 'Ny', 10, 'tol', 1e-6, 'max_iter', 1200);
        cfg.N4 = struct('Nx', 32, 'Ny', 8, 'tol', 1e-7, 'max_iter', 2000);
        cfg.N5tol = 1e-6;  cfg.N5maxit = 1200;
    else
        cfg.N1 = struct('Nx', 50, 'Ny', 20, 'tol', 1e-8, 'max_iter', 5000);
        cfg.alpha_u_vals = [0.3 0.5 0.6 0.7 0.8 0.9];
        cfg.alpha_p_vals = [0.1 0.2 0.3 0.5 0.7 0.8];
        cfg.N3 = struct('Nx', 50, 'Ny', 20, 'tol', 1e-8, 'max_iter', 3000);
        cfg.N4 = struct('Nx', 50, 'Ny', 20, 'tol', 1e-8, 'max_iter', 3000);
        cfg.N5tol = 1e-7;  cfg.N5maxit = 3000;
    end
end

function estudo1_precisao(cfg)
%ESTUDO1_PRECISAO  Erro e custo de cada esquema x algoritmo no Couette base.
    fprintf('\n--- Estudo 1: Precisao (esquema x algoritmo) ---\n');
    algorithms = cfg.algorithms;  esquemas7 = cfg.esquemas7;  N1 = cfg.N1;
    nE = numel(esquemas7);  nA = numel(algorithms);

    c = build_case("CouettePoiseuille", 'Nx', N1.Nx, 'Ny', N1.Ny, 'tol', N1.tol, 'max_iter', N1.max_iter);
    imid = round(N1.Nx/2);

    S1.errL2     = nan(nE, nA);
    S1.errLinf   = nan(nE, nA);
    S1.n_iter    = nan(nE, nA);
    S1.t_cpu     = nan(nE, nA);
    S1.converged = false(nE, nA);
    S1.u_perfil  = cell(nE, nA);   % perfil u em x = Lx/2

    for ie = 1:nE
        for ia = 1:nA
            fprintf('  [%-20s | %-7s] ', esquemas7{ie}, algorithms(ia));
            r = run_case(c, algorithms(ia), esquemas7{ie});
            [errL2, errLinf] = couette_error(c, r);

            S1.errL2(ie,ia)     = errL2;
            S1.errLinf(ie,ia)   = errLinf;
            S1.n_iter(ie,ia)    = r.n_iter;
            S1.t_cpu(ie,ia)     = r.t_cpu;
            S1.converged(ie,ia) = r.converged;
            S1.u_perfil{ie,ia}  = r.u(imid, :);

            fprintf('iter=%4d  t=%.2fs  errL2=%.3e  conv=%d\n', r.n_iter, r.t_cpu, errL2, r.converged);
        end
    end

    y_perfil1 = c.yp;
    u_ref_perfil1 = c.u_analitica(c.yp);
    save(fullfile(cfg.mat_dir, 'study1_precisao.mat'), 'S1', 'esquemas7', 'algorithms', ...
         'y_perfil1', 'u_ref_perfil1', 'N1');

    fig = figure('Visible', 'off', 'Position', [50 50 1400 900]);
    for ie = 1:nE
        subplot(3, 3, ie); hold on;
        plot(u_ref_perfil1, y_perfil1, 'k-', 'LineWidth', 2, 'DisplayName', 'Analitica');
        for ia = 1:nA
            plot(S1.u_perfil{ie,ia}, y_perfil1, '--o', 'Color', cfg.cores(ia,:), ...
                'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', char(algorithms(ia)));
        end
        xlabel('u [m/s]'); ylabel('y [m]');
        title(esquemas7{ie}, 'FontWeight', 'normal', 'Interpreter', 'none');
        legend(['Analitica', cfg.algo_labels], 'Location', 'best', 'FontSize', 7);
        grid on; hold off;
    end
    safe_sgtitle('Estudo 1 - Perfil u em x=L/2 - SIMPLE vs SIMPLEC vs PISO (caso base)');
    print(fig, fullfile(cfg.fig_dir, 'estudo1_perfil_u.png'), '-dpng', '-r150');
    close(fig);

    grouped_bar_plot(esquemas7, S1.errL2, cfg.algo_labels, 'Erro L2 relativo', ...
        'Estudo 1 - Erro L2 vs solucao analitica', fullfile(cfg.fig_dir, 'estudo1_erroL2.png'));
    grouped_bar_plot(esquemas7, S1.errLinf, cfg.algo_labels, 'Erro L_\infty relativo', ...
        'Estudo 1 - Erro L_\infty vs solucao analitica', fullfile(cfg.fig_dir, 'estudo1_erroLinf.png'));
    grouped_bar_plot(esquemas7, S1.t_cpu, cfg.algo_labels, 'Tempo de CPU [s]', ...
        'Estudo 1 - Tempo de CPU', fullfile(cfg.fig_dir, 'estudo1_tempoCPU.png'));

    fprintf('Estudo 1 concluido.\n');
end

function estudo3_subrelaxacao(cfg)
%ESTUDO3_SUBRELAXACAO  Convergência no espaço (alpha_u, alpha_p) por algoritmo.
    fprintf('\n--- Estudo 3: Sub-relaxacao (heatmap alpha_u x alpha_p) ---\n');
    algorithms = cfg.algorithms;  N3 = cfg.N3;
    alpha_u_vals = cfg.alpha_u_vals;  alpha_p_vals = cfg.alpha_p_vals;
    nA = numel(algorithms);  nU = numel(alpha_u_vals);  nP = numel(alpha_p_vals);

    S3.iters     = nan(nU, nP, nA);
    S3.t_cpu     = nan(nU, nP, nA);
    S3.converged = false(nU, nP, nA);

    caso = @(au, ap) build_case("CouettePoiseuille", 'Nx', N3.Nx, 'Ny', N3.Ny, ...
        'tol', N3.tol, 'max_iter', N3.max_iter, 'alpha_u', au, 'alpha_p', ap);

    for ia = 1:nA
        for iu = 1:nU
            for ip = 1:nP
                r = run_case(caso(alpha_u_vals(iu), alpha_p_vals(ip)), algorithms(ia), cfg.esquema_fixo);
                S3.iters(iu,ip,ia)     = r.n_iter;
                S3.t_cpu(iu,ip,ia)     = r.t_cpu;
                S3.converged(iu,ip,ia) = r.converged;
            end
        end
        fprintf('  Algoritmo %-7s: grade %dx%d concluida\n', algorithms(ia), nU, nP);
    end

    save(fullfile(cfg.mat_dir, 'study3_subrelaxacao.mat'), 'S3', 'algorithms', ...
         'alpha_u_vals', 'alpha_p_vals', 'N3');

    fig = figure('Visible', 'off', 'Position', [50 50 1500 450]);
    for ia = 1:nA
        ax = subplot(1, nA, ia);
        iters_plot = S3.iters(:,:,ia);
        diverged   = ~S3.converged(:,:,ia);
        iters_plot(diverged) = N3.max_iter;
        conv_heatmap_plot(ax, iters_plot, diverged, alpha_u_vals, alpha_p_vals, ...
            char(algorithms(ia)), N3.max_iter);
    end
    safe_sgtitle('Estudo 3 - Mapa de convergencia no espaco (\alpha_u, \alpha_p)');
    print(fig, fullfile(cfg.fig_dir, 'estudo3_heatmap_subrelaxacao.png'), '-dpng', '-r150');
    close(fig);

    % SIMPLE e SIMPLEC: menor n_iter entre os convergidos, na linha de alpha_u
    % mais próxima da base. PISO: alpha_p = 1, sem relaxação de pressão.
    [~, iu_base] = min(abs(alpha_u_vals - cfg.alpha_u_base));
    best_ap = nan(1,2); best_iter = nan(1,3); best_time = nan(1,3); best_label = cell(1,3);
    for ia = 1:2
        row_iter = squeeze(S3.iters(iu_base, :, ia));
        row_conv = squeeze(S3.converged(iu_base, :, ia));
        row_iter(~row_conv) = Inf;
        [best_iter(ia), idx_best] = min(row_iter);
        best_ap(ia) = alpha_p_vals(idx_best);
        best_time(ia) = S3.t_cpu(iu_base, idx_best, ia);
        best_label{ia} = sprintf('%s (\\alpha_p=%.1f)', algorithms(ia), best_ap(ia));
    end

    r_piso = run_case(caso(alpha_u_vals(iu_base), 1.0), "PISO", cfg.esquema_fixo);
    best_iter(3) = r_piso.n_iter;
    best_time(3) = r_piso.t_cpu;
    best_label{3} = 'PISO (\alpha_p=1.0)';

    fig = figure('Visible', 'off', 'Position', [100 100 900 450]);
    subplot(1,2,1);
    bar(best_iter); set(gca, 'XTick', 1:3, 'XTickLabel', best_label);
    ylabel('iteracoes ate convergir'); title('Iteracoes no melhor \alpha_p', 'FontWeight', 'normal');
    grid on;
    subplot(1,2,2);
    bar(best_time); set(gca, 'XTick', 1:3, 'XTickLabel', best_label);
    ylabel('tempo de CPU [s]'); title('Tempo no melhor \alpha_p', 'FontWeight', 'normal');
    grid on;
    safe_sgtitle(sprintf('Estudo 3 - Melhor configuracao (\\alpha_u=%.1f)', alpha_u_vals(iu_base)));
    print(fig, fullfile(cfg.fig_dir, 'estudo3_melhor_alpha_p.png'), '-dpng', '-r150');
    close(fig);

    save(fullfile(cfg.mat_dir, 'study3_melhor_alpha_p.mat'), 'best_ap', 'best_iter', 'best_time', 'best_label');
    fprintf('Estudo 3 concluido.\n');
end

function estudo4_parametrico(cfg)
%ESTUDO4_PARAMETRICO  Sensibilidade do Couette a mu, G >= 0, G < 0 e U.
    fprintf('\n--- Estudo 4: Sensibilidade parametrica (Blocos A-D) ---\n');
    algorithms = cfg.algorithms;  N4 = cfg.N4;
    nA = numel(algorithms);
    Ly0 = 1; mu0 = 0.1; rho0 = 1; U0 = 2; G0 = 4;

    % Bloco A - variação de mu (Pe de célula alvo, dy fixo)
    Pe_alvo = [0.1 0.2 0.5 1 2];
    blocos.A = bloco('Bloco A - variacao de \mu', 'Pe_{celula}', Pe_alvo, ...
        rho0*U0*(Ly0/N4.Ny) ./ Pe_alvo, G0*ones(size(Pe_alvo)), U0*ones(size(Pe_alvo)));

    % Bloco B - variação de G >= 0 (parâmetro P alvo, mu e U fixos)
    P_alvo = [0 2.5 5 10 20 40];
    blocos.B = bloco('Bloco B - variacao de G (>=0)', 'P_{param}', P_alvo, ...
        mu0*ones(size(P_alvo)), P_alvo * 2*mu0*U0 / Ly0^2, U0*ones(size(P_alvo)));

    % Bloco C - gradiente adverso
    absG = [0.5 1 2 4 8];
    blocos.C = bloco('Bloco C - gradiente adverso (G<0)', '|G| [Pa/m]', absG, ...
        mu0*ones(size(absG)), -absG, U0*ones(size(absG)));

    % Bloco D - variação de U (G e mu fixos)
    U_alvo = [0.25 0.5 1 2 4 8];
    blocos.D = bloco('Bloco D - variacao de U', 'U [m/s]', U_alvo, ...
        mu0*ones(size(U_alvo)), G0*ones(size(U_alvo)), U_alvo);

    S4 = struct();
    for nb = ["A", "B", "C", "D"]
        b = blocos.(nb);
        npts = numel(b.x);
        errL2_b = nan(npts, nA);
        niter_b = nan(npts, nA);
        conv_b  = false(npts, nA);

        for ipt = 1:npts
            c = build_case("CouettePoiseuille", 'Nx', N4.Nx, 'Ny', N4.Ny, 'tol', N4.tol, ...
                'max_iter', N4.max_iter, 'mu', b.mu(ipt), 'G', b.G(ipt), 'U', b.U(ipt));
            for ia = 1:nA
                r = run_case(c, algorithms(ia), cfg.esquema_fixo);
                errL2_b(ipt,ia) = couette_error(c, r);
                niter_b(ipt,ia) = r.n_iter;
                conv_b(ipt,ia)  = r.converged;
            end
        end

        S4.(nb).errL2 = errL2_b;
        S4.(nb).n_iter = niter_b;
        S4.(nb).converged = conv_b;
        fprintf('  %s concluido (%d pontos)\n', b.label, npts);

        fig = figure('Visible', 'off', 'Position', [100 100 1100 450]);
        painel_bloco(1, b, errL2_b, cfg, 'Erro L2 relativo', 'Erro L2');
        painel_bloco(2, b, niter_b, cfg, 'Iteracoes ate convergir', 'Iteracoes');
        safe_sgtitle(sprintf('Estudo 4 - %s', b.label));
        print(fig, fullfile(cfg.fig_dir, sprintf('estudo4_bloco%s.png', nb)), '-dpng', '-r150');
        close(fig);
    end

    save(fullfile(cfg.mat_dir, 'study4_parametrico.mat'), 'S4', 'blocos', 'algorithms', 'N4');
    fprintf('Estudo 4 concluido.\n');
end

function b = bloco(label, x_label, x, mu, G, U)
%BLOCO  Agrupa rótulos e valores de um bloco de sensibilidade do Estudo 4.
    b = struct('label', label, 'xlabel', x_label, 'x', x, 'mu', mu, 'G', G, 'U', U);
end

function painel_bloco(pos, b, dados, cfg, y_label, titulo)
%PAINEL_BLOCO  Subplot com uma curva por algoritmo de um bloco do Estudo 4.
    subplot(1, 2, pos); hold on;
    for ia = 1:numel(cfg.algorithms)
        plot(b.x, dados(:,ia), '-o', 'Color', cfg.cores(ia,:), ...
            'LineWidth', 1.5, 'DisplayName', char(cfg.algorithms(ia)));
    end
    xlabel(b.xlabel); ylabel(y_label);
    title(titulo, 'FontWeight', 'normal');
    legend('Location', 'best'); grid on; hold off;
end

function estudo5_robustez(cfg)
%ESTUDO5_ROBUSTEZ  Iterações, tempo e campo de |V| nos 4 escoamentos.
    fprintf('\n--- Estudo 5: Robustez nos 4 escoamentos ---\n');
    algorithms = cfg.algorithms;
    nA = numel(algorithms);

    escoamentos5 = ["CouettePoiseuille", "LidDrivenCavity", "FlowOverCavity", "FlowOverStep"];
    parametros = { ...
        {'Nx', 32, 'Ny', 12}, ...
        {'Nx', 16, 'Ny', 16, 'Re', 100}, ...
        {'Nx', 24, 'Ny', 12, 'Re', 100, 'L_cav', 0.5, 'D_cav', 0.5, ...
         'Ly_canal', 0.5, 'Lx_ante', 0.5, 'Lx_post', 1.0}, ...
        {'Nx', 20, 'Ny', 10, 'Re', 50, 'Lx_post', 1.5}};
    nF = numel(escoamentos5);

    S5.n_iter    = nan(nF, nA);
    S5.t_cpu     = nan(nF, nA);
    S5.converged = false(nF, nA);
    S5.campo     = cell(nF, nA);

    for ifl = 1:nF
        c = build_case(escoamentos5(ifl), parametros{ifl}{:}, 'tol', cfg.N5tol, 'max_iter', cfg.N5maxit);
        for ia = 1:nA
            fprintf('  [%-17s | %-7s] ', escoamentos5(ifl), algorithms(ia));
            r = run_case(c, algorithms(ia), cfg.esquema_fixo);
            [u_c, v_c] = cell_center_velocity(c, r.u, r.v);
            S5.n_iter(ifl,ia)    = r.n_iter;
            S5.t_cpu(ifl,ia)     = r.t_cpu;
            S5.converged(ifl,ia) = r.converged;
            S5.campo{ifl,ia}     = struct('u', r.u, 'v', r.v, 'p', r.p, ...
                'umag', sqrt(u_c.^2 + v_c.^2), 'xp', c.xp, 'yp', c.yp, 'Nx', c.Nx, 'Ny', c.Ny);
            fprintf('iter=%4d  t=%.2fs  conv=%d\n', r.n_iter, r.t_cpu, r.converged);
        end
    end

    esc_labels5 = cellstr(escoamentos5);
    save(fullfile(cfg.mat_dir, 'study5_robustez.mat'), 'S5', 'escoamentos5', 'esc_labels5', 'algorithms');

    fig = figure('Visible', 'off', 'Position', [100 100 1300 450]);
    subplot(1,2,1);
    bar(S5.n_iter); set(gca, 'XTick', 1:nF, 'XTickLabel', esc_labels5);
    ylabel('iteracoes ate convergir'); title('Iteracoes por escoamento', 'FontWeight', 'normal');
    legend(cfg.algo_labels, 'Location', 'best'); grid on;
    subplot(1,2,2);
    bar(S5.t_cpu); set(gca, 'XTick', 1:nF, 'XTickLabel', esc_labels5);
    ylabel('tempo de CPU [s]'); title('Tempo de CPU por escoamento', 'FontWeight', 'normal');
    legend(cfg.algo_labels, 'Location', 'best'); grid on;
    safe_sgtitle('Estudo 5 - Robustez nos 4 escoamentos');
    print(fig, fullfile(cfg.fig_dir, 'estudo5_iteracoes_tempo.png'), '-dpng', '-r150');
    close(fig);

    for ifl = 1:nF
        fig = figure('Visible', 'off', 'Position', [50 50 1500 450]);
        for ia = 1:nA
            campo = S5.campo{ifl,ia};
            subplot(1, nA, ia);
            imagesc(campo.xp, campo.yp, campo.umag'); set(gca, 'YDir', 'normal'); colorbar;
            xlabel('x [m]'); ylabel('y [m]');
            title(char(algorithms(ia)), 'FontWeight', 'normal');
        end
        safe_sgtitle(sprintf('Estudo 5 - |V| - %s', escoamentos5(ifl)));
        print(fig, fullfile(cfg.fig_dir, sprintf('estudo5_campo_%s.png', escoamentos5(ifl))), '-dpng', '-r150');
        close(fig);
    end

    fprintf('Estudo 5 concluido.\n');
end
