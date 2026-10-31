%% RUN_COMPARACAO_PARTE1
%
%   Suíte comparativa SIMPLE x SIMPLEC x PISO — Parte 1.
%
%   Roda 4 estudos e salva tudo em ./Comparacao_parte_1/{mat,figuras}:
%
%     1) Precisão            : esquema_convec (7) x algorithm (3), caso base
%                               Couette-Poiseuille — perfil u(y), erro L2/Linf,
%                               tempo de CPU.
%     3) Sub-relaxação       : alpha_u x alpha_p x algorithm — heatmap de
%                               convergência + comparação no melhor alpha_p
%                               de cada algoritmo (PISO sem relaxação de p).
%     4) Sensibilidade       : Blocos A-D (mu, G>=0, G<0, U) x algorithm,
%        paramétrica           esquema fixo.
%     5) Robustez            : 4 escoamentos x algorithm, esquema fixo —
%                               iterações/tempo + campos de velocidade.
%
%   CUSTO COMPUTACIONAL: este script foi desenhado para rodar o mais
%   rápido possível mantendo resultados úteis. O flag FAST_MODE abaixo
%   controla malha/tolerância/nº de pontos de cada estudo:
%     true  (default) -> malhas menores, tolerância mais frouxa, grades
%                         de sub-relaxação mais esparsas. Pensado para
%                         iteração rápida durante a análise.
%     false            -> malhas e critérios próximos aos do relatório
%                         original (mais lento, mais preciso/publicável).
%   Nenhum arquivo de solver (simple_solver.m, simplec_solver.m,
%   piso_solver.m, pressure_correc*.m) é modificado por este script.
%
%   Pré-requisito: rodar a partir da pasta que contém todas as funções
%   do solver (mesma pasta de main.m).

clear; clc; close all;

FAST_MODE = false;   % true = rápido (default) | false = malhas maiores

% ────────────────────────────────────────────────────────────────────
%  ESTRUTURA DE PASTAS
% ────────────────────────────────────────────────────────────────────
base_dir = 'Comparacao_parte_1';
mat_dir  = fullfile(base_dir, 'mat');
fig_dir  = fullfile(base_dir, 'figuras');
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

algorithms = {"SIMPLE", "SIMPLEC", "PISO"};
algo_labels = cellfun(@char, algorithms, 'UniformOutput', false);   % rótulos p/ legend()
esquemas7  = {'Upwind 1a ordem', 'Upwind 2a ordem', 'Diferenças Centrais', ...
              'Power Law', 'QUICK', 'MUSCL', 'HRIC'};
esquema_fixo = 'Upwind 1a ordem';   % usado nos Estudos 3, 4 e 5
interp_pressao = "Linear";

% Caso base Couette-Poiseuille (mesmo do relatório, Seção 4.1)
Lx0 = 5; Ly0 = 1; mu0 = 0.1; rho0 = 1; U0 = 2; G0 = 4;

if FAST_MODE
    % Estudo 1 — precisão
    N1 = struct('Nx', 40, 'Ny', 16, 'tol', 1e-7, 'max_iter', 3000);
    % Estudo 3 — sub-relaxação
    alpha_u_vals = [0.3 0.5 0.7 0.9];
    alpha_p_vals = [0.1 0.3 0.5 0.7 0.9];
    N3 = struct('Nx', 25, 'Ny', 10, 'tol', 1e-6, 'max_iter', 1200);
    % Estudo 4 — sensibilidade paramétrica
    N4 = struct('Nx', 32, 'Ny', 8, 'tol', 1e-7, 'max_iter', 2000);
    % Estudo 5 — robustez nos 4 escoamentos
    N5tol = 1e-6; N5maxit = 1200;
else
    N1 = struct('Nx', 50, 'Ny', 20, 'tol', 1e-8, 'max_iter', 5000);
    alpha_u_vals = [0.3 0.5 0.6 0.7 0.8 0.9];
    alpha_p_vals = [0.1 0.2 0.3 0.5 0.7 0.8];
    N3 = struct('Nx', 50, 'Ny', 20, 'tol', 1e-8, 'max_iter', 3000);
    N4 = struct('Nx', 50, 'Ny', 20, 'tol', 1e-8, 'max_iter', 3000);
    N5tol = 1e-7; N5maxit = 3000;
end

alpha_u_base = 0.8; alpha_p_base = 0.8;

t_suite = tic;
fprintf('=== Suite comparativa Parte 1 (FAST_MODE=%d) ===\n', FAST_MODE);

%% ════════════════════════════════════════════════════════════════════
%  ESTUDO 1 — PRECISÃO (esquema x algoritmo)
%  ════════════════════════════════════════════════════════════════════
fprintf('\n--- Estudo 1: Precisao (esquema x algoritmo) ---\n');

[numerical_param1, physical_param1, D1, u_an1, xp1, yp1] = build_couette_case( ...
    N1.Nx, N1.Ny, Lx0, Ly0, mu0, rho0, U0, G0, alpha_u_base, alpha_p_base, N1.tol, N1.max_iter);

nE = numel(esquemas7); nA = numel(algorithms);
S1.errL2    = nan(nE, nA);
S1.errLinf  = nan(nE, nA);
S1.n_iter   = nan(nE, nA);
S1.t_cpu    = nan(nE, nA);
S1.converged= false(nE, nA);
S1.u_perfil = cell(nE, nA);     % perfil u em x=Lx/2 (para plot)

% Partida do repouso: começar da solução analítica falsearia as iterações
u0_1 = zeros(N1.Nx-1, N1.Ny);
v0_1 = zeros(N1.Nx, N1.Ny-1);
p0_1 = zeros(N1.Nx, N1.Ny);

imid1 = round(N1.Nx/2);

for ie = 1:nE
    for ia = 1:nA
        esq = esquemas7{ie};
        alg = algorithms{ia};
        fprintf('  [%-20s | %-7s] ', esq, alg);
        r = run_case(alg, u0_1, v0_1, p0_1, u_an1, D1, ...
                "CouettePoiseuille", esq, interp_pressao, numerical_param1, physical_param1);

        [errL2, errLinf] = couette_error(r.u, u_an1, xp1, yp1, N1.Nx, N1.Ny);

        S1.errL2(ie,ia)     = errL2;
        S1.errLinf(ie,ia)   = errLinf;
        S1.n_iter(ie,ia)    = r.n_iter;
        S1.t_cpu(ie,ia)     = r.t_cpu;
        S1.converged(ie,ia) = r.converged;
        S1.u_perfil{ie,ia}  = r.u(imid1, :);

        fprintf('iter=%4d  t=%.2fs  errL2=%.3e  conv=%d\n', r.n_iter, r.t_cpu, errL2, r.converged);
    end
end

y_perfil1 = yp1;
u_ref_perfil1 = u_an1(yp1);

save(fullfile(mat_dir, 'study1_precisao.mat'), 'S1', 'esquemas7', 'algorithms', ...
     'y_perfil1', 'u_ref_perfil1', 'N1');

% --- Figura 1a: perfis u(y) — um painel por esquema, 3 curvas (algoritmos) ---
fig = figure('Visible', 'off', 'Position', [50 50 1400 900]);
cores = lines(nA);
for ie = 1:nE
    subplot(3, 3, ie); hold on;
    plot(u_ref_perfil1, y_perfil1, 'k-', 'LineWidth', 2, 'DisplayName', 'Analitica');
    for ia = 1:nA
        plot(S1.u_perfil{ie,ia}, y_perfil1, '--o', 'Color', cores(ia,:), ...
            'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', char(algorithms{ia}));
    end
    xlabel('u [m/s]'); ylabel('y [m]');
    title(esquemas7{ie}, 'FontWeight', 'normal', 'Interpreter', 'none');
    legend(['Analitica', algo_labels], 'Location', 'best', 'FontSize', 7);
    grid on; hold off;
end
safe_sgtitle('Estudo 1 - Perfil u em x=L/2 - SIMPLE vs SIMPLEC vs PISO (caso base)');
print(fig, fullfile(fig_dir, 'estudo1_perfil_u.png'), '-dpng', '-r150');
close(fig);

% --- Figura 1b: barras de erro L2 ---
grouped_bar_plot(esquemas7, S1.errL2, algorithms, 'Erro L2 relativo', ...
    'Estudo 1 - Erro L2 vs solucao analitica', fullfile(fig_dir, 'estudo1_erroL2.png'));

% --- Figura 1c: barras de erro Linf ---
grouped_bar_plot(esquemas7, S1.errLinf, algorithms, 'Erro L_\infty relativo', ...
    'Estudo 1 - Erro L_\infty vs solucao analitica', fullfile(fig_dir, 'estudo1_erroLinf.png'));

% --- Figura 1d: barras de tempo de CPU ---
grouped_bar_plot(esquemas7, S1.t_cpu, algorithms, 'Tempo de CPU [s]', ...
    'Estudo 1 - Tempo de CPU', fullfile(fig_dir, 'estudo1_tempoCPU.png'));

fprintf('Estudo 1 concluido.\n');

%% ════════════════════════════════════════════════════════════════════
%  ESTUDO 3 — SUB-RELAXAÇÃO (heatmap alpha_u x alpha_p x algoritmo)
%  ════════════════════════════════════════════════════════════════════
fprintf('\n--- Estudo 3: Sub-relaxacao (heatmap alpha_u x alpha_p) ---\n');

nU = numel(alpha_u_vals); nP = numel(alpha_p_vals);
S3.iters     = nan(nU, nP, nA);
S3.t_cpu     = nan(nU, nP, nA);
S3.converged = false(nU, nP, nA);

[~, physical_param3, D3, u_an3, xp3, yp3] = build_couette_case( ...
    N3.Nx, N3.Ny, Lx0, Ly0, mu0, rho0, U0, G0, alpha_u_base, alpha_p_base, N3.tol, N3.max_iter);

u0_3 = zeros(N3.Nx-1, N3.Ny);
v0_3 = zeros(N3.Nx, N3.Ny-1);
p0_3 = zeros(N3.Nx, N3.Ny);

for ia = 1:nA
    alg = algorithms{ia};
    for iu = 1:nU
        for ip = 1:nP
            num = {alpha_u_vals(iu), alpha_p_vals(ip), N3.tol, N3.max_iter, ...
                   N3.Nx, N3.Ny, Lx0/N3.Nx, Ly0/N3.Ny, xp3, yp3};

            r = run_case(alg, u0_3, v0_3, p0_3, u_an3, D3, ...
                    "CouettePoiseuille", esquema_fixo, interp_pressao, num, physical_param3);

            S3.iters(iu,ip,ia)     = r.n_iter;
            S3.t_cpu(iu,ip,ia)     = r.t_cpu;
            S3.converged(iu,ip,ia) = r.converged;
        end
    end
    fprintf('  Algoritmo %-7s: grade %dx%d concluida\n', alg, nU, nP);
end

save(fullfile(mat_dir, 'study3_subrelaxacao.mat'), 'S3', 'algorithms', ...
     'alpha_u_vals', 'alpha_p_vals', 'N3');

% --- Figura 3a: heatmap de convergência, 1 painel por algoritmo ---
fig = figure('Visible', 'off', 'Position', [50 50 1500 450]);
for ia = 1:nA
    ax = subplot(1, nA, ia);
    iters_mat = S3.iters(:,:,ia);
    diverged  = ~S3.converged(:,:,ia);
    iters_plot = iters_mat; iters_plot(diverged) = N3.max_iter;
    conv_heatmap_plot(ax, iters_plot, diverged, alpha_u_vals, alpha_p_vals, ...
        char(algorithms{ia}), N3.max_iter);
end
safe_sgtitle('Estudo 3 - Mapa de convergencia no espaco (\alpha_u, \alpha_p)');
print(fig, fullfile(fig_dir, 'estudo3_heatmap_subrelaxacao.png'), '-dpng', '-r150');
close(fig);

% --- Figura 3b: melhor alpha_p de cada algoritmo (PISO sem relaxacao de p) ---
% Para SIMPLE e SIMPLEC: melhor alpha_p = menor n_iter (dentre convergidos)
% na linha alpha_u = alpha_u_base (ou a mais proxima disponivel).
[~, iu_base] = min(abs(alpha_u_vals - alpha_u_base));

best_ap = nan(1,2); best_iter = nan(1,3); best_time = nan(1,3); best_label = cell(1,3);
for ia = 1:2   % SIMPLE, SIMPLEC
    row_iter = squeeze(S3.iters(iu_base, :, ia));
    row_conv = squeeze(S3.converged(iu_base, :, ia));
    row_iter(~row_conv) = Inf;
    [best_iter(ia), idx_best] = min(row_iter);
    best_ap(ia) = alpha_p_vals(idx_best);
    best_time(ia) = S3.t_cpu(iu_base, idx_best, ia);
    best_label{ia} = sprintf('%s (\\alpha_p=%.1f)', char(algorithms{ia}), best_ap(ia));
end

% PISO com alpha_p = 1.0 (sem relaxacao de pressao) no mesmo alpha_u_base
num_piso = {alpha_u_vals(iu_base), 1.0, N3.tol, N3.max_iter, N3.Nx, N3.Ny, ...
            Lx0/N3.Nx, Ly0/N3.Ny, xp3, yp3};
r_piso = run_case("PISO", u0_3, v0_3, p0_3, u_an3, D3, ...
        "CouettePoiseuille", esquema_fixo, interp_pressao, num_piso, physical_param3);
best_iter(3) = r_piso.n_iter;
best_time(3) = r_piso.t_cpu;
best_label{3} = 'PISO (\alpha_p=1.0)';

fig = figure('Visible', 'off', 'Position', [100 100 900 450]);
subplot(1,2,1);
bar(best_iter); set(gca,'XTickLabel', best_label, 'XTick', 1:3);
ylabel('iteracoes ate convergir'); title('Iteracoes no melhor \alpha_p', 'FontWeight','normal');
grid on;
subplot(1,2,2);
bar(best_time); set(gca,'XTickLabel', best_label, 'XTick', 1:3);
ylabel('tempo de CPU [s]'); title('Tempo no melhor \alpha_p', 'FontWeight','normal');
grid on;
safe_sgtitle(sprintf('Estudo 3 - Melhor configuracao (\\alpha_u=%.1f)', alpha_u_vals(iu_base)));
print(fig, fullfile(fig_dir, 'estudo3_melhor_alpha_p.png'), '-dpng', '-r150');
close(fig);

save(fullfile(mat_dir, 'study3_melhor_alpha_p.mat'), 'best_ap', 'best_iter', 'best_time', 'best_label');

fprintf('Estudo 3 concluido.\n');

%% ════════════════════════════════════════════════════════════════════
%  ESTUDO 4 — SENSIBILIDADE PARAMÉTRICA (Blocos A-D x algoritmo)
%  ════════════════════════════════════════════════════════════════════
% fprintf('\n--- Estudo 4: Sensibilidade parametrica (Blocos A-D) ---\n');
% 
% blocos = struct();
% 
% % Bloco A — variacao de mu (Pe_celula alvo, dy fixo)
% dy4 = Ly0/N4.Ny;
% Pe_alvo_A = [0.1 0.2 0.5 1 2];
% blocos.A.label = 'Bloco A - variacao de \mu';
% blocos.A.xlabel = 'Pe_{celula}';
% blocos.A.x = Pe_alvo_A;
% blocos.A.mu = rho0*U0*dy4 ./ Pe_alvo_A;
% blocos.A.G  = G0*ones(size(Pe_alvo_A));
% blocos.A.U  = U0*ones(size(Pe_alvo_A));
% 
% % Bloco B — variacao de G>=0 (Pparam alvo, mu/U fixos)
% Pparam_alvo_B = [0 2.5 5 10 20 40];
% blocos.B.label = 'Bloco B - variacao de G (>=0)';
% blocos.B.xlabel = 'P_{param}';
% blocos.B.x = Pparam_alvo_B;
% blocos.B.mu = mu0*ones(size(Pparam_alvo_B));
% blocos.B.G  = Pparam_alvo_B * 2*mu0*U0 / Ly0^2;
% blocos.B.U  = U0*ones(size(Pparam_alvo_B));
% 
% % Bloco C — gradiente adverso (G<0)
% absG_alvo_C = [0.5 1 2 4 8];
% blocos.C.label = 'Bloco C - gradiente adverso (G<0)';
% blocos.C.xlabel = '|G| [Pa/m]';
% blocos.C.x = absG_alvo_C;
% blocos.C.mu = mu0*ones(size(absG_alvo_C));
% blocos.C.G  = -absG_alvo_C;
% blocos.C.U  = U0*ones(size(absG_alvo_C));
% 
% % Bloco D — variacao de U (G, mu fixos)
% U_alvo_D = [0.25 0.5 1 2 4 8];
% blocos.D.label = 'Bloco D - variacao de U';
% blocos.D.xlabel = 'U [m/s]';
% blocos.D.x = U_alvo_D;
% blocos.D.mu = mu0*ones(size(U_alvo_D));
% blocos.D.G  = G0*ones(size(U_alvo_D));
% blocos.D.U  = U_alvo_D;
% 
% nomes_bloco = {'A','B','C','D'};
% S4 = struct();
% 
% for ib = 1:numel(nomes_bloco)
%     nb = nomes_bloco{ib};
%     npts = numel(blocos.(nb).x);
%     errL2_b  = nan(npts, nA);
%     niter_b  = nan(npts, nA);
%     conv_b   = false(npts, nA);
% 
%     for ipt = 1:npts
%         mu_pt = blocos.(nb).mu(ipt);
%         G_pt  = blocos.(nb).G(ipt);
%         U_pt  = blocos.(nb).U(ipt);
% 
%         [num4, phys4, D4, u_an4, xp4, yp4] = build_couette_case( ...
%             N4.Nx, N4.Ny, Lx0, Ly0, mu_pt, rho0, U_pt, G_pt, ...
%             alpha_u_base, alpha_p_base, N4.tol, N4.max_iter);
% 
%         u0_4 = zeros(N4.Nx-1, N4.Ny);
%         v0_4 = zeros(N4.Nx, N4.Ny-1);
%         p0_4 = zeros(N4.Nx, N4.Ny);
% 
%         for ia = 1:nA
%             alg = algorithms{ia};
%             r = run_case(alg, u0_4, v0_4, p0_4, u_an4, D4, ...
%                     "CouettePoiseuille", esquema_fixo, interp_pressao, num4, phys4);
%             [errL2_b(ipt,ia), ~] = couette_error(r.u, u_an4, xp4, yp4, N4.Nx, N4.Ny);
%             niter_b(ipt,ia) = r.n_iter;
%             conv_b(ipt,ia)  = r.converged;
%         end
%     end
% 
%     S4.(nb).errL2 = errL2_b;
%     S4.(nb).n_iter = niter_b;
%     S4.(nb).converged = conv_b;
% 
%     fprintf('  %s concluido (%d pontos)\n', blocos.(nb).label, npts);
% 
%     % --- Figura: erro L2 e iteracoes, 2 subplots, 3 curvas (algoritmos) ---
%     fig = figure('Visible', 'off', 'Position', [100 100 1100 450]);
%     subplot(1,2,1); hold on;
%     for ia = 1:nA
%         plot(blocos.(nb).x, errL2_b(:,ia), '-o', 'Color', cores(ia,:), ...
%             'LineWidth', 1.5, 'DisplayName', char(algorithms{ia}));
%     end
%     xlabel(blocos.(nb).xlabel); ylabel('Erro L2 relativo');
%     title('Erro L2', 'FontWeight', 'normal');
%     legend('Location', 'best'); grid on; hold off;
% 
%     subplot(1,2,2); hold on;
%     for ia = 1:nA
%         plot(blocos.(nb).x, niter_b(:,ia), '-o', 'Color', cores(ia,:), ...
%             'LineWidth', 1.5, 'DisplayName', char(algorithms{ia}));
%     end
%     xlabel(blocos.(nb).xlabel); ylabel('Iteracoes ate convergir');
%     title('Iteracoes', 'FontWeight', 'normal');
%     legend('Location', 'best'); grid on; hold off;
% 
%     safe_sgtitle(sprintf('Estudo 4 - %s', blocos.(nb).label));
%     print(fig, fullfile(fig_dir, sprintf('estudo4_bloco%s.png', nb)), '-dpng', '-r150');
%     close(fig);
% end
% 
% save(fullfile(mat_dir, 'study4_parametrico.mat'), 'S4', 'blocos', 'algorithms', 'N4');
% 
% fprintf('Estudo 4 concluido.\n');

%% ════════════════════════════════════════════════════════════════════
%  ESTUDO 5 — ROBUSTEZ NOS 4 ESCOAMENTOS (escoamento x algoritmo)
%  ════════════════════════════════════════════════════════════════════
fprintf('\n--- Estudo 5: Robustez nos 4 escoamentos ---\n');

escoamentos5 = {"CouettePoiseuille", "LidDrivenCavity", "FlowOverCavity", "FlowOverStep"};
nF = numel(escoamentos5);

S5.n_iter    = nan(nF, nA);
S5.t_cpu     = nan(nF, nA);
S5.converged = false(nF, nA);
S5.campo     = cell(nF, nA);   % struct com u,v,p,xp,yp,Nx,Ny p/ plot

for ifl = 1:nF
    esc = escoamentos5{ifl};

    switch esc
        case "CouettePoiseuille"
            Nx5=32; Ny5=12;
            dx5=Lx0/Nx5; dy5=Ly0/Ny5;
            xp5=((1:Nx5)-0.5)'*dx5; yp5=((1:Ny5)-0.5)*dy5;
            phys5=[Lx0,Ly0,mu0,rho0,U0,G0];
            [~,u_an5]=analytical_Couette_Pouiseuille(phys5);
            D5=dif_coeff(mu0,dx5,dy5);
            num5={alpha_u_base,alpha_p_base,N5tol,N5maxit,Nx5,Ny5,dx5,dy5,xp5,yp5};
            u05=zeros(Nx5-1,Ny5); v05=zeros(Nx5,Ny5-1); p05=zeros(Nx5,Ny5);

        case "LidDrivenCavity"
            L=1; Re=100; Uc=1; rhoc=1; muc=rhoc*Uc*L/Re;
            Nx5=16; Ny5=16;
            dx5=L/Nx5; dy5=L/Ny5;
            xp5=((1:Nx5)-0.5)'*dx5; yp5=((1:Ny5)-0.5)*dy5;
            phys5=[L,L,muc,rhoc,Uc,0];
            u_an5=@(y) zeros(size(y));
            D5=dif_coeff(muc,dx5,dy5);
            num5={0.5,0.5,N5tol,N5maxit,Nx5,Ny5,dx5,dy5,xp5,yp5};
            u05=zeros(Nx5-1,Ny5); v05=zeros(Nx5,Ny5-1); p05=zeros(Nx5,Ny5);

        case "FlowOverCavity"
            L_cav=0.5; D_cav=0.5; Ly_canal=0.5;
            Lx_ante=0.5; Lx_post=1.0;
            Lx=Lx_ante+L_cav+Lx_post; Ly=Ly_canal+D_cav;
            Re=100; Uc=1; rhoc=1; muc=rhoc*Uc*D_cav/Re;
            Nx5=24; Ny5=12;
            dx5=Lx/Nx5; dy5=Ly/Ny5;
            xp5=((1:Nx5)-0.5)'*dx5; yp5=((1:Ny5)-0.5)*dy5;
            ic_s=find(xp5>Lx_ante,1,'first');
            ic_e=find(xp5>Lx_ante+L_cav,1,'first')-1;
            jc_t=find(yp5<=D_cav,1,'last');
            phys5=[Lx,Ly,muc,rhoc,Uc,0, ic_s, ic_e, jc_t];
            u_an5=@(y) 6*Uc.*((y-D_cav)./Ly_canal).*(1-(y-D_cav)./Ly_canal).*(y>D_cav);
            D5=dif_coeff(muc,dx5,dy5);
            iw=false(Nx5,Ny5);
            iw([1:ic_s-1, ic_e+1:Nx5], 1:jc_t) = true;   % sólido é o entorno da cavidade
            num5={0.5,0.5,N5tol,N5maxit,Nx5,Ny5,dx5,dy5,xp5,yp5,iw};
            u05=repmat(u_an5(yp5),Nx5-1,1); v05=zeros(Nx5,Ny5-1); p05=zeros(Nx5,Ny5);

        case "FlowOverStep"
            H=0.5; Ly_canal=1; Lx_ante=1; Lx_post=1.5;
            Lx=Lx_ante+Lx_post; Ly=Ly_canal+H;
            Re=50; Uc=1; rhoc=1; muc=rhoc*Uc*H/Re;
            Nx5=20; Ny5=10;
            dx5=Lx/Nx5; dy5=Ly/Ny5;
            xp5=((1:Nx5)-0.5)'*dx5; yp5=((1:Ny5)-0.5)*dy5;
            i_st=find(xp5>Lx_ante,1,'first');
            j_st=find(yp5<=H,1,'last');
            phys5=[Lx,Ly,muc,rhoc,Uc,0, i_st, j_st];
            u_an5=@(y) 6*Uc.*((y-H)./Ly_canal).*(1-(y-H)./Ly_canal).*(y>H);
            D5=dif_coeff(muc,dx5,dy5);
            iw=false(Nx5,Ny5);
            for ii=1:i_st-1, for jj=1:j_st, iw(ii,jj)=true; end, end
            num5={0.5,0.5,N5tol,N5maxit,Nx5,Ny5,dx5,dy5,xp5,yp5,iw};
            u05=repmat(u_an5(yp5),Nx5-1,1); v05=zeros(Nx5,Ny5-1); p05=zeros(Nx5,Ny5);
    end

    for ia = 1:nA
        alg = algorithms{ia};
        fprintf('  [%-17s | %-7s] ', esc, alg);
        r = run_case(alg, u05, v05, p05, u_an5, D5, esc, esquema_fixo, ...
                interp_pressao, num5, phys5);
        S5.n_iter(ifl,ia)    = r.n_iter;
        S5.t_cpu(ifl,ia)     = r.t_cpu;
        S5.converged(ifl,ia) = r.converged;
        S5.campo{ifl,ia}     = struct('u', r.u, 'v', r.v, 'p', r.p, ...
            'xp', xp5, 'yp', yp5, 'Nx', Nx5, 'Ny', Ny5);
        fprintf('iter=%4d  t=%.2fs  conv=%d\n', r.n_iter, r.t_cpu, r.converged);
    end
end

esc_labels5 = cellfun(@char, escoamentos5, 'UniformOutput', false);
save(fullfile(mat_dir, 'study5_robustez.mat'), 'S5', 'escoamentos5', 'esc_labels5', 'algorithms');

% --- Figura 5a: barras de iteracoes e tempo, agrupadas por escoamento ---
fig = figure('Visible', 'off', 'Position', [100 100 1300 450]);
subplot(1,2,1);
bar(S5.n_iter); set(gca,'XTickLabel', esc_labels5, 'XTick', 1:nF);
ylabel('iteracoes ate convergir'); title('Iteracoes por escoamento', 'FontWeight','normal');
legend(algo_labels, 'Location', 'best'); grid on;
subplot(1,2,2);
bar(S5.t_cpu); set(gca,'XTickLabel', esc_labels5, 'XTick', 1:nF);
ylabel('tempo de CPU [s]'); title('Tempo de CPU por escoamento', 'FontWeight','normal');
legend(algo_labels, 'Location', 'best'); grid on;
safe_sgtitle('Estudo 5 - Robustez nos 4 escoamentos');
print(fig, fullfile(fig_dir, 'estudo5_iteracoes_tempo.png'), '-dpng', '-r150');
close(fig);

% --- Figuras 5b-e: campo de magnitude de velocidade, 1 fig por escoamento,
%     3 paineis (algoritmos) ---
for ifl = 1:nF
    esc = escoamentos5{ifl};
    fig = figure('Visible', 'off', 'Position', [50 50 1500 450]);
    for ia = 1:nA
        c = S5.campo{ifl,ia};
        subplot(1, nA, ia);
        % Velocidade nos centros do VC escalar (interp simples u,v -> centro)
        umag = zeros(c.Nx, c.Ny);
        for i = 1:c.Nx
            for j = 1:c.Ny
                if i==1, ue=c.u(1,j); uw=ue; elseif i==c.Nx, uw=c.u(c.Nx-1,j); ue=uw; else, uw=c.u(i-1,j); ue=c.u(i,j); end
                if j==1, vn=c.v(i,1); vs=vn; elseif j==c.Ny, vs=c.v(i,c.Ny-1); vn=vs; else, vs=c.v(i,j-1); vn=c.v(i,j); end
                uc = 0.5*(uw+ue); vc = 0.5*(vs+vn);
                umag(i,j) = sqrt(uc^2+vc^2);
            end
        end
        imagesc(c.xp, c.yp, umag'); set(gca,'YDir','normal'); colorbar;
        xlabel('x [m]'); ylabel('y [m]');
        title(char(algorithms{ia}), 'FontWeight', 'normal');
    end
    safe_sgtitle(sprintf('Estudo 5 - |V| - %s', char(esc)));
    print(fig, fullfile(fig_dir, sprintf('estudo5_campo_%s.png', char(esc))), '-dpng', '-r150');
    close(fig);
end

fprintf('Estudo 5 concluido.\n');

%% ════════════════════════════════════════════════════════════════════
t_total = toc(t_suite);
fprintf('\n=== Suite completa em %.1f s (%.1f min) ===\n', t_total, t_total/60);
fprintf('Resultados em: %s\n', fullfile(pwd, base_dir));
