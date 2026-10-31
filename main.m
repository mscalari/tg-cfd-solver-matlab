clc; 
clear; 
close all;

%% =========================================================
%  DEFINIÇÕES 
% =========================================================

escoamento     = "CouettePoiseuille"; % "CouettePoiseuille" | "LidDrivenCavity" | "FlowOverCavity" | "FlowOverStep"
interp_pressao = "Linear";
esquema_convec = 'Upwind 1a ordem'; % "Upwind 1a ordem" | "Diferenças Centrais" | "Power Law" | "Upwind 2a ordem" | "QUICK" | "MUSCL" | "HRIC"
algorithm      = "PISO"; % "SIMPLE" | "SIMPLEC" | "PISO"

%% =========================================================
%  PARÂMETROS FÍSICOS
% =========================================================

switch escoamento

    case "CouettePoiseuille"
        Lx  = 5;      % comprimento do canal [m]
        Ly  = 1;      % altura do canal [m]
        mu  = 1e-1;   % viscosidade dinâmica [Pa·s]
        rho = 1.0;    % massa específica [kg/m³]
        U   = 2;      % velocidade da placa superior [m/s]
        G   = 4;      % gradiente de pressão imposto (-dp/dx) [Pa/m]

    case "LidDrivenCavity"
        L   = 1.0;    % lado da cavidade [m]
        Re  = 1000;   % número de Reynolds
        U   = 1.0;    % velocidade da tampa [m/s]
        rho = 1.0;    % massa específica [kg/m³]
        mu  = rho * U * L / Re;
        Lx  = L;
        Ly  = L;
        G   = 0;

    case "FlowOverCavity"
        % Geometria total do domínio
        Ly_canal = 1.0;   % altura da região do canal (acima da cavidade) [m]
        L_cav    = 1.0;   % comprimento da cavidade [m]
        D_cav    = 0.5;   % profundidade da cavidade [m]
        Lx_ante  = 1.0;   % comprimento do canal antes da cavidade [m]
        Lx_post  = 1.0;   % comprimento do canal depois da cavidade [m]

        % Domínio total
        Lx = Lx_ante + L_cav + Lx_post;  % comprimento total [m]
        Ly = Ly_canal + D_cav;            % altura total (canal + cavidade) [m]

        % Propriedades do fluido
        rho = 1.0;    % massa específica [kg/m³]
        U   = 1.0;    % velocidade média de entrada [m/s]
        Re  = 147;    % Reynolds baseado em Ly_canal e U
        mu  = rho * U * Ly_canal / Re;
        G   = 0;      % sem gradiente imposto - escoamento dirigido pela entrada

        fprintf('\n Flow Over Cavity\n');
        fprintf(' Re = %.0f  |  L/D = %.2f  |  mu = %.5f\n\n', ...
            Re, L_cav/D_cav, mu);

    case "FlowOverStep"
        % Escoamento sobre degrau backward-facing.
        % O degrau ocupa x em [0, Lx_ante] e y em [0, H_step].
        % O canal acima do degrau tem altura Ly_canal = Ly - H_step.
        % A partir de x = Lx_ante o escoamento se expande para Ly inteiro.
        H_step   = 0.5;   % altura do degrau [m]
        Ly_canal = 1.0;   % altura do canal acima do degrau [m]
        Lx_ante  = 1.0;   % comprimento do trecho com o degrau (antes da queda) [m]
        Lx_post  = 4.0;   % comprimento após o degrau (região de reattachment) [m]

        % Domínio total
        Lx = Lx_ante + Lx_post;   % comprimento total [m]
        Ly = Ly_canal + H_step;    % altura total [m]

        % Propriedades do fluido
        rho = 1.0;    % massa específica [kg/m³]
        U   = 1.0;    % velocidade média de entrada no canal [m/s]
        Re  = 100;    % Reynolds baseado em H_step e U
        mu  = rho * U * H_step / Re;
        G   = 0;      % sem gradiente imposto

        fprintf('\n Flow Over Step\n');
        fprintf(' Re = %.0f  |  H_step = %.2f  |  Ly_canal = %.2f  |  mu = %.5f\n\n', ...
            Re, H_step, Ly_canal, mu);

end

%% =========================================================
%  physical_param — vetor de parâmetros físicos
%
%  Posições fixas para todos os escoamentos:
%    (1) Lx   (2) Ly   (3) mu   (4) rho   (5) U   (6) G
%
%  FlowOverCavity adiciona índices da cavidade:
%    (7) i_cav_start   (8) i_cav_end   (9) j_cav_top
%  Esses índices são calculados após a malha ser definida.
% =========================================================

physical_param = [Lx, Ly, mu, rho, U, G];   % é estendido abaixo para FlowOverCavity

%% =========================================================
%  PARÂMETROS NUMÉRICOS
% =========================================================

switch escoamento

    case "CouettePoiseuille"
        alpha_u  = 0.8;
        alpha_p  = 0.8;
        tol      = 1e-8;
        max_iter = 5000;
        Nx = 50;
        Ny = 20;

    case "LidDrivenCavity"
        alpha_u  = 0.5;
        alpha_p  = 0.5;
        tol      = 1e-6;
        max_iter = 5000;
        Nx = 64;
        Ny = 64;

    case "FlowOverCavity"
        alpha_u  = 0.5;
        alpha_p  = 0.5;
        tol      = 1e-6;
        max_iter = 8000;
        Nx = 60;   % células em x (domínio total Lx)
        Ny = 30;   % células em y (domínio total Ly)

    case "FlowOverStep"
        alpha_u  = 0.5;
        alpha_p  = 0.5;
        tol      = 1e-6;
        max_iter = 8000;
        Nx = 80;   % células em x (domínio total Lx)
        Ny = 30;   % células em y (domínio total Ly)

end

dx = Lx / Nx;
dy = Ly / Ny;
xp = ((1:Nx) - 0.5)' * dx;   % centros em x — Nx x 1
yp = ((1:Ny) - 0.5)  * dy;   % centros em y — 1  x Ny

numerical_param = {alpha_u, alpha_p, tol, max_iter, Nx, Ny, dx, dy, xp, yp};

%% =========================================================
%  ÍNDICES DA CAVIDADE (só para FlowOverCavity)
%
%  A malha cobre o domínio total [0, Lx] x [0, Ly].
%  A cavidade ocupa x em [Lx_ante, Lx_ante+L_cav] e y em [0, D_cav].
%
%  i_cav_start: primeiro índice i cujo centro xp(i) está dentro da cavidade
%  i_cav_end  : último  índice i cujo centro xp(i) está dentro da cavidade
%  j_cav_top  : último  índice j cujo centro yp(j) está dentro da cavidade
%               (a interface canal-cavidade fica entre j_cav_top e j_cav_top+1)
% =========================================================

if escoamento == "FlowOverCavity"
    i_cav_start = find(xp >= Lx_ante,          1, 'first');
    i_cav_end   = find(xp <= Lx_ante + L_cav,  1, 'last');
    j_cav_top   = find(yp <= D_cav,            1, 'last');

    % Estender physical_param com os indices da cavidade
    physical_param = [physical_param, i_cav_start, i_cav_end, j_cav_top];

    % Mascara de paredes solidas - VCs que sao solidos (fundo do canal fora da cavidade)
    % is_wall(i,j) = true  se o VC escalar (i,j) for parede solida
    % Regra: solido quando j <= j_cav_top E (i < i_cav_start OU i > i_cav_end)
    is_wall = false(Nx, Ny);
    for ii = 1:Nx
        for jj = 1:j_cav_top
            if ii < i_cav_start || ii > i_cav_end
                is_wall(ii,jj) = true;
            end
        end
    end

    % Adicionar is_wall como 11o elemento de numerical_param
    numerical_param{11} = is_wall;

    fprintf(' Malha: %dx%d  |  dx=%.4f  dy=%.4f\n', Nx, Ny, dx, dy);
    fprintf(' Cavidade: i=%d..%d  j=1..%d  (%d x %d celulas)\n', ...
        i_cav_start, i_cav_end, j_cav_top, ...
        i_cav_end-i_cav_start+1, j_cav_top);
    fprintf(' VCs solidos: %d de %d\n\n', sum(is_wall(:)), Nx*Ny);
end

if escoamento == "FlowOverStep"
    % i_step: primeira coluna de VC escalar após o degrau (x > Lx_ante)
    % j_step_top: último nível de VC dentro da região do degrau (y <= H_step)
    i_step     = find(xp > Lx_ante, 1, 'first');
    j_step_top = find(yp <= H_step, 1, 'last');

    % Estender physical_param com os índices do degrau
    % Posição 7: i_step   Posição 8: j_step_top
    physical_param = [physical_param, i_step, j_step_top];

    % Máscara de paredes sólidas — o degrau ocupa i < i_step e j <= j_step_top
    % is_wall(i,j) = true quando o VC está dentro do bloco do degrau
    is_wall = false(Nx, Ny);
    for ii = 1:i_step-1
        for jj = 1:j_step_top
            is_wall(ii,jj) = true;
        end
    end

    % Adicionar is_wall como 11o elemento de numerical_param
    numerical_param{11} = is_wall;

    fprintf(' Malha: %dx%d  |  dx=%.4f  dy=%.4f\n', Nx, Ny, dx, dy);
    fprintf(' Degrau: i=1..%d  j=1..%d  (%d x %d celulas solidas)\n', ...
        i_step-1, j_step_top, i_step-1, j_step_top);
    fprintf(' VCs solidos: %d de %d\n\n', sum(is_wall(:)), Nx*Ny);
end

%% =========================================================
%  SOLUÇÃO ANALÍTICA / REFERÊNCIA
% =========================================================

switch escoamento

    case "CouettePoiseuille"
        [P_param, u_analitica] = analytical_Couette_Pouiseuille(physical_param);
        fprintf('\n Parâmetro de Couette-Poiseuille P = %.4f\n', P_param);
        fprintf(' Re = %.1f\n\n', rho*U*Ly/mu);

    case "LidDrivenCavity"
        u_analitica = @(y) zeros(size(y));
        fprintf('\n Re = %.0f  |  mu = %.5f  |  U_lid = %.2f m/s\n\n', Re, mu, U);

    case "FlowOverCavity"
        % Perfil de entrada: Poiseuille no canal (altura Ly_canal),
        % com velocidade média U.
        % u(y) = 6*U * (y - D_cav)/Ly_canal * (1 - (y - D_cav)/Ly_canal)
        % onde (y - D_cav) é a coordenada local dentro do canal.
        % Para y <= D_cav (dentro da cavidade), u_entrada = 0.
        u_analitica = @(y) (y > D_cav) .* 6 .* U .* ((y - D_cav)./Ly_canal) .* (1 - (y - D_cav)./Ly_canal);

    case "FlowOverStep"
        % Perfil de entrada: Poiseuille no canal acima do degrau (altura Ly_canal).
        % u(y) = 6*U * (y - H_step)/Ly_canal * (1 - (y - H_step)/Ly_canal)
        % Para y <= H_step (dentro do degrau), u_analitica retorna 0.
        u_analitica = @(y) (y > H_step) .* 6 .* U .* ((y - H_step)./Ly_canal) .* (1 - (y - H_step)./Ly_canal);

end

%% =========================================================
%  COEFICIENTES DE DIFUSÃO
% =========================================================

D = dif_coeff(mu, dx, dy);

%% =========================================================
%  INICIALIZAÇÃO DOS CAMPOS
% =========================================================

switch escoamento

    case "CouettePoiseuille"
        % Partida do repouso: começar da solução analítica falsearia as iterações
        p = zeros(Nx, Ny);
        u = zeros(Nx-1, Ny);
        v = zeros(Nx, Ny-1);

    case "LidDrivenCavity"
        p = zeros(Nx, Ny);
        u = zeros(Nx-1, Ny);
        v = zeros(Nx, Ny-1);

    case "FlowOverCavity"
        % Chute inicial: perfil de entrada em toda a extensão do canal.
        % Dentro da cavidade (j <= j_cav_top), u_analitica retorna 0.
        u = repmat(u_analitica(yp), Nx-1, 1);   % (Nx-1) x Ny
        v = zeros(Nx, Ny-1);
        p = zeros(Nx, Ny);

    case "FlowOverStep"
        % Chute inicial: perfil de entrada em toda a extensão do domínio.
        % Dentro do degrau (i < i_step e j <= j_step_top), u_analitica
        % retorna 0 pois y <= H_step. Os VCs sólidos serão zerados pelo
        % is_wall no momentum_predictor.
        u = repmat(u_analitica(yp), Nx-1, 1);   % (Nx-1) x Ny
        v = zeros(Nx, Ny-1);
        p = zeros(Nx, Ny);

end

%% =========================================================
%  SOLVER SIMPLE
% =========================================================

fprintf(' Iter      Residuo de massa\n');
fprintf('--------------------------------\n');

switch algorithm
    case "SIMPLE"
        [u, v, p, res_hist] = simple_solver( ...
            u, v, p, ...
            u_analitica, ...
            D, ...
            escoamento, esquema_convec, interp_pressao, ...
            numerical_param, physical_param);

    case "SIMPLEC"
        [u, v, p, res_hist] = simplec_solver( ...
            u, v, p, ...
            u_analitica, ...
            D, ...
            escoamento, esquema_convec, interp_pressao, ...
            numerical_param, physical_param);

    case "PISO"
        [u, v, p, res_hist] = piso_solver( ...
            u, v, p, ...
            u_analitica, ...
            D, ...
            escoamento, esquema_convec, interp_pressao, ...
            numerical_param, physical_param);
end

fprintf('\n Iteracoes: %d  |  Residuo final: %.4e\n\n', length(res_hist), res_hist(end));

%% =========================================================
%  PÓS-PROCESSAMENTO
% =========================================================

switch escoamento

    % ──────────────────────────────────────────────────────
    case "CouettePoiseuille"
    % ──────────────────────────────────────────────────────
        imid   = round(Nx/2);
        y_plot = [0, yp, Ly];
        u_num  = [0; u(imid,:)'; U];
        u_ref  = [0; u_analitica(yp)'; U];

        u_centro = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
                if i == 1
                    ue = u(1,j); uw = u_analitica(yp(j));
                elseif i == Nx
                    ue = u(Nx-1,j); uw = u(Nx-1,j);
                else
                    ue = u(i,j); uw = u(i-1,j);
                end
                u_centro(i,j) = 0.5*(uw + ue);
            end
        end
        u_ref_grid = repmat(u_analitica(yp), Nx, 1);
        erro_L2    = norm(u_centro - u_ref_grid,'fro') / norm(u_ref_grid,'fro');
        fprintf(' Erro L2 relativo: %.4e\n\n', erro_L2);

        figure('Name','Perfil u - Couette-Poiseuille')
        hold on; grid on; box on;
        plot(u_ref, y_plot, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analitica')
        plot(u_num, y_plot, 'b-',  'LineWidth', 2,   'DisplayName', esquema_convec)
        xlabel('u  [m/s]', 'FontSize', 12)
        ylabel('y  [m]',   'FontSize', 12)
        title(sprintf('Perfil u - %s\nErro L2 = %.3e', esquema_convec, erro_L2), ...
            'FontSize', 11, 'FontWeight', 'normal')
        legend('Location', 'northwest', 'FontSize', 10)

        figure('Name','Convergencia')
        semilogy(res_hist, 'b-', 'LineWidth', 1.5)
        xlabel('Iteracao',                             'FontSize', 12)
        ylabel('Residuo de massa (norma Frobenius)',   'FontSize', 12)
        title(sprintf('Convergencia %s - Couette-Poiseuille', algorithm), 'FontSize', 11, 'FontWeight', 'normal')
        grid on; box on;

        v_centro = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
                if j == 1
                    vs = 0; vn = v(i,1);
                elseif j == Ny
                    vs = v(i,Ny-1); vn = 0;
                else
                    vs = v(i,j-1); vn = v(i,j);
                end
                v_centro(i,j) = 0.5*(vs + vn);
            end
        end

        [X, Y] = meshgrid(xp, yp);

        figure('Name','Campo u - Couette-Poiseuille')
        contourf(xp, yp, u_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); axis tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de velocidade u  [m/s] - %s', esquema_convec), ...
            'FontSize', 11, 'FontWeight', 'normal')

        figure('Name','Campo v - Couette-Poiseuille')
        contourf(xp, yp, v_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); axis tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de velocidade v  [m/s] - %s', esquema_convec), ...
            'FontSize', 11, 'FontWeight', 'normal')

        figure('Name','Vetores velocidade - Couette-Poiseuille')
        speed = sqrt(u_centro'.^2 + v_centro'.^2);
        contourf(xp, yp, speed, 20, 'LineStyle', 'none'); colorbar; colormap('jet'); hold on;
        skip_x = max(1, round(Nx/20));
        skip_y = max(1, round(Ny/10));
        ix = 1:skip_x:Nx;  iy = 1:skip_y:Ny;
        quiver(X(iy,ix), Y(iy,ix), u_centro(ix,iy)', v_centro(ix,iy)', 0.6, 'w', 'LineWidth', 1)
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Magnitude da velocidade [m/s] + vetores - %s', esquema_convec), ...
            'FontSize', 11, 'FontWeight', 'normal')
        axis tight

    % ──────────────────────────────────────────────────────
    case "LidDrivenCavity"
    % ──────────────────────────────────────────────────────
        u_centro = zeros(Nx, Ny);
        v_centro = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
                if i == 1
                    uw = 0; ue = u(1,j);
                elseif i == Nx
                    uw = u(Nx-1,j); ue = 0;
                else
                    uw = u(i-1,j); ue = u(i,j);
                end
                u_centro(i,j) = 0.5*(uw + ue);

                if j == 1
                    vs = 0; vn = v(i,1);
                elseif j == Ny
                    vs = v(i,Ny-1); vn = 0;
                else
                    vs = v(i,j-1); vn = v(i,j);
                end
                v_centro(i,j) = 0.5*(vs + vn);
            end
        end

        imid = round(Nx/2);
        jmid = round(Ny/2);

        y_vert  = [0, yp, Ly];
        u_vert  = [0; u_centro(imid,:)'; U];
        x_horiz = [0, xp', Lx];
        v_horiz = [0, v_centro(:,jmid)', 0];

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

        u_ref_int  = interp1(ace_y, ace_u, yp/Ly, 'linear', 'extrap');
        u_mid_norm = u_centro(imid,:) / U;
        erro_ref   = norm(u_mid_norm - u_ref_int) / norm(u_ref_int);
        fprintf(' Erro L2 vs Acenumerics benchmark (perfil u): %.4e\n\n', erro_ref);

        figure('Name','Convergencia')
        semilogy(res_hist, 'b-', 'LineWidth', 1.5)
        xlabel('Iteracao', 'FontSize', 12)
        ylabel('Residuo de massa (norma Frobenius)', 'FontSize', 12)
        title(sprintf('Convergencia %s - Cavity Re=%d, %s', algorithm, Re, esquema_convec), ...
            'FontSize', 11, 'FontWeight', 'normal')
        grid on; box on;

        figure('Name','Perfil u')
        hold on; grid on; box on;
        plot(u_vert/U, y_vert/Ly, 'b-', 'LineWidth', 2, 'DisplayName', esquema_convec)
        plot(ace_u, ace_y, 'ko', 'MarkerSize', 7, 'DisplayName', 'Acenumerics benchmark')
        xlabel('u / U_{lid}', 'FontSize', 12); ylabel('y / L', 'FontSize', 12)
        title(sprintf('Perfil u em x=L/2 - Re=%d\nErro L2 = %.3e', Re, erro_ref), ...
            'FontSize', 11, 'FontWeight', 'normal')
        legend('Location', 'northwest', 'FontSize', 10)

        figure('Name','Perfil v')
        hold on; grid on; box on;
        plot(x_horiz/Lx, v_horiz/U, 'r-', 'LineWidth', 2, 'DisplayName', esquema_convec)
        plot(ace_x, ace_v, 'ko', 'MarkerSize', 7, 'DisplayName', 'Acenumerics benchmark')
        xlabel('x / L', 'FontSize', 12); ylabel('v / U_{lid}', 'FontSize', 12)
        title(sprintf('Perfil v em y=L/2 - Re=%d', Re), 'FontSize', 11, 'FontWeight', 'normal')
        legend('Location', 'northeast', 'FontSize', 10)

        figure('Name','Streamlines')
        [X, Y] = meshgrid(xp, yp);
        streamslice(X, Y, u_centro', v_centro', 2)
        axis equal tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Linhas de corrente - Re=%d, %s, %dx%d', Re, esquema_convec, Nx, Ny), ...
            'FontSize', 11, 'FontWeight', 'normal')

        figure('Name','Pressao')
        contourf(xp, yp, p', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); axis equal tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de pressao - Re=%d', Re), 'FontSize', 11, 'FontWeight', 'normal')

        [X, Y] = meshgrid(xp, yp);

        figure('Name','Campo u - Cavity')
        contourf(xp, yp, u_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); axis equal tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de velocidade u  [m/s] - Re=%d', Re), ...
            'FontSize', 11, 'FontWeight', 'normal')

        figure('Name','Campo v - Cavity')
        contourf(xp, yp, v_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); axis equal tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de velocidade v  [m/s] - Re=%d', Re), ...
            'FontSize', 11, 'FontWeight', 'normal')

        figure('Name','Vetores velocidade - Cavity')
        speed = sqrt(u_centro'.^2 + v_centro'.^2);
        contourf(xp, yp, speed, 20, 'LineStyle', 'none'); colorbar; colormap('jet');
        hold on; axis equal tight
        skip = max(1, round(Nx/16));
        ix = 1:skip:Nx;  iy = 1:skip:Ny;
        quiver(X(iy,ix), Y(iy,ix), u_centro(ix,iy)', v_centro(ix,iy)', 0.6, 'w', 'LineWidth', 1)
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Magnitude da velocidade [m/s] + vetores - Re=%d', Re), ...
            'FontSize', 11, 'FontWeight', 'normal')

    % ──────────────────────────────────────────────────────
    case "FlowOverCavity"
    % ──────────────────────────────────────────────────────

        % Recuperar índices da cavidade
        ic_s = physical_param(7);
        ic_e = physical_param(8);
        jc_t = physical_param(9);

        % Recuperar mascara de paredes
        is_wall = numerical_param{11};

        % Interpolar u e v para centros dos VCs escalares
        % VCs solidos recebem NaN para nao aparecer nos plots
        % Para VCs fluidos adjacentes a paredes, as faces que tocam
        % solidos contribuem com zero (no-slip).
        u_centro = zeros(Nx, Ny);
        v_centro = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
                if is_wall(i,j)
                    u_centro(i,j) = NaN;
                    v_centro(i,j) = NaN;
                    continue
                end

                % ── u_centro ──────────────────────────────────────────────
                % Face oeste: parede W externa ou VC vizinho sólido → uw = 0
                if i == 1
                    uw = u_analitica(yp(j));
                elseif is_wall(i-1,j)
                    uw = 0;
                else
                    uw = u(i-1,j);
                end
                % Face leste: saída livre, parede E externa ou VC vizinho sólido → ue = 0
                if i == Nx
                    ue = u(Nx-1,j);
                elseif is_wall(i+1,j)
                    ue = 0;
                else
                    ue = u(i,j);
                end
                u_centro(i,j) = 0.5*(uw + ue);

                % ── v_centro ──────────────────────────────────────────────
                % Face sul: fundo externo ou VC vizinho sólido → vs = 0
                if j == 1
                    vs = 0;
                elseif is_wall(i,j-1)
                    vs = 0;
                else
                    vs = v(i,j-1);
                end
                % Face norte: teto externo ou VC vizinho sólido → vn = 0
                if j == Ny
                    vn = 0;
                elseif is_wall(i,j+1)
                    vn = 0;
                else
                    vn = v(i,j);
                end
                v_centro(i,j) = 0.5*(vs + vn);
            end
        end

        % Mascarar pressão nos VCs sólidos para o plot
        p_plot = p;
        p_plot(is_wall) = NaN;

        [X, Y] = meshgrid(xp, yp);

        % Figura 1 — Convergencia
        figure('Name','Convergencia - FlowOverCavity')
        semilogy(res_hist, 'b-', 'LineWidth', 1.5)
        xlabel('Iteracao', 'FontSize', 12)
        ylabel('Residuo de massa (norma Frobenius)', 'FontSize', 12)
        title(sprintf('Convergencia %s - Flow Over Cavity Re=%d, %s', algorithm, Re, esquema_convec), ...
            'FontSize', 11, 'FontWeight', 'normal')
        grid on; box on;

        % Figura 2 — Streamlines (principal figura de analise)
        figure('Name','Streamlines - FlowOverCavity')
        streamslice(X, Y, u_centro', v_centro', 2)
        hold on
        % Contorno da cavidade (linhas das paredes internas)
        x_cav_s = xp(ic_s) - 0.5*dx;
        x_cav_e = xp(ic_e) + 0.5*dx;
        y_cav_t = yp(jc_t) + 0.5*dy;
        plot([x_cav_s x_cav_s], [0 y_cav_t], 'k-', 'LineWidth', 2)  % parede esquerda
        plot([x_cav_e x_cav_e], [0 y_cav_t], 'k-', 'LineWidth', 2)  % parede direita
        plot([x_cav_s x_cav_e], [0 0],        'k-', 'LineWidth', 2)  % fundo
        axis tight
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Linhas de corrente - Re=%d, L/D=%.1f, %s', ...
            Re, L_cav/D_cav, esquema_convec), 'FontSize', 11, 'FontWeight', 'normal')

        % Figura 3 — Campo u
        figure('Name','Campo u - FlowOverCavity')
        contourf(xp, yp, u_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); hold on; axis tight
        plot([x_cav_s x_cav_s], [0 y_cav_t], 'w-', 'LineWidth', 1.5)
        plot([x_cav_e x_cav_e], [0 y_cav_t], 'w-', 'LineWidth', 1.5)
        plot([x_cav_s x_cav_e], [0 0],        'w-', 'LineWidth', 1.5)
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de velocidade u  [m/s] - Re=%d', Re), ...
            'FontSize', 11, 'FontWeight', 'normal')

        % Figura 4 — Campo v
        figure('Name','Campo v - FlowOverCavity')
        contourf(xp, yp, v_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); hold on; axis tight
        plot([x_cav_s x_cav_s], [0 y_cav_t], 'w-', 'LineWidth', 1.5)
        plot([x_cav_e x_cav_e], [0 y_cav_t], 'w-', 'LineWidth', 1.5)
        plot([x_cav_s x_cav_e], [0 0],        'w-', 'LineWidth', 1.5)
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de velocidade v  [m/s] - Re=%d', Re), ...
            'FontSize', 11, 'FontWeight', 'normal')

        % Figura 5 — Vetores + magnitude
        figure('Name','Vetores velocidade - FlowOverCavity')
        speed = sqrt(u_centro'.^2 + v_centro'.^2);
        contourf(xp, yp, speed, 20, 'LineStyle', 'none'); colorbar; colormap('jet');
        hold on; axis tight
        skip_x = max(1, round(Nx/24));
        skip_y = max(1, round(Ny/12));
        ix = 1:skip_x:Nx;  iy = 1:skip_y:Ny;
        quiver(X(iy,ix), Y(iy,ix), u_centro(ix,iy)', v_centro(ix,iy)', 0.6, 'w', 'LineWidth', 1)
        plot([x_cav_s x_cav_s], [0 y_cav_t], 'k-', 'LineWidth', 2)
        plot([x_cav_e x_cav_e], [0 y_cav_t], 'k-', 'LineWidth', 2)
        plot([x_cav_s x_cav_e], [0 0],        'k-', 'LineWidth', 2)
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Magnitude da velocidade [m/s] + vetores - Re=%d', Re), ...
            'FontSize', 11, 'FontWeight', 'normal')

        % Figura 6 — Campo de pressao
        figure('Name','Pressao - FlowOverCavity')
        contourf(xp, yp, p_plot', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); hold on; axis tight
        plot([x_cav_s x_cav_s], [0 y_cav_t], 'w-', 'LineWidth', 1.5)
        plot([x_cav_e x_cav_e], [0 y_cav_t], 'w-', 'LineWidth', 1.5)
        plot([x_cav_s x_cav_e], [0 0],        'w-', 'LineWidth', 1.5)
        xlabel('x  [m]', 'FontSize', 12); ylabel('y  [m]', 'FontSize', 12)
        title(sprintf('Campo de pressao - Re=%d', Re), 'FontSize', 11, 'FontWeight', 'normal')


    % ──────────────────────────────────────────────────────
    case "FlowOverStep"
    % ──────────────────────────────────────────────────────

        % Recuperar índices do degrau
        i_st = physical_param(7);
        j_st = physical_param(8);

        % Recuperar máscara de paredes
        is_wall = numerical_param{11};

        % Interpolar u e v para centros dos VCs escalares
        % VCs sólidos recebem NaN; VCs fluidos adjacentes à parede → face = 0
        u_centro = zeros(Nx, Ny);
        v_centro = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
                if is_wall(i,j)
                    u_centro(i,j) = NaN;
                    v_centro(i,j) = NaN;
                    continue
                end

                % ── u_centro ──────────────────────────────────────────────
                if i == 1
                    uw = u_analitica(yp(j));
                elseif is_wall(i-1,j)
                    uw = 0;
                else
                    uw = u(i-1,j);
                end
                if i == Nx
                    ue = u(Nx-1,j);
                elseif is_wall(i+1,j)
                    ue = 0;
                else
                    ue = u(i,j);
                end
                u_centro(i,j) = 0.5*(uw + ue);

                % ── v_centro ──────────────────────────────────────────────
                if j == 1
                    vs = 0;
                elseif is_wall(i,j-1)
                    vs = 0;
                else
                    vs = v(i,j-1);
                end
                if j == Ny
                    vn = 0;
                elseif is_wall(i,j+1)
                    vn = 0;
                else
                    vn = v(i,j);
                end
                v_centro(i,j) = 0.5*(vs + vn);
            end
        end

        % Mascarar pressão nos VCs sólidos
        p_plot = p;
        p_plot(is_wall) = NaN;

        [X, Y] = meshgrid(xp, yp);

        % Coordenadas da parede do degrau para o plot
        x_step = xp(i_st) - 0.5*dx;   % posição x da parede vertical
        y_step = yp(j_st) + 0.5*dy;   % altura do degrau

        figure('Name','FlowOverStep')

        % ── Subplot 1: Convergência ───────────────────────────────────────
        subplot(2,3,1)
        semilogy(res_hist, 'b-', 'LineWidth', 1.5)
        xlabel('Iteracao', 'FontSize', 10)
        ylabel('Residuo de massa', 'FontSize', 10)
        title(sprintf('Convergencia %s', algorithm), 'FontSize', 9, 'FontWeight', 'normal')
        grid on; box on;

        % ── Subplot 2: Streamlines ────────────────────────────────────────
        subplot(2,3,2)
        streamslice(X, Y, u_centro', v_centro', 2)
        hold on
        % Contorno do degrau: parede vertical e fundo até o degrau
        plot([0 x_step],     [y_step y_step], 'k-', 'LineWidth', 2)  % topo do degrau (horizontal)
        plot([x_step x_step],[0 y_step],       'k-', 'LineWidth', 2)  % parede vertical
        axis tight
        xlabel('x  [m]', 'FontSize', 10); ylabel('y  [m]', 'FontSize', 10)
        title(sprintf('Linhas de corrente - Re=%d, %s', Re, esquema_convec), ...
            'FontSize', 9, 'FontWeight', 'normal')

        % ── Subplot 3: Campo u ────────────────────────────────────────────
        subplot(2,3,3)
        contourf(xp, yp, u_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); hold on; axis tight
        plot([0 x_step],     [y_step y_step], 'w-', 'LineWidth', 1.5)
        plot([x_step x_step],[0 y_step],       'w-', 'LineWidth', 1.5)
        xlabel('x  [m]', 'FontSize', 10); ylabel('y  [m]', 'FontSize', 10)
        title(sprintf('Campo u  [m/s] - Re=%d', Re), 'FontSize', 9, 'FontWeight', 'normal')

        % ── Subplot 4: Campo v ────────────────────────────────────────────
        subplot(2,3,4)
        contourf(xp, yp, v_centro', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); hold on; axis tight
        plot([0 x_step],     [y_step y_step], 'w-', 'LineWidth', 1.5)
        plot([x_step x_step],[0 y_step],       'w-', 'LineWidth', 1.5)
        xlabel('x  [m]', 'FontSize', 10); ylabel('y  [m]', 'FontSize', 10)
        title(sprintf('Campo v  [m/s] - Re=%d', Re), 'FontSize', 9, 'FontWeight', 'normal')

        % ── Subplot 5: Magnitude + vetores ───────────────────────────────
        subplot(2,3,5)
        speed = sqrt(u_centro'.^2 + v_centro'.^2);
        contourf(xp, yp, speed, 20, 'LineStyle', 'none'); colorbar; colormap('jet');
        hold on; axis tight
        skip_x = max(1, round(Nx/24));
        skip_y = max(1, round(Ny/12));
        ix = 1:skip_x:Nx;  iy = 1:skip_y:Ny;
        quiver(X(iy,ix), Y(iy,ix), u_centro(ix,iy)', v_centro(ix,iy)', 0.6, 'w', 'LineWidth', 1)
        plot([0 x_step],     [y_step y_step], 'k-', 'LineWidth', 2)
        plot([x_step x_step],[0 y_step],       'k-', 'LineWidth', 2)
        xlabel('x  [m]', 'FontSize', 10); ylabel('y  [m]', 'FontSize', 10)
        title('Magnitude [m/s] + vetores', 'FontSize', 9, 'FontWeight', 'normal')

        % ── Subplot 6: Pressão ────────────────────────────────────────────
        subplot(2,3,6)
        contourf(xp, yp, p_plot', 20, 'LineStyle', 'none')
        colorbar; colormap('jet'); hold on; axis tight
        plot([0 x_step],     [y_step y_step], 'w-', 'LineWidth', 1.5)
        plot([x_step x_step],[0 y_step],       'w-', 'LineWidth', 1.5)
        xlabel('x  [m]', 'FontSize', 10); ylabel('y  [m]', 'FontSize', 10)
        title(sprintf('Pressao - Re=%d', Re), 'FontSize', 9, 'FontWeight', 'normal')

        sgtitle(sprintf('Flow Over Step — Re=%d, %s', Re, esquema_convec), 'FontSize', 12)

end
