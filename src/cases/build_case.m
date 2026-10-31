function c = build_case(escoamento, varargin)
%BUILD_CASE  Monta um caso completo: geometria, malha, parâmetros e campos iniciais.
%
%   c = BUILD_CASE(escoamento) usa os valores padrão do escoamento.
%   c = BUILD_CASE(escoamento, 'nome', valor, ...) sobrescreve parâmetros.
%
%   ENTRADAS
%     escoamento : "CouettePoiseuille" | "LidDrivenCavity" | "FlowOverCavity"
%                  | "FlowOverStep"
%     'nome'     : qualquer campo padrão do escoamento (ver a tabela no
%                  README); nome desconhecido gera erro.
%
%   SAÍDA - struct c com:
%     escoamento, Nx, Ny, dx, dy, Lx, Ly, U, Re, xp, yp
%     fechado         : true se o domínio não tem entrada nem saída
%     is_wall (Nx,Ny) : máscara dos VCs escalares sólidos
%     outline (n,4)   : segmentos [x1 x2 y1 y2] do contorno sólido, para plot
%     u_analitica     : solução analítica (Couette) ou perfil de entrada
%     u0, v0, p0      : campos iniciais
%     D               : condutâncias difusivas (ver dif_coeff)
%     physical_param  : [Lx, Ly, mu, rho, U, G]
%     numerical_param : {alpha_u, alpha_p, tol, max_iter, Nx, Ny, dx, dy,
%                        xp, yp, is_wall}
%     ic_s, ic_e, jc_t : índices da cavidade (só em FlowOverCavity)
%
%   Os kernels numéricos leem physical_param e numerical_param por posição;
%   essas ordens não podem mudar.

    escoamento = string(escoamento);
    p = defaults(escoamento);

    if mod(numel(varargin), 2) ~= 0
        error('build_case: os parâmetros devem vir em pares nome/valor.');
    end
    for k = 1:2:numel(varargin)
        nome = char(varargin{k});
        if ~isfield(p, nome)
            error('build_case: parâmetro "%s" não existe em %s. Válidos: %s.', ...
                nome, escoamento, strjoin(fieldnames(p), ', '));
        end
        p.(nome) = varargin{k+1};
    end
    validateattributes(p.Nx, {'numeric'}, {'scalar', 'integer', '>=', 3}, 'build_case', 'Nx');
    validateattributes(p.Ny, {'numeric'}, {'scalar', 'integer', '>=', 3}, 'build_case', 'Ny');
    validateattributes(p.alpha_u, {'numeric'}, {'scalar', '>', 0, '<=', 1}, 'build_case', 'alpha_u');
    validateattributes(p.alpha_p, {'numeric'}, {'scalar', '>', 0, '<=', 1}, 'build_case', 'alpha_p');

    % Dimensões e viscosidade dependem do escoamento; o resto é comum
    switch escoamento
        case "CouettePoiseuille"
            Lx = p.Lx;  Ly = p.Ly;  mu = p.mu;  G = p.G;
            Re = p.rho * p.U * Ly / mu;
        case "LidDrivenCavity"
            Lx = p.L;  Ly = p.L;  Re = p.Re;  G = 0;
            mu = p.rho * p.U * p.L / Re;
        case "FlowOverCavity"
            Lx = p.Lx_ante + p.L_cav + p.Lx_post;
            Ly = p.Ly_canal + p.D_cav;
            Re = p.Re;  G = 0;
            mu = p.rho * p.U * p.Ly_canal / Re;
        case "FlowOverStep"
            Lx = p.Lx_ante + p.Lx_post;
            Ly = p.Ly_canal + p.H_step;
            Re = p.Re;  G = 0;
            mu = p.rho * p.U * p.H_step / Re;
    end

    Nx = p.Nx;  Ny = p.Ny;  U = p.U;
    dx = Lx / Nx;
    dy = Ly / Ny;
    xp = ((1:Nx) - 0.5)' * dx;   % centros em x - Nx x 1
    yp = ((1:Ny) - 0.5)  * dy;   % centros em y - 1 x Ny

    c.is_wall = false(Nx, Ny);
    c.outline = zeros(0, 4);

    switch escoamento
        case "CouettePoiseuille"
            % Perfil analítico em qualquer y; também é o perfil de entrada
            c.P_param = G * Ly^2 / (2 * mu * U);
            ua = @(y) U*(y/Ly) + (G*Ly/(2*mu)).*y.*(1 - y/Ly);

        case "LidDrivenCavity"
            ua = @(y) zeros(size(y));

        case "FlowOverCavity"
            % O sólido é o entorno da cavidade, abaixo do nível do canal
            ic_s = find(xp >= p.Lx_ante,           1, 'first');
            ic_e = find(xp <= p.Lx_ante + p.L_cav, 1, 'last');
            jc_t = find(yp <= p.D_cav,             1, 'last');
            c.is_wall([1:ic_s-1, ic_e+1:Nx], 1:jc_t) = true;
            c.ic_s = ic_s;  c.ic_e = ic_e;  c.jc_t = jc_t;

            x_s = xp(ic_s) - 0.5*dx;
            x_e = xp(ic_e) + 0.5*dx;
            y_t = yp(jc_t) + 0.5*dy;
            c.outline = [x_s x_s 0 y_t; x_e x_e 0 y_t; x_s x_e 0 0];
            ua = poiseuille_inlet(U, p.D_cav, p.Ly_canal);

        case "FlowOverStep"
            % O degrau ocupa as colunas antes de i_step, até j_step_top
            i_step     = find(xp > p.Lx_ante, 1, 'first');
            j_step_top = find(yp <= p.H_step, 1, 'last');
            c.is_wall(1:i_step-1, 1:j_step_top) = true;

            x_step = xp(i_step) - 0.5*dx;
            y_step = yp(j_step_top) + 0.5*dy;
            c.outline = [0 x_step y_step y_step; x_step x_step 0 y_step];
            ua = poiseuille_inlet(U, p.H_step, p.Ly_canal);
    end

    % Escoamentos abertos partem do perfil de entrada; os demais, do repouso,
    % porque começar da solução analítica falsearia as iterações.
    if any(escoamento == ["FlowOverCavity", "FlowOverStep"])
        c.u0 = repmat(ua(yp), Nx-1, 1);
    else
        c.u0 = zeros(Nx-1, Ny);
    end
    c.v0 = zeros(Nx, Ny-1);
    c.p0 = zeros(Nx, Ny);

    c.escoamento = escoamento;
    c.fechado = escoamento == "LidDrivenCavity";
    c.Nx = Nx;  c.Ny = Ny;  c.dx = dx;  c.dy = dy;
    c.Lx = Lx;  c.Ly = Ly;  c.U = U;    c.Re = Re;
    c.xp = xp;  c.yp = yp;
    c.u_analitica = ua;
    c.D = dif_coeff(mu, dx, dy);
    c.physical_param = [Lx, Ly, mu, p.rho, U, G];
    c.numerical_param = {p.alpha_u, p.alpha_p, p.tol, p.max_iter, ...
                         Nx, Ny, dx, dy, xp, yp, c.is_wall};
end

function p = defaults(escoamento)
%DEFAULTS  Parâmetros padrão de cada escoamento (unidades SI).
    switch escoamento
        case "CouettePoiseuille"
            p = struct('Lx', 5, 'Ly', 1, 'mu', 1e-1, 'rho', 1, 'U', 2, 'G', 4, ...
                'alpha_u', 0.8, 'alpha_p', 0.8, 'tol', 1e-8, 'max_iter', 5000, ...
                'Nx', 50, 'Ny', 20);
        case "LidDrivenCavity"
            p = struct('L', 1, 'Re', 1000, 'U', 1, 'rho', 1, ...
                'alpha_u', 0.5, 'alpha_p', 0.5, 'tol', 1e-6, 'max_iter', 5000, ...
                'Nx', 64, 'Ny', 64);
        case "FlowOverCavity"
            % Re baseado em Ly_canal e U
            p = struct('Ly_canal', 1, 'L_cav', 1, 'D_cav', 0.5, 'Lx_ante', 1, ...
                'Lx_post', 1, 'Re', 147, 'U', 1, 'rho', 1, ...
                'alpha_u', 0.5, 'alpha_p', 0.5, 'tol', 1e-6, 'max_iter', 8000, ...
                'Nx', 60, 'Ny', 30);
        case "FlowOverStep"
            % Re baseado em H_step e U
            p = struct('H_step', 0.5, 'Ly_canal', 1, 'Lx_ante', 1, 'Lx_post', 4, ...
                'Re', 100, 'U', 1, 'rho', 1, ...
                'alpha_u', 0.5, 'alpha_p', 0.5, 'tol', 1e-6, 'max_iter', 8000, ...
                'Nx', 80, 'Ny', 30);
        otherwise
            error('build_case: escoamento desconhecido "%s".', escoamento);
    end
end

function ua = poiseuille_inlet(U, y0, h)
%POISEUILLE_INLET  Perfil parabólico de velocidade média U no canal [y0, y0+h].
%   Vale 0 para y <= y0, isto é, abaixo do nível do canal.
    ua = @(y) (y > y0) .* 6 .* U .* ((y - y0)./h) .* (1 - (y - y0)./h);
end
