function [u_star, v_star, aPU, aU, aPV, aV, is_wall, is_wall_u, is_wall_v] = ...
    momentum_predictor(u, v, p, u_analitica, D, escoamento, esquema_convec, ...
                        interp_pressao, numerical_param, physical_param)
%MOMENTUM_PREDICTOR  Passo preditor comum ao SIMPLE, SIMPLEC e PISO.
%
%   Resolve as equacoes de quantidade de movimento u*, v* a partir dos
%   campos u, v, p da iteracao anterior e devolve, por no de velocidade,
%   o coeficiente central aP e o vetor de coeficientes vizinhos
%   [aE,aW,aN,aS], necessarios para:
%
%     - SIMPLEC: recalcular
%           d = A / (aP/alpha_u - sum(anb))
%     - PISO   : recalcular a "correcao de vizinhanca"
%           sum(anb .* (u**_nb - u*_nb)) / (aP/alpha_u)
%       do segundo corretor.
%
%   O passo preditor combina mass_flux, espac_disc, boundary_conditions,
%   pressure_interp e sub-relaxacao.
%
%   ENTRADAS
%     u, v, p        : campos da iteracao anterior
%     u_analitica    : perfil de entrada
%     D              : condutancias de dif_coeff
%     escoamento     : nome do escoamento
%     esquema_convec : esquema aceito por espac_disc
%     interp_pressao : interpolacao aceita por pressure_interp
%     numerical_param, physical_param : ver build_case
%
%   SAIDAS
%     u_star, v_star : campos de velocidade preditos
%     aPU   (Nx-1,Ny)    : coeficiente central aP de cada no de u
%     aU    (Nx-1,Ny,4)  : [aE,aW,aN,aS] de cada no de u
%     aPV   (Nx,Ny-1)    : coeficiente central aP de cada no de v
%     aV    (Nx,Ny-1,4)  : [aE,aW,aN,aS] de cada no de v
%     is_wall, is_wall_u, is_wall_v : mascaras de parede solida (VC
%           escalar, no de u, no de v)

U        = physical_param(5);
alpha_u  = numerical_param{1};
Nx       = numerical_param{5};  Ny = numerical_param{6};
dx       = numerical_param{7};  dy = numerical_param{8};

is_wall = numerical_param{11};   % máscara de sólidos nos VCs escalares

% Máscara para nós de u(i,j) - sólido se VC(i,j) ou VC(i+1,j) for parede
is_wall_u = false(Nx-1, Ny);
for i = 1:Nx-1
    for j = 1:Ny
        is_wall_u(i,j) = is_wall(i,j) || is_wall(i+1,j);
    end
end

% Máscara para nós de v(i,j) - sólido se VC(i,j) ou VC(i,j+1) for parede
is_wall_v = false(Nx, Ny-1);
for i = 1:Nx
    for j = 1:Ny-1
        is_wall_v(i,j) = is_wall(i,j) || is_wall(i,j+1);
    end
end

% aP infinito nos nós mascarados anula d: paredes ficam sem fluxo em p'
aPU = Inf(Nx-1, Ny);
aPV = Inf(Nx,   Ny-1);
aU  = zeros(Nx-1, Ny, 4);   % [aE,aW,aN,aS]
aV  = zeros(Nx,   Ny-1, 4);

% ────────────────────────────────────────────────────
%  PASSO 1a - Momento u*
% ────────────────────────────────────────────────────
u_star = u;
u_prev = u;
for i = 1:Nx-1
    for j = 1:Ny
        if is_wall_u(i,j)
            u_star(i,j) = 0;
            continue
        end
        F     = mass_flux("u", escoamento, i, j, u, v, u_analitica, u_star, numerical_param, physical_param);
        [a, aP, Su_dc] = espac_disc("u", escoamento, esquema_convec, D, F, i, j, u_prev, u_analitica, is_wall_u, numerical_param, physical_param);
        s_vc  = boundary_conditions("u", escoamento, i, j, u, u_analitica, a, numerical_param, U);
        src_p = pressure_interp("u", interp_pressao, p, dy, dx, i, j);
        rhs   = s_vc + src_p + Su_dc + (1 - alpha_u)/alpha_u * aP * u(i,j);
        u_star(i,j) = rhs / (aP / alpha_u);

        aPU(i,j)  = aP;
        aU(i,j,:) = a;
    end
end

% ────────────────────────────────────────────────────
%  PASSO 1b - Momento v*
% ────────────────────────────────────────────────────
v_star = v;
v_prev = v;
for i = 1:Nx
    for j = 1:Ny-1
        if is_wall_v(i,j)
            v_star(i,j) = 0;
            continue
        end
        F     = mass_flux("v", escoamento, i, j, u, v, u_analitica, u_star, numerical_param, physical_param);
        [a, aP, Su_dc] = espac_disc("v", escoamento, esquema_convec, D, F, i, j, v_prev, u_analitica, is_wall_v, numerical_param, physical_param);
        % boundary_conditions recebe v_prev no argumento "u" - internamente
        % o case "v" usa esse argumento como campo de velocidade v
        s_vc  = boundary_conditions("v", escoamento, i, j, v_prev, u_analitica, a, numerical_param, U);
        src_p = pressure_interp("v", interp_pressao, p, dy, dx, i, j);
        rhs   = s_vc + src_p + Su_dc + (1 - alpha_u)/alpha_u * aP * v(i,j);
        v_star(i,j) = rhs / (aP / alpha_u);

        aPV(i,j)  = aP;
        aV(i,j,:) = a;
    end
end

end
