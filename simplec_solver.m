function [u, v, p, res_hist] = simplec_solver( ...
    u, v, p, ...
    u_analitica, ...
    D, ...
    escoamento, esquema_convec, ...
    interp_pressao, ...
    numerical_param, physical_param)

alpha_u  = numerical_param{1};  alpha_p = numerical_param{2};
max_iter = numerical_param{4};
Nx       = numerical_param{5};  Ny      = numerical_param{6};
dx       = numerical_param{7};  dy      = numerical_param{8};

res_hist = zeros(max_iter, 1);

for iter = 1:max_iter

    % ────────────────────────────────────────────────────
    %  PASSO 1 — Preditor (u*, v*) + coeficientes de momento
    % ────────────────────────────────────────────────────
    [u_star, v_star, aPU, aU, aPV, aV, is_wall, is_wall_u, is_wall_v] = ...
        momentum_predictor(u, v, p, u_analitica, D, escoamento, ...
            esquema_convec, interp_pressao, numerical_param, physical_param);

    % ────────────────────────────────────────────────────
    %  PASSO 2 — Resíduo de massa b'
    % ────────────────────────────────────────────────────
    bp = zeros(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            if is_wall(i,j)
                bp(i,j) = 0;
                continue
            end
            bp(i,j) = mass_residual(u_star, v_star, u_analitica, i, j, ...
                                    numerical_param, physical_param, escoamento);
        end
    end

    % ────────────────────────────────────────────────────
    %  PASSO 3 — Campo d (fórmula SIMPLEC) e correção de pressão p'
    % ────────────────────────────────────────────────────
    sumnb_u = sum(aU, 3);   % aE+aW+aN+aS por nó de u
    sumnb_v = sum(aV, 3);   % aE+aW+aN+aS por nó de v

    d_u = compute_d_field(aPU, sumnb_u, alpha_u, dy, "SIMPLEC");
    d_v = compute_d_field(aPV, sumnb_v, alpha_u, dx, "SIMPLEC");

    p_prime = pressure_correc_field(d_u, d_v, numerical_param, physical_param, bp, escoamento);

    % ────────────────────────────────────────────────────
    %  PASSO 4 — Corrigir p, u, v
    % ────────────────────────────────────────────────────
    p = p + alpha_p * p_prime;
    for i = 1:Nx-1
        for j = 1:Ny
            if ~is_wall_u(i,j)
                u(i,j) = u_star(i,j) + d_u(i,j)*(p_prime(i,j) - p_prime(i+1,j));
            end
        end
    end
    for i = 1:Nx
        for j = 1:Ny-1
            if ~is_wall_v(i,j)
                v(i,j) = v_star(i,j) + d_v(i,j)*(p_prime(i,j) - p_prime(i,j+1));
            end
        end
    end

    % ────────────────────────────────────────────────────
    %  CONVERGÊNCIA
    % ────────────────────────────────────────────────────
    [~, res_hist, flag] = convergence(bp, iter, res_hist, numerical_param);
    if flag, break; end
end

end
