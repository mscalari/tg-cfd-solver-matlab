function [u, v, p, res_hist] = piso_solver( ...
    u, v, p, ...
    u_analitica, ...
    D, ...
    escoamento, esquema_convec, ...
    interp_pressao, ...
    numerical_param, physical_param)

U        = physical_param(5);
alpha_u  = numerical_param{1};  alpha_p = numerical_param{2};
max_iter = numerical_param{4};
Nx       = numerical_param{5};  Ny      = numerical_param{6};
dx       = numerical_param{7};  dy      = numerical_param{8};

res_hist = zeros(max_iter, 1);

for iter = 1:max_iter

    % ────────────────────────────────────────────────────
    %  PREDITOR — u*, v* + coeficientes de momento (aP, aE, aW, aN, aS)
    % ────────────────────────────────────────────────────
    [u_star, v_star, aPU, aU, aPV, aV, is_wall, is_wall_u, is_wall_v] = ...
        momentum_predictor(u, v, p, u_analitica, D, escoamento, ...
            esquema_convec, interp_pressao, numerical_param, physical_param);

    % d usado nos DOIS corretores — fórmula do SIMPLE, pois os
    % coeficientes de momento não mudam dentro da mesma iteração externa
    d_u = compute_d_field(aPU, sum(aU,3), alpha_u, dy, "SIMPLE");
    d_v = compute_d_field(aPV, sum(aV,3), alpha_u, dx, "SIMPLE");

    % ────────────────────────────────────────────────────
    %  CORRETOR 1 — idêntico ao SIMPLE
    % ────────────────────────────────────────────────────
    bp1 = zeros(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            if is_wall(i,j), continue; end
            bp1(i,j) = mass_residual(u_star, v_star, u_analitica, i, j, ...
                                      numerical_param, physical_param, escoamento);
        end
    end
    p_prime = pressure_correc_field(d_u, d_v, numerical_param, physical_param, bp1, escoamento);

    u_dstar = u_star;   % u**
    for i = 1:Nx-1
        for j = 1:Ny
            if ~is_wall_u(i,j)
                u_dstar(i,j) = u_star(i,j) + d_u(i,j)*(p_prime(i,j) - p_prime(i+1,j));
            end
        end
    end
    v_dstar = v_star;   % v**
    for i = 1:Nx
        for j = 1:Ny-1
            if ~is_wall_v(i,j)
                v_dstar(i,j) = v_star(i,j) + d_v(i,j)*(p_prime(i,j) - p_prime(i,j+1));
            end
        end
    end

    % ────────────────────────────────────────────────────
    %  CORRETOR 2 — correção de vizinhança (o passo extra do PISO)
    % ────────────────────────────────────────────────────
    u_til = u_dstar;
    for i = 1:Nx-1
        for j = 1:Ny
            if is_wall_u(i,j), continue; end
            a = reshape(aU(i,j,:), 1, 4);
            s_vc_star  = boundary_conditions("u", escoamento, i, j, u_star,  u_analitica, a, numerical_param, U);
            s_vc_dstar = boundary_conditions("u", escoamento, i, j, u_dstar, u_analitica, a, numerical_param, U);
            u_til(i,j) = u_dstar(i,j) + (s_vc_dstar - s_vc_star) / (aPU(i,j)/alpha_u);
        end
    end
    v_til = v_dstar;
    for i = 1:Nx
        for j = 1:Ny-1
            if is_wall_v(i,j), continue; end
            a = reshape(aV(i,j,:), 1, 4);
            s_vc_star  = boundary_conditions("v", escoamento, i, j, v_star,  u_analitica, a, numerical_param, U);
            s_vc_dstar = boundary_conditions("v", escoamento, i, j, v_dstar, u_analitica, a, numerical_param, U);
            v_til(i,j) = v_dstar(i,j) + (s_vc_dstar - s_vc_star) / (aPV(i,j)/alpha_u);
        end
    end

    bp2 = zeros(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            if is_wall(i,j), continue; end
            bp2(i,j) = mass_residual(u_til, v_til, u_analitica, i, j, ...
                                      numerical_param, physical_param, escoamento);
        end
    end
    p_prime2 = pressure_correc_field(d_u, d_v, numerical_param, physical_param, bp2, escoamento);

    u_triple = u_til;
    for i = 1:Nx-1
        for j = 1:Ny
            if ~is_wall_u(i,j)
                u_triple(i,j) = u_til(i,j) + d_u(i,j)*(p_prime2(i,j) - p_prime2(i+1,j));
            end
        end
    end
    v_triple = v_til;
    for i = 1:Nx
        for j = 1:Ny-1
            if ~is_wall_v(i,j)
                v_triple(i,j) = v_til(i,j) + d_v(i,j)*(p_prime2(i,j) - p_prime2(i,j+1));
            end
        end
    end

    % ────────────────────────────────────────────────────
    %  ATUALIZAR CAMPOS
    % ────────────────────────────────────────────────────
    p = p + alpha_p * (p_prime + p_prime2);
    u = u_triple;
    v = v_triple;

    % ────────────────────────────────────────────────────
    %  CONVERGÊNCIA — resíduo do preditor, mesma medida do SIMPLE/SIMPLEC
    % ────────────────────────────────────────────────────
    [~, res_hist, flag] = convergence(bp1, iter, res_hist, numerical_param);
    if flag, break; end
end

end
