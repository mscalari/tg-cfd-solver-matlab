function result = run_case(c, algorithm, esquema_convec)
%RUN_CASE  Resolve um caso com SIMPLE, SIMPLEC ou PISO, partindo de c.u0/v0/p0.
%
%   Os três algoritmos compartilham o laço externo:
%     preditor -> resíduo de massa -> p' -> correção de velocidade.
%   O SIMPLEC só muda a fórmula de d; o PISO acrescenta um 2o corretor com a
%   correção de vizinhança.
%
%   ENTRADAS
%     c              : caso montado por build_case
%     algorithm      : "SIMPLE" | "SIMPLEC" | "PISO"
%     esquema_convec : nome do esquema aceito por espac_disc
%
%   SAÍDA - struct com u, v, p, res_hist, n_iter, t_cpu, res_final, converged

    algorithm = string(algorithm);
    if ~any(algorithm == ["SIMPLE", "SIMPLEC", "PISO"])
        error('run_case: algoritmo desconhecido "%s".', algorithm);
    end

    num  = c.numerical_param;
    phys = c.physical_param;
    alpha_u = num{1};  alpha_p = num{2};  tol = num{3};  max_iter = num{4};
    dx = c.dx;  dy = c.dy;
    u = c.u0;  v = c.v0;  p = c.p0;

    % d no PISO usa a fórmula do SIMPLE: os coeficientes não mudam na iteração
    formula_d = "SIMPLE";
    if algorithm == "SIMPLEC", formula_d = "SIMPLEC"; end

    res_hist = zeros(max_iter, 1);
    t0 = tic;

    for iter = 1:max_iter
        [u_star, v_star, aPU, aU, aPV, aV, ~, wall_u, wall_v] = ...
            momentum_predictor(u, v, p, c.u_analitica, c.D, c.escoamento, ...
                esquema_convec, "Linear", num, phys);

        bp  = residual_field(u_star, v_star, c);
        d_u = compute_d_field(aPU, sum(aU,3), alpha_u, dy, formula_d);
        d_v = compute_d_field(aPV, sum(aV,3), alpha_u, dx, formula_d);

        p_prime = pressure_correc_field(d_u, d_v, num, phys, bp, c.escoamento);
        [u, v]  = correct_velocity(u_star, v_star, d_u, d_v, p_prime, wall_u, wall_v);

        if algorithm == "PISO"
            u_dstar = u;   v_dstar = v;   % u**, v**
            u_til = neighbour_correction("u", c, u_star, u_dstar, aU, aPU, wall_u, alpha_u);
            v_til = neighbour_correction("v", c, v_star, v_dstar, aV, aPV, wall_v, alpha_u);

            bp2 = residual_field(u_til, v_til, c);
            p_prime2 = pressure_correc_field(d_u, d_v, num, phys, bp2, c.escoamento);
            [u, v] = correct_velocity(u_til, v_til, d_u, d_v, p_prime2, wall_u, wall_v);
            p = p + alpha_p * (p_prime + p_prime2);
        else
            p = p + alpha_p * p_prime;
        end

        % Resíduo do preditor: mesma medida nos três algoritmos
        [~, res_hist, flag] = convergence(bp, iter, res_hist, num);
        if flag, break; end
    end

    result.u         = u;
    result.v         = v;
    result.p         = p;
    result.res_hist  = res_hist;
    result.n_iter    = numel(res_hist);
    result.t_cpu     = toc(t0);
    result.res_final = res_hist(end);
    result.converged = ~isnan(result.res_final) && (result.res_final < tol);
end

function bp = residual_field(u, v, c)
%RESIDUAL_FIELD  Resíduo de massa em todos os VCs escalares (0 nos sólidos).
    bp = zeros(c.Nx, c.Ny);
    for i = 1:c.Nx
        for j = 1:c.Ny
            if c.is_wall(i,j), continue; end
            bp(i,j) = mass_residual(u, v, c.u_analitica, i, j, ...
                                    c.numerical_param, c.physical_param, c.escoamento);
        end
    end
end

function [u, v] = correct_velocity(u, v, d_u, d_v, p_prime, wall_u, wall_v)
%CORRECT_VELOCITY  Aplica u' = d*(p'_P - p'_E) e análogo em v, fora dos sólidos.
    for i = 1:size(u,1)
        for j = 1:size(u,2)
            if ~wall_u(i,j)
                u(i,j) = u(i,j) + d_u(i,j)*(p_prime(i,j) - p_prime(i+1,j));
            end
        end
    end
    for i = 1:size(v,1)
        for j = 1:size(v,2)
            if ~wall_v(i,j)
                v(i,j) = v(i,j) + d_v(i,j)*(p_prime(i,j) - p_prime(i,j+1));
            end
        end
    end
end

function f_til = neighbour_correction(direction, c, f_star, f_dstar, aNb, aP, wall, alpha_u)
%NEIGHBOUR_CORRECTION  Termo de vizinhança do 2o corretor do PISO.
%   f_til = f** + sum(anb*(f**_nb - f*_nb)) / (aP/alpha_u), para u ou v.
    U = c.U;
    f_til = f_dstar;
    for i = 1:size(f_star,1)
        for j = 1:size(f_star,2)
            if wall(i,j), continue; end
            a = reshape(aNb(i,j,:), 1, 4);
            s_star  = boundary_conditions(direction, c.escoamento, i, j, f_star,  c.u_analitica, a, c.numerical_param, U);
            s_dstar = boundary_conditions(direction, c.escoamento, i, j, f_dstar, c.u_analitica, a, c.numerical_param, U);
            f_til(i,j) = f_dstar(i,j) + (s_dstar - s_star) / (aP(i,j)/alpha_u);
        end
    end
end
