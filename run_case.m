function result = run_case(algorithm, u0, v0, p0, u_analitica, D, ...
    escoamento, esquema_convec, interp_pressao, numerical_param, physical_param)
%RUN_CASE  Roda o algoritmo pedido (SIMPLE/SIMPLEC/PISO), cronometra e
%   empacota os resultados num struct único - usado por todos os estudos
%   da suíte comparativa para evitar repetir o switch/tic/toc em cada um.

    tol = numerical_param{3};

    t0 = tic;
    switch algorithm
        case "SIMPLE"
            [u, v, p, res_hist] = simple_solver(u0, v0, p0, u_analitica, D, ...
                escoamento, esquema_convec, interp_pressao, numerical_param, physical_param);
        case "SIMPLEC"
            [u, v, p, res_hist] = simplec_solver(u0, v0, p0, u_analitica, D, ...
                escoamento, esquema_convec, interp_pressao, numerical_param, physical_param);
        case "PISO"
            [u, v, p, res_hist] = piso_solver(u0, v0, p0, u_analitica, D, ...
                escoamento, esquema_convec, interp_pressao, numerical_param, physical_param);
        otherwise
            error('run_case: algoritmo desconhecido "%s".', algorithm);
    end
    t_cpu = toc(t0);

    result.u         = u;
    result.v         = v;
    result.p         = p;
    result.res_hist  = res_hist;
    result.n_iter    = numel(res_hist);
    result.t_cpu     = t_cpu;
    result.res_final = res_hist(end);
    result.converged = ~isnan(result.res_final) && (result.res_final < tol);
end
