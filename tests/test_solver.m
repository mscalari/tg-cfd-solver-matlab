%TEST_SOLVER  Verificações rápidas do solver. Rodar com: cfd_run test

% Couette-Poiseuille: partindo do repouso, os 3 algoritmos chegam à analítica
% e o erro cai com 2a ordem ao refinar em y.
Ny_vals = [10 20];
alg = ["SIMPLE", "SIMPLEC", "PISO"];
for k = 1:3
    err = zeros(size(Ny_vals));
    for n = 1:numel(Ny_vals)
        c = build_case("CouettePoiseuille", 'Lx', 2, 'Nx', 10, 'Ny', Ny_vals(n), ...
            'alpha_u', 0.7, 'alpha_p', 0.3, 'tol', 1e-7, 'max_iter', 20000);
        r = run_case(c, alg(k), 'Upwind 1a ordem');
        assert(r.converged, '%s nao convergiu', alg(k));
        err(n) = couette_error(c, r);
    end
    fprintf('%-7s erros %.3e %.3e  ordem %.2f\n', alg(k), err, log2(err(1)/err(2)));
    assert(log2(err(1)/err(2)) > 1.7, 'ordem de Couette abaixo de 2 em %s', alg(k));
end

% Cavidade aberta: massa que entra sai, inclusive com esquema de alta ordem.
c = build_case("FlowOverCavity", 'Nx', 24, 'Ny', 12, 'Lx_ante', 0.5, 'L_cav', 0.5, ...
    'Lx_post', 1, 'Ly_canal', 0.5, 'D_cav', 0.5, 'Re', 25, 'tol', 1e-6, 'max_iter', 3000);
r = run_case(c, "SIMPLEC", 'QUICK');
assert(r.converged, 'cavidade nao convergiu');
Q_in = sum(c.u_analitica(c.yp))*c.dy; Q_out = sum(r.u(end,:))*c.dy;
fprintf('Cavidade: iter=%d Q_in=%.6f Q_out=%.6f\n', r.n_iter, Q_in, Q_out);
assert(abs(Q_in - Q_out) < 1e-4, 'massa nao conservada');
assert(any(r.u(c.ic_s:c.ic_e-1, 1:c.jc_t) < 0, 'all'), 'cavidade sem recirculacao');

% Entradas inválidas precisam falhar com erro, não rodar um caso errado.
falhou = false;
try, build_case("LidDrivenCavity", 'Nxx', 10); catch, falhou = true; end
assert(falhou, 'build_case aceitou parametro inexistente');
falhou = false;
try, run_case(c, "FOO", 'QUICK'); catch, falhou = true; end
assert(falhou, 'run_case aceitou algoritmo inexistente');

fprintf('OK\n');
