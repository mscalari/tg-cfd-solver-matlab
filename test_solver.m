%% TEST_SOLVER  Verificações rápidas do solver. Rodar a partir desta pasta.

% Couette-Poiseuille: partindo do repouso, os 3 algoritmos chegam à analítica
% e o erro cai com 2a ordem ao refinar em y.
Lx = 2; Ly = 1; mu = 0.1; rho = 1; U = 2; G = 4;
Ny_vals = [10 20];
for k = 1:3
    alg = ["SIMPLE", "SIMPLEC", "PISO"];
    err = zeros(size(Ny_vals));
    for n = 1:numel(Ny_vals)
        Ny = Ny_vals(n); Nx = 10;
        [num, phys, D, ua, xp, yp] = build_couette_case(Nx, Ny, Lx, Ly, mu, rho, U, G, 0.7, 0.3, 1e-7, 20000);
        r = run_case(alg(k), zeros(Nx-1,Ny), zeros(Nx,Ny-1), zeros(Nx,Ny), ua, D, ...
            "CouettePoiseuille", 'Upwind 1a ordem', "Linear", num, phys);
        assert(r.converged, '%s nao convergiu', alg(k));
        err(n) = couette_error(r.u, ua, xp, yp, Nx, Ny);
    end
    fprintf('%-7s erros %.3e %.3e  ordem %.2f\n', alg(k), err, log2(err(1)/err(2)));
    assert(log2(err(1)/err(2)) > 1.7, 'ordem de Couette abaixo de 2 em %s', alg(k));
end

% Cavidade aberta: massa que entra sai, inclusive com esquema de alta ordem.
Nx = 24; Ny = 12; Lx = 2; Ly = 1; dx = Lx/Nx; dy = Ly/Ny; mu = 0.02;
xp = ((1:Nx)-0.5)'*dx; yp = ((1:Ny)-0.5)*dy;
ic_s = find(xp > 0.5, 1, 'first'); ic_e = find(xp <= 1, 1, 'last'); jc_t = find(yp <= 0.5, 1, 'last');
iw = false(Nx,Ny); iw([1:ic_s-1, ic_e+1:Nx], 1:jc_t) = true;
ua = @(y) 6*((y-0.5)/0.5).*(1-(y-0.5)/0.5).*(y > 0.5);
num = {0.5, 0.5, 1e-6, 3000, Nx, Ny, dx, dy, xp, yp, iw};
phys = [Lx, Ly, mu, 1, 1, 0];
r = run_case("SIMPLEC", repmat(ua(yp),Nx-1,1), zeros(Nx,Ny-1), zeros(Nx,Ny), ua, dif_coeff(mu,dx,dy), ...
    "FlowOverCavity", 'QUICK', "Linear", num, phys);
assert(r.converged, 'cavidade nao convergiu');
Q_in = sum(ua(yp))*dy; Q_out = sum(r.u(end,:))*dy;
fprintf('Cavidade: iter=%d Q_in=%.6f Q_out=%.6f\n', r.n_iter, Q_in, Q_out);
assert(abs(Q_in - Q_out) < 1e-4, 'massa nao conservada');
assert(any(r.u(ic_s:ic_e-1, 1:jc_t) < 0, 'all'), 'cavidade sem recirculacao');

fprintf('OK\n');
