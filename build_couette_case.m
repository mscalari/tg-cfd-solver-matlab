function [numerical_param, physical_param, D, u_analitica, xp, yp] = ...
    build_couette_case(Nx, Ny, Lx, Ly, mu, rho, U, G, alpha_u, alpha_p, tol, max_iter)
%   Monta numerical_param/physical_param/D/u_analitica
%   para um caso de Couette-Poiseuille, evitando repetir esse bloco em
%   cada estudo da suíte comparativa.

    dx = Lx/Nx;  dy = Ly/Ny;
    xp = ((1:Nx)-0.5)'*dx;
    yp = ((1:Ny)-0.5)*dy;

    physical_param = [Lx, Ly, mu, rho, U, G];
    [~, u_analitica] = analytical_Couette_Pouiseuille(physical_param);
    D = dif_coeff(mu, dx, dy);

    numerical_param = {alpha_u, alpha_p, tol, max_iter, Nx, Ny, dx, dy, xp, yp};
end
