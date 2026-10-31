function [errL2, errLinf] = couette_error(c, r)
%COUETTE_ERROR  Erros relativos L2 e Linf de Couette-Poiseuille frente à analítica.
%
%   ENTRADAS
%     c : caso de build_case
%     r : resultado de run_case
%
%   SAÍDAS
%     errL2, errLinf : erros de u nos centros dos VCs, relativos à analítica
    u_centro = cell_center_velocity(c, r.u, r.v);
    u_ref = repmat(c.u_analitica(c.yp), c.Nx, 1);

    errL2   = norm(u_centro - u_ref, 'fro') / norm(u_ref, 'fro');
    errLinf = max(abs(u_centro(:) - u_ref(:))) / max(abs(u_ref(:)));
end
