function [errL2, errLinf] = couette_error(u, u_analitica, xp, yp, Nx, Ny)

    u_centro = zeros(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            if i == 1
                ue = u(1,j);      uw = u_analitica(yp(j));
            elseif i == Nx
                ue = u(Nx-1,j);   uw = u(Nx-1,j);
            else
                ue = u(i,j);      uw = u(i-1,j);
            end
            u_centro(i,j) = 0.5*(uw + ue);
        end
    end

    u_ref = repmat(u_analitica(yp), Nx, 1);

    errL2   = norm(u_centro - u_ref, 'fro') / norm(u_ref, 'fro');
    errLinf = max(abs(u_centro(:) - u_ref(:))) / max(abs(u_ref(:)));
end
