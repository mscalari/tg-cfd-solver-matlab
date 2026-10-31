function [u_centro, v_centro] = cell_center_velocity(c, u, v)
%CELL_CENTER_VELOCITY  Interpola u e v das faces para os centros dos VCs escalares.
%
%   ENTRADAS
%     c    : caso de build_case
%     u, v : campos na malha deslocada
%
%   SAÍDAS
%     u_centro, v_centro (Nx,Ny) : velocidades centradas; NaN nos VCs sólidos
%
%   Faces sobre paredes valem 0, a face de entrada usa c.u_analitica e a de
%   saída repete o último nó de u (gradiente nulo).
    Nx = c.Nx;  Ny = c.Ny;  is_wall = c.is_wall;

    u_centro = NaN(Nx, Ny);
    v_centro = NaN(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            if is_wall(i,j), continue; end

            if i == 1
                uw = c.u_analitica(c.yp(j));
            elseif is_wall(i-1,j)
                uw = 0;
            else
                uw = u(i-1,j);
            end
            if i == Nx
                if c.fechado, ue = 0; else, ue = u(Nx-1,j); end
            elseif is_wall(i+1,j)
                ue = 0;
            else
                ue = u(i,j);
            end
            u_centro(i,j) = 0.5*(uw + ue);

            if j == 1 || is_wall(i,j-1)
                vs = 0;
            else
                vs = v(i,j-1);
            end
            if j == Ny || is_wall(i,j+1)
                vn = 0;
            else
                vn = v(i,j);
            end
            v_centro(i,j) = 0.5*(vs + vn);
        end
    end
end
