function F = mass_flux(direction, escoamento, i, j, u, v, u_analitica, u_star, numerical_param, physical_param)
%MASS_FLUX  Fluxos de massa [Fe,Fw,Fn,Fs] nas faces do VC deslocado.
%   As faces ficam entre nós de velocidade, logo o fluxo normal é a média
%   de dois nós. Nós mascarados valem 0, então paredes internas não exigem
%   tratamento próprio.
    rho = physical_param(4);
    Nx  = numerical_param{5};  Ny = numerical_param{6};
    dx  = numerical_param{7};  dy = numerical_param{8};
    yp  = numerical_param{10};

    fechado = escoamento == "LidDrivenCavity";

    if direction == "u"
        if i < Nx-1
            ue = 0.5*(u(i,j) + u(i+1,j));
        elseif fechado
            ue = 0.5*u(i,j);
        else
            ue = u(i,j);   % saída: gradiente nulo
        end
        if i == 1
            uw = 0.5*(u_analitica(yp(j)) + u(1,j));
        else
            uw = 0.5*(u(i-1,j) + u(i,j));
        end
        if j == Ny
            vn = 0;
        else
            vn = 0.5*(v(i,j) + v(i+1,j));
        end
        if j == 1
            vs = 0;
        else
            vs = 0.5*(v(i,j-1) + v(i+1,j-1));
        end
        Fe = rho*ue*dy;  Fw = rho*uw*dy;
        Fn = rho*vn*dx;  Fs = rho*vs*dx;

    elseif direction == "v"
        if i < Nx
            ue_v = 0.5*(u_star(i,j) + u_star(i,j+1));
        elseif fechado
            ue_v = 0;
        else
            ue_v = 0.5*(u_star(Nx-1,j) + u_star(Nx-1,j+1));
        end
        if i == 1
            uw_v = 0.5*(u_analitica(yp(j)) + u_analitica(yp(j+1)));
        else
            uw_v = 0.5*(u_star(i-1,j) + u_star(i-1,j+1));
        end
        if j == Ny-1
            vn_v = 0.5*v(i,j);
        else
            vn_v = 0.5*(v(i,j) + v(i,j+1));
        end
        if j == 1
            vs_v = 0.5*v(i,j);
        else
            vs_v = 0.5*(v(i,j-1) + v(i,j));
        end
        Fe = rho*ue_v*dy;  Fw = rho*uw_v*dy;
        Fn = rho*vn_v*dx;  Fs = rho*vs_v*dx;
    end

    F = [Fe, Fw, Fn, Fs];
end
