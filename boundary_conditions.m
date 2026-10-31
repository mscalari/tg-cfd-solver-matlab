function s_vc = boundary_conditions(direction, escoamento, i, j, u, u_analitica, a, numerical_param, U)
%BOUNDARY_CONDITIONS  Soma aE*uE + aW*uW + aN*uN + aS*uS com os contornos.
%   'u' é o campo da componente em 'direction'. Nós mascarados valem 0,
%   então paredes internas não exigem tratamento próprio.
    aE = a(1);  aW = a(2);  aN = a(3);  aS = a(4);

    Nx  = numerical_param{5};  Ny = numerical_param{6};
    yp  = numerical_param{10};

    fechado = escoamento == "LidDrivenCavity";
    U_topo  = U * any(escoamento == ["CouettePoiseuille", "LidDrivenCavity"]);

    if direction == "u"
        if j == 1,    uS = 0;       else, uS = u(i,j-1); end
        if j == Ny,   uN = U_topo;  else, uN = u(i,j+1); end
        if i == 1,    uW = u_analitica(yp(j)); else, uW = u(i-1,j); end
        if i < Nx-1
            uE = u(i+1,j);
        elseif fechado
            uE = 0;
        else
            uE = u(i,j);   % saída: gradiente nulo
        end

    elseif direction == "v"
        if j == 1,    uS = 0;  else, uS = u(i,j-1); end
        if j == Ny-1, uN = 0;  else, uN = u(i,j+1); end
        if i == 1,    uW = 0;  else, uW = u(i-1,j); end
        if i < Nx
            uE = u(i+1,j);
        elseif fechado
            uE = 0;
        else
            uE = u(i,j);   % saída: gradiente nulo
        end
    end

    s_vc = aW*uW + aE*uE + aN*uN + aS*uS;
end
