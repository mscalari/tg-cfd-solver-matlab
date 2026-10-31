function bp = mass_residual(u_star, v_star, u_analitica, i, j, numerical_param, physical_param, escoamento)
%  Resíduo de continuidade integrado no VC escalar (i,j):
%    b'(i,j) = rho * [(uw - ue)*dy + (vs - vn)*dx]
%  Nós mascarados valem 0, então paredes internas têm fluxo nulo.

    rho = physical_param(4);
    Nx  = numerical_param{5};  Ny = numerical_param{6};
    dx  = numerical_param{7};  dy = numerical_param{8};
    yp  = numerical_param{10};

    if i < Nx
        ue_vc = u_star(i,j);
    elseif escoamento == "LidDrivenCavity"
        ue_vc = 0;
    else
        ue_vc = u_star(Nx-1,j);   % saída: gradiente nulo
    end
    if i > 1,  uw_vc = u_star(i-1,j);  else, uw_vc = u_analitica(yp(j)); end
    if j < Ny, vn_vc = v_star(i,j);    else, vn_vc = 0; end
    if j > 1,  vs_vc = v_star(i,j-1);  else, vs_vc = 0; end

    bp = rho * ((uw_vc - ue_vc)*dy + (vs_vc - vn_vc)*dx);
end
