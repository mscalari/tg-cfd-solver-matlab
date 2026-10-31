function p_prime = pressure_correc_field(d_u, d_v, numerical_param, physical_param, bp, escoamento)
% Resolve p' por Gauss-Seidel. d = 0 nos nós mascarados isola as paredes.

    rho = physical_param(4);
    Nx = numerical_param{5}; Ny = numerical_param{6};
    dx = numerical_param{7}; dy = numerical_param{8};
    fechado = escoamento == "LidDrivenCavity";

    p_prime = zeros(Nx, Ny);

    for gs = 1:200   % iterações internas de Gauss-Seidel para p'

        p_old = p_prime;

        for i = 1:Nx
            for j = 1:Ny

                aP_p  = 0;
                rhs_p = bp(i,j);

                % Leste — face entre P(i,j) e E(i+1,j), usa d_u do nó (i,j)
                if i < Nx
                    aE_p  = rho * d_u(i,j) * dy;
                    aP_p  = aP_p + aE_p;
                    rhs_p = rhs_p + aE_p * p_prime(i+1,j);
                elseif ~fechado
                    % Saída: p' = 0 na face x = Lx, a dx/2 do centro
                    aP_p = aP_p + 2 * rho * d_u(Nx-1,j) * dy;
                end
                % Oeste — face entre W(i-1,j) e P(i,j), usa d_u do nó (i-1,j)
                if i > 1
                    aW_p  = rho * d_u(i-1,j) * dy;
                    aP_p  = aP_p + aW_p;
                    rhs_p = rhs_p + aW_p * p_prime(i-1,j);
                end
                % Norte — usa d_v do nó (i,j)
                if j < Ny
                    aN_p  = rho * d_v(i,j) * dx;
                    aP_p  = aP_p + aN_p;
                    rhs_p = rhs_p + aN_p * p_prime(i,j+1);
                end
                % Sul — usa d_v do nó (i,j-1)
                if j > 1
                    aS_p  = rho * d_v(i,j-1) * dx;
                    aP_p  = aP_p + aS_p;
                    rhs_p = rhs_p + aS_p * p_prime(i,j-1);
                end

                if aP_p > 0
                    p_prime(i,j) = rhs_p / aP_p;
                end

            end
        end

        % Só o domínio fechado precisa de referência; os abertos a têm na saída
        if fechado, p_prime = p_prime - p_prime(1,1); end

        if norm(p_prime - p_old, 'fro') < 1e-12, break; end

    end
end
