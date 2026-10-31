function d = compute_d_field(aP, sumnb, alpha_u, A, algorithm)
% Coeficiente d da correção de velocidade
% u' = d*(p'_W - p'_P), calculado nó a nó.

    aP_relax = aP / alpha_u;

    switch algorithm
        case "SIMPLE"
            d = A ./ aP_relax;
        case "SIMPLEC"
            denom      = aP_relax - sumnb;
            denom_min  = 0.1 * aP_relax;            % limiar seguro (>0 sempre)
            unsafe     = denom < denom_min;
            denom_safe = denom;
            denom_safe(unsafe) = aP_relax(unsafe);  % fallback pontual = fórmula SIMPLE
            d = A ./ denom_safe;
        otherwise
            error('compute_d_field: algorithm deve ser "SIMPLE" ou "SIMPLEC".');
    end
end
