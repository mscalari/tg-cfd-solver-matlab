function d = compute_d_field(aP, sumnb, alpha_u, A, algorithm)
%COMPUTE_D_FIELD  Coeficiente d da correção de velocidade, nó a nó.
%   u' = d*(p'_W - p'_P).
%
%   ENTRADAS
%     aP        : coeficiente central de momento, sem sub-relaxação
%     sumnb     : soma dos coeficientes vizinhos (só usada no SIMPLEC)
%     alpha_u   : sub-relaxação de velocidade
%     A         : área da face (dy para u, dx para v)
%     algorithm : "SIMPLE" (d = A/aP) ou "SIMPLEC" (d = A/(aP - sum(anb)))
%
%   No SIMPLEC, se o denominador ficar abaixo de 10% de aP, o nó volta para a
%   fórmula do SIMPLE.

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
