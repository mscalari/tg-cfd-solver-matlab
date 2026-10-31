function [a, aP, Su_dc] = espac_disc(direction, escoamento, esquema_convec, D, F, i, j, u_prev, u_analitica, mask, numerical_param, physical_param)
%ESPAC_DISC  Coeficientes convectivo-difusivos de um nó de velocidade.
%
%   ENTRADAS
%     direction      : "u" ou "v"
%     esquema_convec : 'Upwind 1a ordem' | 'Diferenças Centrais' | 'Power Law'
%                      | 'Upwind 2a ordem' | 'QUICK' | 'MUSCL' | 'HRIC'
%     D              : condutâncias de dif_coeff
%     F              : fluxos [Fe,Fw,Fn,Fs] de mass_flux
%     u_prev         : campo da componente na iteração anterior
%     mask           : nós dessa componente sobre sólidos
%
%   SAÍDAS
%     a     : [aE,aW,aN,aS]
%     aP    : coeficiente central, sem sub-relaxação
%     Su_dc : correção diferida. Os esquemas de alta ordem mantêm na matriz os
%             coeficientes de upwind (estáveis) e levam a diferença para o RHS,
%             avaliada com u_prev. Fora deles, Su_dc = 0.
    Nx = numerical_param{5};  Ny = numerical_param{6};
    dx = numerical_param{7};  dy = numerical_param{8};
    yp = numerical_param{10};

    fechado = escoamento == "LidDrivenCavity";
    U_topo  = physical_param(5) * any(escoamento == ["CouettePoiseuille", "LidDrivenCavity"]);

    if direction == "u"
        De_u = D(1); Dw_u = D(2); Dn_u = D(3); Ds_u = D(4);
        % Parede tangencial fica a dy/2 do nó: condutância dobra
        if j == Ny || mask(i,j+1), Dn_u = 2*Dn_u; end
        if j == 1  || mask(i,j-1), Ds_u = 2*Ds_u; end
    elseif direction == "v"
        De_u = D(5); Dw_u = D(6); Dn_u = D(7); Ds_u = D(8);
        if i == 1 || mask(i-1,j), Dw_u = 2*Dw_u; end
        if (i == Nx && fechado) || (i < Nx && mask(i+1,j)), De_u = 2*De_u; end
    end

    Fe = F(1); Fw = F(2); Fn = F(3); Fs = F(4);
    Su_dc = 0;

    switch esquema_convec
        case 'Diferenças Centrais'
            aE = De_u - Fe/2;
            aW = Dw_u + Fw/2;
            aN = Dn_u - Fn/2;
            aS = Ds_u + Fs/2;

        case 'Power Law'
            pl_A = @(Pe) max(0, (1 - 0.1*abs(Pe))^5);
            aE = De_u * pl_A(Fe / De_u) + max(-Fe, 0);
            aW = Dw_u * pl_A(Fw / Dw_u) + max( Fw, 0);
            aN = Dn_u * pl_A(Fn / Dn_u) + max(-Fn, 0);
            aS = Ds_u * pl_A(Fs / Ds_u) + max( Fs, 0);

        otherwise
            aE = De_u + max(-Fe, 0);
            aW = Dw_u + max( Fw, 0);
            aN = Dn_u + max(-Fn, 0);
            aS = Ds_u + max( Fs, 0);
    end

    % Esquemas de correção diferida precisam de dois vizinhos em cada direção
    if ~any(strcmp(esquema_convec, {'Upwind 1a ordem', 'Diferenças Centrais', 'Power Law'}))
        uP = u_prev(i,j);
        if direction == "u"
            if fechado, tipo_E = "no"; else, tipo_E = "saida"; end
            [uW,uWW] = vizinhos(u_prev, mask, i-1,j, i-2,j, uP, "no",    u_analitica(yp(j)), "no");
            [uE,uEE] = vizinhos(u_prev, mask, i+1,j, i+2,j, uP, tipo_E,  0,                  "no");
            [uN,uNN] = vizinhos(u_prev, mask, i,j+1, i,j+2, uP, "ghost", U_topo,             "ghost");
            [uS,uSS] = vizinhos(u_prev, mask, i,j-1, i,j-2, uP, "ghost", 0,                  "ghost");
        else
            if fechado, tipo_E = "ghost"; else, tipo_E = "saida"; end
            [uW,uWW] = vizinhos(u_prev, mask, i-1,j, i-2,j, uP, "ghost", 0, "ghost");
            [uE,uEE] = vizinhos(u_prev, mask, i+1,j, i+2,j, uP, tipo_E,  0, "ghost");
            [uN,uNN] = vizinhos(u_prev, mask, i,j+1, i,j+2, uP, "no",    0, "no");
            [uS,uSS] = vizinhos(u_prev, mask, i,j-1, i,j-2, uP, "no",    0, "no");
        end

        % Peso do HRIC: alinhamento do gradiente local com cada direção
        grad_x = (uE - uW) / (2*dx);
        grad_y = (uN - uS) / (2*dy);
        grad_mag = sqrt(grad_x^2 + grad_y^2) + 1e-14;
        cos2_x = (grad_x / grad_mag)^2;
        cos2_y = (grad_y / grad_mag)^2;

        % Nós [C A U] de cada face: fluxo positivo, depois negativo.
        % Faces leste/norte entram com -F; oeste/sul, com +F.
        Su_dc = face_correction(esquema_convec, -Fe, Fe >= 0, [uP uE uW], [uE uP uEE], cos2_x) ...
              + face_correction(esquema_convec,  Fw, Fw >= 0, [uW uP uWW], [uP uW uE], cos2_x) ...
              + face_correction(esquema_convec, -Fn, Fn >= 0, [uP uN uS], [uN uP uNN], cos2_y) ...
              + face_correction(esquema_convec,  Fs, Fs >= 0, [uS uP uSS], [uP uS uN], cos2_y);
    end

    a  = [aE, aW, aN, aS];
    aP = aE + aW + aN + aS + (Fe - Fw + Fn - Fs);
end

function S = face_correction(esquema_convec, F_s, positivo, nos_pos, nos_neg, cos2)
%FACE_CORRECTION  Correção diferida de uma face, relativa ao upwind de 1a ordem.
%
%   F_s      : fluxo da face já com o sinal da sua contribuição
%   positivo : true se o fluxo físico da face é >= 0
%   nos_pos, nos_neg : nós [C A U] de cada sentido do fluxo, onde
%                      C = central (upwind), A = a jusante, U = dois a montante
%   cos2     : quadrado do cosseno entre o gradiente e a direção (HRIC)
%
%   LUD   : phi_f = C + (C-U)/2
%   QUICK : phi_f = C + (3A - 2C - U)/8
%   MUSCL : kappa = 1/3 com limitador de van Albada
%   HRIC  : mistura do valor limitado e de diferenças centrais, pesada por cos2
    if positivo, nos = nos_pos; else, nos = nos_neg; end
    C = nos(1);  A = nos(2);  U = nos(3);

    switch esquema_convec
        case 'Upwind 2a ordem'
            S = (F_s/2) * (C - U);

        case 'QUICK'
            S = (F_s/8) * (3*A - 2*C - U);

        case 'MUSCL'
            d = A - C;
            if abs(d) < 1e-14, r = 0; else, r = (C - U) / d; end
            psi = (r > 0) .* (r.^2 + r) ./ (r.^2 + 1);   % psi = 0 em extremos
            S = (F_s/6) * psi * ((C - U) + 2*(A - C));

        case 'HRIC'
            dDE = A - U;
            if abs(dDE) < 1e-14
                S = 0;
                return
            end
            phiC = (C - U) / dDE;
            if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
            elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
            else,                             phif_star = phiC;
            end
            phif_cd = (1 + phiC) / 2;
            phif = phif_star * cos2 + phif_cd * (1 - cos2);
            S = F_s * (phif - phiC) * dDE;

        otherwise
            error('espac_disc: esquema desconhecido "%s".', esquema_convec);
    end
end

function [n1, n2] = vizinhos(f, mask, a, b, c, d, uP, tipo, w, tipo_mask)
%VIZINHOS  Valores no 1o (a,b) e 2o (c,d) nó vizinhos, com contorno de valor w.
%   tipo vale para o contorno do domínio; tipo_mask, para paredes internas:
%     "no"    : contorno sobre um nó da malha deslocada
%     "ghost" : parede a meia célula (nó fantasma espelhado)
%     "saida" : gradiente nulo
    [Mx, My] = size(f);
    dentro = @(p,q) p >= 1 && p <= Mx && q >= 1 && q <= My;

    if dentro(a,b) && ~mask(a,b)
        n1 = f(a,b);
        if dentro(c,d) && ~mask(c,d)
            n2 = f(c,d);
            return
        end
        if dentro(c,d), tipo = tipo_mask; w = 0; end
        switch tipo
            case "no",    n2 = w;
            case "ghost", n2 = 2*w - n1;
            case "saida", n2 = n1;
        end
        return
    end

    if dentro(a,b), tipo = tipo_mask; w = 0; end
    switch tipo
        case "no",    n1 = w;         n2 = 2*w - uP;
        case "ghost", n1 = 2*w - uP;  n2 = n1;   % só alimenta faces de fluxo nulo
        case "saida", n1 = uP;        n2 = uP;
    end
end
