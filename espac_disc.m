function [a, aP, Su_dc] = espac_disc(direction, escoamento, esquema_convec, D, F, i, j, u_prev, u_analitica, mask, numerical_param, physical_param)
%ESPAC_DISC  Coeficientes convectivo-difusivos de um nó de velocidade.
%   direction = "u" ou "v"; mask marca os nós dessa componente sobre sólidos.
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
    end

    Fe = F(1); Fw = F(2); Fn = F(3); Fs = F(4);

    switch esquema_convec
        case 'Upwind 1a ordem'
        % 1ª ordem
        % aW = Dw + max(Fw,0),  aE = De + max(-Fe,0)
    
            aE = De_u + max(-Fe, 0);
            aW = Dw_u + max( Fw, 0);
            aN = Dn_u + max(-Fn, 0);
            aS = Ds_u + max( Fs, 0);
            Su_dc = 0;

            a = [aE,aW,aN,aS];

        case 'Diferenças Centrais'

            aE = De_u - Fe/2;
            aW = Dw_u + Fw/2;
            aN = Dn_u - Fn/2;
            aS = Ds_u + Fs/2;
            Su_dc = 0;

            a = [aE,aW,aN,aS];
         
        case 'Power Law'
            pl_A = @(Pe) max(0, (1 - 0.1*abs(Pe))^5); % Função auxiliar para power-law

            Pee = Fe / De_u;   Pew = Fw / Dw_u;
            Pen = Fn / Dn_u;   Pes = Fs / Ds_u;
    
            aE = De_u * pl_A(Pee) + max(-Fe, 0);
            aW = Dw_u * pl_A(Pew) + max( Fw, 0);
            aN = Dn_u * pl_A(Pen) + max(-Fn, 0);
            aS = Ds_u * pl_A(Pes) + max( Fs, 0);
            Su_dc = 0;

            a = [aE,aW,aN,aS];
        
        case 'Upwind 2a ordem'
            % LUD — Linear Upwind Differencing, 2ª ordem
            %
            % Extrapolação linear com dois nós upstream:
            %   Fe>0: phi_e = phi_P + (1/2)*(phi_P - phi_W)
            %   Fe<0: phi_e = phi_E + (1/2)*(phi_E - phi_EE)
            %
            % Coeficientes principais = upwind1 (estável).
            % Correção diferida S_DC carrega o termo extra de 2a ordem,
            % avaliado com u_prev (valores da iteração anterior).
            %
            % Derivação do S_DC para face leste (Fe>0):
            %   phi_e_LUD  = uP + (1/2)*(uP - uW)
            %   phi_e_UD   = uP
            %   diferença  = (1/2)*(uP - uW)
            %   S_east = -Fe * (1/2)*(uP - uW)   [entra no RHS]
            %
            % Para face oeste (Fw>0) — phi_w entra com sinal +Fw:
            %   phi_w_LUD  = uW + (1/2)*(uW - uWW)
            %   phi_w_UD   = uW
            %   S_west = +Fw * (1/2)*(uW - uWW)
            %
            % Análogo para norte/sul com Fn, Fs.
            % Fronteiras: mirror node  phi_mirror = 2*phi_BC - phi_P.
            aE = De_u + max(-Fe, 0);
            aW = Dw_u + max( Fw, 0);
            aN = Dn_u + max(-Fn, 0);
            aS = Ds_u + max( Fs, 0);
    
    
            % ── Correção diferida S_DC ─────────────────────────
    
            % Face leste
            if Fe >= 0
                S_east = -(Fe/2)*(uP - uW);
            else
                S_east = -(Fe/2)*(uE - uEE);
            end
    
            % Face oeste
            if Fw >= 0
                S_west =  (Fw/2)*(uW - uWW);
            else
                S_west =  (Fw/2)*(uP - uE);
            end
    
            % Face norte
            if Fn >= 0
                S_north = -(Fn/2)*(uP - uS);
            else
                S_north = -(Fn/2)*(uN - uNN);
            end
    
            % Face sul
            if Fs >= 0
                S_south =  (Fs/2)*(uS - uSS);
            else
                S_south =  (Fs/2)*(uP - uN);
            end
    
            Su_dc = S_east + S_west + S_north + S_south;

            a = [aE,aW,aN,aS];

        case "QUICK"
            aE = De_u + max(-Fe, 0);
            aW = Dw_u + max( Fw, 0);
            aN = Dn_u + max(-Fn, 0);
            aS = Ds_u + max( Fs, 0);
        
        
            % S_DC — face leste
            if Fe >= 0
                S_east = -(Fe/8)*(3*uE - 2*uP - uW);
            else
                S_east = -(Fe/8)*(3*uP - 2*uE - uEE);
            end
        
            % S_DC — face oeste
            if Fw >= 0
                S_west =  (Fw/8)*(3*uP - 2*uW - uWW);
            else
                S_west =  (Fw/8)*(3*uW - 2*uP - uE);
            end
        
            % S_DC — face norte
            if Fn >= 0
                S_north = -(Fn/8)*(3*uN - 2*uP - uS);
            else
                S_north = -(Fn/8)*(3*uP - 2*uN - uNN);
            end
        
            % S_DC — face sul
            if Fs >= 0
                S_south =  (Fs/8)*(3*uP - 2*uS - uSS);
            else
                S_south =  (Fs/8)*(3*uS - 2*uP - uN);
            end
        
            Su_dc = S_east + S_west + S_north + S_south;
            a = [aE,aW,aN,aS];

        case "MUSCL"
            aE = De_u + max(-Fe, 0);
            aW = Dw_u + max( Fw, 0);
            aN = Dn_u + max(-Fn, 0);
            aS = Ds_u + max( Fs, 0);
        
        
            % ── Limitador de van Albada ────────────────────────────────────
            % psi(r) = (r^2 + r) / (r^2 + 1)
            % kappa = 1/3  →  S_DC = -(F/6)*psi*[(uC-uU) + 2*(uA-uC)]
            vanalbada = @(r) (r > 0) .* (r.^2 + r) ./ (r.^2 + 1);  % psi = 0 em extremos (r <= 0)
        
            % ── Face leste ────────────────────────────────────────────────
            if Fe >= 0   % C=P, A=E, U=W
                d = uE - uP;
                if abs(d) < 1e-14, r = 0; else, r = (uP - uW) / d; end
                psi = vanalbada(r);
                S_east = -(Fe/6) * psi * ((uP - uW) + 2*(uE - uP));
            else         % C=E, A=P, U=EE
                d = uP - uE;
                if abs(d) < 1e-14, r = 0; else, r = (uE - uEE) / d; end
                psi = vanalbada(r);
                S_east = -(Fe/6) * psi * ((uE - uEE) + 2*(uP - uE));
            end
        
            % ── Face oeste ────────────────────────────────────────────────
            if Fw >= 0   % C=W, A=P, U=WW
                d = uP - uW;
                if abs(d) < 1e-14, r = 0; else, r = (uW - uWW) / d; end
                psi = vanalbada(r);
                S_west = (Fw/6) * psi * ((uW - uWW) + 2*(uP - uW));
            else         % C=P, A=W, U=E
                d = uW - uP;
                if abs(d) < 1e-14, r = 0; else, r = (uP - uE) / d; end
                psi = vanalbada(r);
                S_west = (Fw/6) * psi * ((uP - uE) + 2*(uW - uP));
            end
        
            % ── Face norte ────────────────────────────────────────────────
            if Fn >= 0   % C=P, A=N, U=S
                d = uN - uP;
                if abs(d) < 1e-14, r = 0; else, r = (uP - uS) / d; end
                psi = vanalbada(r);
                S_north = -(Fn/6) * psi * ((uP - uS) + 2*(uN - uP));
            else         % C=N, A=P, U=NN
                d = uP - uN;
                if abs(d) < 1e-14, r = 0; else, r = (uN - uNN) / d; end
                psi = vanalbada(r);
                S_north = -(Fn/6) * psi * ((uN - uNN) + 2*(uP - uN));
            end
        
            % ── Face sul ──────────────────────────────────────────────────
            if Fs >= 0   % C=S, A=P, U=SS
                d = uP - uS;
                if abs(d) < 1e-14, r = 0; else, r = (uS - uSS) / d; end
                psi = vanalbada(r);
                S_south = (Fs/6) * psi * ((uS - uSS) + 2*(uP - uS));
            else         % C=P, A=S, U=N
                d = uS - uP;
                if abs(d) < 1e-14, r = 0; else, r = (uP - uN) / d; end
                psi = vanalbada(r);
                S_south = (Fs/6) * psi * ((uP - uN) + 2*(uS - uP));
            end
        
            Su_dc = S_east + S_west + S_north + S_south;
            a = [aE,aW,aN,aS];
            
        case "HRIC"
            aE = De_u + max(-Fe, 0);
            aW = Dw_u + max( Fw, 0);
            aN = Dn_u + max(-Fn, 0);
            aS = Ds_u + max( Fs, 0);
        
        
            % ── Gradiente local para blending com ângulo ──────────────────
            % Gradiente de u em x e y no nó P (diferenças centrais com u_prev)
            grad_x = (uE - uW) / (2*dx);
            grad_y = (uN - uS) / (2*dy);
            grad_mag = sqrt(grad_x^2 + grad_y^2) + 1e-14;
        
            % Face leste/oeste: 
            cos2_x = (grad_x / grad_mag)^2;
            % Face norte/sul: 
            cos2_y = (grad_y / grad_mag)^2;
                
            % ── Face leste ────────────────────────────────────────────────
            if Fe >= 0   % C=P, A=E, U=W
                dDE = uE - uW;
                if abs(dDE) < 1e-14
                    S_east = 0;
                else
                    phiC = (uP - uW) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_x + phif_cd * (1 - cos2_x);
                    S_east = -Fe * (phif - phiC) * dDE;
                end
            else         % C=E, A=P, U=EE
                dDE = uP - uEE;
                if abs(dDE) < 1e-14
                    S_east = 0;
                else
                    phiC = (uE - uEE) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_x + phif_cd * (1 - cos2_x);
                    S_east = -Fe * (phif - phiC) * dDE;
                end
            end
        
            % ── Face oeste ────────────────────────────────────────────────
            if Fw >= 0   % C=W, A=P, U=WW
                dDE = uP - uWW;
                if abs(dDE) < 1e-14
                    S_west = 0;
                else
                    phiC = (uW - uWW) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_x + phif_cd * (1 - cos2_x);
                    S_west = Fw * (phif - phiC) * dDE;
                end
            else         % C=P, A=W, U=E
                dDE = uW - uE;
                if abs(dDE) < 1e-14
                    S_west = 0;
                else
                    phiC = (uP - uE) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_x + phif_cd * (1 - cos2_x);
                    S_west = Fw * (phif - phiC) * dDE;
                end
            end
        
            % ── Face norte ────────────────────────────────────────────────
            if Fn >= 0   % C=P, A=N, U=S
                dDE = uN - uS;
                if abs(dDE) < 1e-14
                    S_north = 0;
                else
                    phiC = (uP - uS) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_y + phif_cd * (1 - cos2_y);
                    S_north = -Fn * (phif - phiC) * dDE;
                end
            else         % C=N, A=P, U=NN
                dDE = uP - uNN;
                if abs(dDE) < 1e-14
                    S_north = 0;
                else
                    phiC = (uN - uNN) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_y + phif_cd * (1 - cos2_y);
                    S_north = -Fn * (phif - phiC) * dDE;
                end
            end
        
            % ── Face sul ──────────────────────────────────────────────────
            if Fs >= 0   % C=S, A=P, U=SS
                dDE = uP - uSS;
                if abs(dDE) < 1e-14
                    S_south = 0;
                else
                    phiC = (uS - uSS) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_y + phif_cd * (1 - cos2_y);
                    S_south = Fs * (phif - phiC) * dDE;
                end
            else         % C=P, A=S, U=N
                dDE = uS - uN;
                if abs(dDE) < 1e-14
                    S_south = 0;
                else
                    phiC = (uP - uN) / dDE;
                    if     phiC >= 0   && phiC < 0.5, phif_star = 2*phiC;
                    elseif phiC >= 0.5 && phiC < 1,   phif_star = 1;
                    else,                               phif_star = phiC;
                    end
                    phif_cd = (1 + phiC) / 2;
                    phif = phif_star * cos2_y + phif_cd * (1 - cos2_y);
                    S_south = Fs * (phif - phiC) * dDE;
                end
            end
        
            Su_dc = S_east + S_west + S_north + S_south;
            a = [aE,aW,aN,aS];
    end

    aP = aE + aW + aN + aS + (Fe - Fw + Fn - Fs);
end

function [n1, n2] = vizinhos(f, mask, a, b, c, d, uP, tipo, w, tipo_mask)
% Valores no 1o (a,b) e 2o (c,d) nó vizinhos, com contorno de valor w:
%   "no"    : contorno sobre um nó da malha deslocada
%   "ghost" : parede a meia célula (nó fantasma espelhado)
%   "saida" : gradiente nulo
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
