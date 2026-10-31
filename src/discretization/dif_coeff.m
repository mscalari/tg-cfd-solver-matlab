function D = dif_coeff(mu,dx,dy)
    %DIF_COEFF  Condutâncias difusivas de u e v em malha uniforme com mu constante.
    %  Devolve D = [De_u, Dw_u, Dn_u, Ds_u, De_v, Dw_v, Dn_v, Ds_v].
    %    De = mu * dy/dx   (difusão u na dir. x entre P e E)
    %    Dw = mu * dy/dx
    %    Dn = mu * dx/dy   (difusão u na dir. y entre P e N)
    %    Ds = mu * dx/dy
    
    %  Para v, as condutâncias são análogas mas com dx<->dy trocados.
    
    De_u = mu * dy / dx;
    Dw_u = mu * dy / dx;
    Dn_u = mu * dx / dy;
    Ds_u = mu * dx / dy;
    
    De_v = mu * dy / dx;
    Dw_v = mu * dy / dx;
    Dn_v = mu * dx / dy;
    Ds_v = mu * dx / dy;
        
    D = [De_u, Dw_u, Dn_u,Ds_u,De_v,Dw_v,Dn_v,Ds_v];
end