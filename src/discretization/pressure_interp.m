function src_p = pressure_interp(direction, interp_pressao, p, dy, dx, i, j)
%PRESSURE_INTERP  Termo fonte de pressão no nó de velocidade (i,j).
%   direction : "u" ou "v"; interp_pressao : hoje só "Linear".
%   src_p = (p_P - p_E)*dy para u e (p_P - p_N)*dx para v.
    switch interp_pressao
        case "Linear"
            if direction == "u"
                src_p = (p(i,j) - p(i+1,j)) * dy;
            elseif direction == "v"
                src_p = (p(i,j) - p(i,j+1)) * dx;
            end
    end
end