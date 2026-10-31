function src_p = pressure_interp(direction, interp_pressao, p, dy, dx, i, j)
    switch interp_pressao
        case "Linear"
            if direction == "u"
                src_p = (p(i,j) - p(i+1,j)) * dy;
            elseif direction == "v"
                src_p = (p(i,j) - p(i,j+1)) * dx;
            end
    end
end