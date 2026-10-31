function [P_param, u_analitica] = analytical_Couette_Pouiseuille(physical_param)
    Ly = physical_param(2); mu = physical_param(3); 
    U = physical_param(5); G = physical_param(6);

    P_param = G * Ly^2 / (2 * mu * U); % parâmetro de Couette-Poiseuille
    % Perfil analítico em qualquer y (vetor coluna)
    u_analitica = @(y) U*(y/Ly) + (G*Ly/(2*mu)).*y.*(1 - y/Ly);
end