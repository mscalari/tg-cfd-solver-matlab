function safe_sgtitle(txt)
%SAFE_SGTITLE  Chama sgtitle (MATLAB R2018b+) se disponível; caso
%   contrário, cai para um título de figura via annotation de texto no
%   topo, evitando que o script quebre em versões mais antigas do
%   MATLAB ou no Octave (onde sgtitle pode não existir).

    try
        sgtitle(txt);
    catch
        try
            annotation('textbox', [0 0.95 1 0.05], 'String', txt, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'FontSize', 11);
        catch
            % último recurso: não falha o script por causa do título
        end
    end
end
