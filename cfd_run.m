function cfd_run(alvo, varargin)
%CFD_RUN  Entrada única do projeto: simulação, comparação ou testes.
%
%   Pelo terminal do sistema:
%     matlab -batch "cfd_run LidDrivenCavity algoritmo SIMPLEC esquema quick Re 400 Nx 64 Ny 64"
%     matlab -batch "cfd_run compare estudos 135 fast 1"
%     matlab -batch "cfd_run test"
%
%   No MATLAB, as mesmas chamadas funcionam digitando-as direto, ou como função:
%     cfd_run("FlowOverStep", "algoritmo", "PISO", "Re", 50)
%
%   ENTRADAS
%     alvo : escoamento (CouettePoiseuille, LidDrivenCavity, FlowOverCavity,
%            FlowOverStep, ou um prefixo deles), "compare" ou "test"
%     Pares nome/valor de uma simulação:
%       algoritmo : SIMPLE | SIMPLEC | PISO (padrão PISO)
%       esquema   : upwind1 | upwind2 | central | powerlaw | quick | muscl |
%                   hric (padrão upwind1)
%       demais    : qualquer parâmetro aceito por build_case
%     Pares nome/valor de "compare":
%       estudos : um dígito por estudo, entre 1, 3, 4 e 5 (padrão 135)
%       fast    : 1 para malhas menores (padrão 0)
%
%   Cada simulação grava PNGs e resultado.mat em results/, numa pasta nova.

    raiz = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(raiz, 'src')), fullfile(raiz, 'tests'));

    if nargin < 1
        help cfd_run
        return
    end
    if mod(numel(varargin), 2) ~= 0
        error('cfd_run: os parâmetros devem vir em pares nome/valor.');
    end
    nomes   = cellfun(@char, varargin(1:2:end), 'UniformOutput', false);
    valores = varargin(2:2:end);

    switch lower(char(alvo))
        case 'test'
            test_solver
        case 'compare'
            comparar(nomes, valores, fullfile(raiz, 'results', 'comparacao_parte_1'));
        otherwise
            simular(alvo, nomes, valores, fullfile(raiz, 'results'));
    end
end

function simular(alvo, nomes, valores, results)
%SIMULAR  Roda um caso, imprime o resumo e salva figuras e campos.
    escoamento = string(validatestring(char(alvo), ...
        {'CouettePoiseuille', 'LidDrivenCavity', 'FlowOverCavity', 'FlowOverStep'}, 'cfd_run', 'escoamento'));

    [algoritmo, nomes, valores] = extrair('algoritmo', nomes, valores, 'PISO');
    algoritmo = string(validatestring(char(algoritmo), {'SIMPLE', 'SIMPLEC', 'PISO'}, 'cfd_run', 'algoritmo'));
    [apelido, nomes, valores] = extrair('esquema', nomes, valores, 'upwind1');

    % Apelidos evitam espaços e acentos na linha de comando
    apelidos = {'upwind1', 'upwind2', 'central', 'powerlaw', 'quick', 'muscl', 'hric'};
    esquemas = {'Upwind 1a ordem', 'Upwind 2a ordem', 'Diferenças Centrais', ...
                'Power Law', 'QUICK', 'MUSCL', 'HRIC'};
    apelido = validatestring(char(apelido), apelidos, 'cfd_run', 'esquema');
    esquema = esquemas{strcmp(apelido, apelidos)};

    args = cell(1, 2*numel(nomes));
    args(1:2:end) = nomes;
    args(2:2:end) = cellfun(@numero, nomes, valores, 'UniformOutput', false);
    c = build_case(escoamento, args{:});

    fprintf('\n %s | %s | %s | malha %dx%d | Re = %.1f\n', ...
        escoamento, algoritmo, esquema, c.Nx, c.Ny, c.Re);
    fprintf(' Iter      Residuo de massa\n--------------------------------\n');
    r = run_case(c, algoritmo, esquema);

    if r.converged, estado = 'convergiu'; else, estado = 'NAO convergiu'; end
    fprintf('\n %s: %d iteracoes, residuo %.4e, %.1f s\n', estado, r.n_iter, r.res_final, r.t_cpu);

    instante = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    pasta = fullfile(results, sprintf('%s_%s_%s_%s', escoamento, algoritmo, apelido, instante));
    mkdir(pasta);
    if any(isnan(r.u(:)))
        fprintf(' Solucao divergiu (NaN): figuras omitidas.\n');
    else
        plot_results(c, r, algoritmo, esquema, pasta);
    end
    save(fullfile(pasta, 'resultado.mat'), 'c', 'r', 'algoritmo', 'esquema');
    fprintf(' Resultados em: %s\n', pasta);
end

function comparar(nomes, valores, pasta)
%COMPARAR  Lê as opções de "compare" e dispara run_comparison.
    [estudos, nomes, valores] = extrair('estudos', nomes, valores, [1 3 5]);
    [fast, nomes, ~] = extrair('fast', nomes, valores, 0);
    if ~isempty(nomes)
        error('cfd_run: opção desconhecida de compare: %s.', strjoin(nomes, ', '));
    end
    % Texto vira um estudo por dígito: a vírgula encerraria a sintaxe de comando
    if ~isnumeric(estudos)
        texto = char(estudos);
        estudos = texto(isstrprop(texto, 'digit')) - '0';
    end
    run_comparison(estudos, numero('fast', fast) ~= 0, pasta);
end

function [valor, nomes, valores] = extrair(nome, nomes, valores, padrao)
%EXTRAIR  Retira o par nome/valor da lista; devolve padrao se ele não existir.
    k = find(strcmp(nome, nomes), 1);
    valor = padrao;
    if ~isempty(k)
        valor = valores{k};
        nomes(k) = [];
        valores(k) = [];
    end
end

function v = numero(nome, v)
%NUMERO  Converte para número um valor recebido como texto do terminal.
    if ~isnumeric(v)
        texto = v;
        v = str2double(texto);
        if isnan(v)
            error('cfd_run: valor inválido para %s: "%s".', nome, texto);
        end
    end
end
