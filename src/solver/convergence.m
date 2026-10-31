function [res, res_hist, flag] = convergence(vetor, iter, res_hist, numerical_param)
%CONVERGENCE  Registra o resíduo da iteração e decide se o laço deve parar.
%
%   ENTRADAS
%     vetor    : campo de resíduo de massa
%     iter     : iteração atual
%     res_hist : histórico; é cortado em iter ao convergir ou virar NaN
%     numerical_param : ver build_case (usa tol e max_iter)
%
%   SAÍDAS
%     res  : norma de Frobenius do resíduo
%     flag : true se convergiu, atingiu max_iter ou divergiu para NaN
    tol = numerical_param{3};
    max_iter = numerical_param{4};

    flag = false;
    res = norm(vetor, 'fro');
    res_hist(iter) = res;

    if mod(iter,100)==0 || iter==1
        fprintf(' %4d    %.4e\n', iter, res);
    end

    if res < tol
        fprintf(' %4d    %.4e  <-- CONVERGIU\n', iter, res);
        res_hist = res_hist(1:iter);
        flag = true;
    end

    if iter == max_iter 
        flag = true;
    end

    if isnan(res)
        flag = true;
        res_hist = res_hist(1:iter);
    end
end