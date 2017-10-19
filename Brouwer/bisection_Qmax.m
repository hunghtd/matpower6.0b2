function p = bisection_Qmax(cpf, step, Y, genlist, Qmax, Qmin)
a = 1;
b = step;

if maxQ_feas_Q(Y, cpf.cpf.V(:, a), genlist, Qmax, Qmin) -  maxQ_feas_Q(Y, cpf.cpf.V(:, b), genlist, Qmax, Qmin) == 0 
    disp('The base power flow solution is not feasible')
    p=0;
else
    p = floor((a + b)/2);
    err = b-a;
    while err >1
        if (maxQ_feas_Q(Y, cpf.cpf.V(:, a), genlist, Qmax, Qmin) -  maxQ_feas_Q(Y, cpf.cpf.V(:, p), genlist, Qmax, Qmin) > 0)
             b = p;
        else
             a = p;          
        end
        p = floor((a + b)/2); 
        err = b-a;
    end  
end

end