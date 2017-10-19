function p = bisection_Vmax(cpf, step, V0, dV)
a = 1;
b = step;

if maxV_feas_Q(cpf.cpf.V(:, a), V0, dV) ==  maxV_feas_Q(cpf.cpf.V(:, b), V0, dV)
    disp('The base power flow solution is not feasible')
    p = 0;
else
    p = floor((a + b)/2);
    err = b-a;
    while err >1
        if (maxV_feas_Q(cpf.cpf.V(:, a), V0, dV) -  maxV_feas_Q(cpf.cpf.V(:, p), V0, dV) > 0)
             b = p;
        else
             a = p;          
        end
        p = floor((a + b)/2);
        err = b-a;
    end  
end

end