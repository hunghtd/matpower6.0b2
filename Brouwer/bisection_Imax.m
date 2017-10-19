function p = bisection_Imax(cpf, Yf, Imax, step)
a = 1;
b = step;

if maxI_feas_Q(Yf, cpf.cpf.V(:, a), Imax) -  maxI_feas_Q(Yf, cpf.cpf.V(:, b), Imax) == 0 
    disp('The base power flow solution is not feasible')
    p = 0; %take the original point
else
    p = floor((a + b)/2);
    err = b-a;
    while err >1
        if (maxI_feas_Q(Yf, cpf.cpf.V(:, a), Imax) -  maxI_feas_Q(Yf, cpf.cpf.V(:, p), Imax) > 0)
             b = p;
        else
             a = p;          
        end
        p = floor((a + b)/2); 
        err = b-a;
    end  
end

end