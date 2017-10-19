function lam_max = lamax_Imax(cpf, Yf, Imax, step_max_I)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here    
if step_max_I > 0
    [checknextstep, vline] = maxI_feas_Q( Yf, cpf.cpf.V(:,step_max_I+1), Imax);

    if checknextstep
       error('the critial step is wrong. Check bisection iteration')
    else
        cacl_iline = @(Yf, V) abs(Yf(vline, :)*V);

        Iv_flr = cacl_iline(Yf, cpf.cpf.V(:,step_max_I));
        Iv_ceil = cacl_iline(Yf, cpf.cpf.V(:,step_max_I+1));
        fac = (Imax(vline) - Iv_flr) / (Iv_ceil - Iv_flr);
    end

    lam_max = fac * (cpf.cpf.lam(step_max_I+1) - cpf.cpf.lam(step_max_I)) + cpf.cpf.lam(step_max_I);  
else
    lam_max = 0;
end

end

