function lam_max = calc_lamax(cpf, Yf, Imax, V0, dV, Y, genlist, Qmax, Qmin, limitQQ, limitIQ)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
step_max = length(cpf.cpf.steps); %number of steps of cpf

%initate the maximum step whose the current satisfies the respected limits
step_max_I = step_max;
step_max_V = step_max;
step_max_Q = step_max;

checkImax = maxI_feas_Q(Yf, cpf.cpf.V(:,step_max), Imax);
checkVmax = maxV_feas_Q(cpf.cpf.V(:,step_max), V0, dV);
checkQmax = maxQ_feas_Q(Y, cpf.cpf.V(:,step_max), genlist, Qmax , Qmin);

violateQ = 0;%if any constraint is violated

if ~min([checkImax; checkVmax; checkQmax])
    violateQ = 1;
    if ~checkImax
        if limitIQ
            step_max_I = bisection_Imax(cpf, Yf, Imax, step_max);
        else
            step_max_I = 1e+100;
        end
    end

    if ~checkVmax
        step_max_V = bisection_Vmax(cpf, step_max, V0, dV);
    end

    if ~checkQmax
        if limitQQ
            step_max_Q = bisection_Qmax(cpf, step_max, Y, genlist, Qmax, Qmin);
        else
            step_max_Q = 1e+100;
        end
    end     
end

if violateQ
  [~, violation] = min([step_max_I; step_max_V; step_max_Q ]) ;          
 lam_max = find_lamax(violation, cpf, Yf, Imax, V0, dV, Y, genlist, Qmax, Qmin, step_max_I, step_max_V, step_max_Q);
else
 lam_max = cpf.cpf.max_lam;
end

end

