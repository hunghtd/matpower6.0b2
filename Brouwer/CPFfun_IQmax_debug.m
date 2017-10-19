clear all; close all;
testcase = 'case9'
busp = [7 9]
nump = 10
dV= 0.05
Angle = pi
rate = 2
stepsize = 0.01
adaptive = 0

define_constants;
mpc = loadcase(testcase);

for k = 1:numel(busp)
    busp(k) = find(mpc.bus(:,1)==busp(k));
end

mpc = ext2int(mpc);
% busp = mpc.order.bus.e2i(busp);%index new buses after arranging

%if the chosen buses are not pq buses, reselect those
pq = find((mpc.bus(:, BUS_TYPE) == 1) & (mpc.bus(:, PD) >0));
if sum(ismember(busp, pq)) ~= 2
    busp = pq(1:2);
end

SBASE = mpc.baseMVA;
genlist = mpc.gen(:, GEN_BUS);

[Y, Yf, ~] = makeYbus(mpc);

% base solution
pfbase = runpf(mpc);
Vpfbase = pfbase.bus(:,VM).*(exp(1i*pfbase.bus(:,VA)*pi/180));
s0 = Vpfbase.*conj(Y*Vpfbase);
V0 = abs(Vpfbase);

% calculate branch currents & Imax
Imax = calc_Imax(mpc, Yf, Vpfbase, rate);

% add qmax, qmin
[Qmax, Qmin] = calc_Qlim(mpc, s0, SBASE, genlist);

%% CPF
targetcase = mpc;
basecase = mpc;

mpopt = mpoption('cpf.stop_at', 'NOSE', 'cpf.step', stepsize, 'cpf.adapt_step', adaptive);

resR = zeros(nump,2); %save the loadability limits

%define # of quadrants to plot
% if nargin < 8
%     phi = linspace(0, pi/2, nump);
% elseif nargin == 8
    phi = linspace(0, Angle, nump);
% else
%     error('requires at most 8 optional inputs')
% end

% now run the CPF
% for k = nump-1:nump-1
       k = nump-1
        targetcase.bus(busp, PD) = (1+ [cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD);
        targetcase.bus(busp, QD) = (1+ [cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, QD);

        cpf = runcpf(basecase, targetcase, mpopt);
        
       % find the loadability level satisfying thermal limits, voltage limits, and reactive power limits
       lam_max = calc_lamax(cpf, Yf, Imax, V0, dV, Y, genlist, Qmax, Qmin)
%%
       step_max = length(cpf.cpf.steps); %number of steps of cpf

%initate the maximum step whose the current satisfies the respected limits
step_max_I = step_max;
step_max_V = step_max;
step_max_Q = step_max;

checkImax = maxI_feas_Q(Yf, cpf.cpf.V(:,step_max), Imax);
checkVmax = maxV_feas_Q(cpf.cpf.V(:,step_max), V0, dV);
checkQmax = maxQ_feas_Q(Y, cpf.cpf.V(:,step_max), genlist, Qmax , Qmin);

violate_Q = 0;%if any constraint is violated

if ~min([checkImax; checkVmax; checkQmax])
    violate_Q = 1;
    if ~checkImax
        step_max_I = bisection_Imax(cpf, Yf, Imax, step_max);
    end

    if ~checkVmax
        step_max_V = bisection_Vmax(cpf, step_max, V0, dV);
    end

    if ~checkQmax
        step_max_Q = bisection_Qmax(cpf, step_max, Y, genlist, Qmax, Qmin);
    end     
end

if violate_Q
  [~, violation] = min([step_max_I; step_max_V; step_max_Q ]) ;          
 lam_max = find_lamax(violation, cpf, Yf, Imax, V0, dV, Y, genlist, Qmax, Qmin, step_max_I, step_max_V, step_max_Q);
else
 lam_max = cpf.cpf.max_lam;
end
       
%%
 if step_max_V > 0
    [checknextstep, vbus] = maxV_feas_Q( cpf.cpf.V(:,step_max_V+1), V0, dV );

    if checknextstep
       error('the critial step is wrong. Check bisection iteration')
    else
        Vv_flr = abs(cpf.cpf.V(vbus,step_max_V));
        Vv_ceil = abs(cpf.cpf.V(vbus,step_max_V+1));

        if Vv_flr > Vv_ceil %increase power -> voltage decrease, lower voltage violation
           Vmin = V0(vbus)-dV;
           fac = (Vv_flr - Vmin) / (Vv_flr - Vv_ceil);
        else
             Vmax = V0(vbus)+dV;
             fac = (Vmax - Vv_flr) / (Vv_ceil - Vv_flr);
        end 
    end

    lammax = fac * (cpf.cpf.lam(step_max_V+1) - cpf.cpf.lam(step_max_V)) + cpf.cpf.lam(step_max_V);  
else
    lammax = 0;
end      
%%
       
       
       
       resR(k, :) = (1+lam_max*[cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD)/SBASE; 
%        resR(k, :) = (1+lam_max*[cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD)/SBASE; 
       busp
% end
plot(resR(:,1),resR(:,2))



%%
% define_constants;
% mpc = loadcase(testcase);
% 
% for k = 1:numel(busp)
%     busp(k) = find(mpc.bus(:,1)==busp(k));
% end
% 
% mpc = ext2int(mpc);
% SBASE = mpc.baseMVA;
% 
% [Y,Yf,~] = makeYbus(mpc);
% 
% % base solution
% pfbase = runpf(mpc);
% Vpfbase = pfbase.bus(:,VM).*(exp(1i*pfbase.bus(:,VA)*pi/180));
% s0 = Vpfbase.*conj(Y*Vpfbase);
% 
% % calculate branch currents & Imax
% I = Yf*Vpfbase;
% Imax = abs(I)*2;
% 
% on = find(mpc.branch(:, BR_STATUS) > 0);  %% which lines are on?
% 
% if length(on) == length(Imax)
%     mpc.branch(on, end+1) = Imax; %% append Imax values to the branch data
% else
%     disp('Something is wrong with the lines')
% end
% 
% %add qmax, qmin
% Qmax = mpc.gen(:, QMAX)/SBASE;
% Qmin = mpc.gen(:, QMIN)/SBASE;
% genlist = mpc.gen(:, GEN_BUS);
% 
% q = imag(s0);
% qmax_vgen = find(q(genlist) > Qmax);
% qmin_vgen = find(q(genlist) < Qmin);
% Qmax(qmax_vgen ) = 1e+7; %if reactive power maximum limits are violated, raise the limits
% Qmin(qmin_vgen ) = -1e+7; %if reactive power minimum limits are violated, lower the limits
% 
% 
% %if the chosen buses are not pq buses, reselect those
% pq = find((mpc.bus(:, BUS_TYPE) == 1) & (mpc.bus(:, PD) >0));
% if sum(ismember(busp, pq)) ~= 2
%     busp = pq(1:2);
% end
% 
% %% CPF
% targetcase = mpc;
% basecase = mpc;
% 
% mpopt = mpoption('cpf.stop_at', 'NOSE', 'cpf.step', 0.1) ;%, 'cpf.adapt_step', 1);
% 
% resR = zeros(nump,2); %save the loadability limits
% 
% %define # of quadrants to plot
% % if nargin < 5
%     phi = linspace(0, pi/2, nump);
% % elseif nargin == 5
% %     phi = linspace(0, Angle, nump);
% % else
% %     error('requires at most 4 optional inputs')
% % end
% 
% % now run the CPF
% for k = 1:nump
%         targetcase.bus(busp, PD) = (1+ [cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD);
%         targetcase.bus(busp, QD) = (1+ [cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, QD);
% 
%         cpf = runcpf(basecase, targetcase, mpopt);
%         
%        %% check thermal limits, voltage limits, and reactive power limits
%         step_max = length(cpf.cpf.steps); %number of steps of cpf
%         
%         %initate the maximum step whose the current satisfies the respected limits
%         step_max_I = step_max;
%         step_max_V = step_max;
%         step_max_Q = step_max;
%         
%         checkImax = maxI_feas_Q(Yf, cpf.cpf.V(:,step_max), Imax);
%         checkVmax = maxV_feas_Q(cpf.cpf.V(:,step_max), dV);
%         checkQmax = maxQ_feas_Q(Y, cpf.cpf.V(:,step_max), genlist, Qmax , Qmin);
%         
%         violate_Q = 0;
%         
%         if ~min([checkImax; checkVmax; checkQmax])
%             violate_Q = 1;
%             if ~checkImax
%                 step_max_I = bisection_Imax(cpf, Yf, Imax, step_max);
%             end
% 
%             if ~checkVmax
%                 step_max_V = bisection_Vmax(cpf, step_max, dV);
%             end
% 
%             if ~checkQmax
%                 step_max_Q = bisection_Qmax(cpf, step_max, Y, genlist, Qmax, Qmin);
%             end     
%         end
%         
%         if violate_Q
%           [~, violation] = min([step_max_I; step_max_V; step_max_Q ]) ;          
%          lam_max = find_lamax(violation, cpf, Yf, Imax, dV, Y, genlist, Qmax, Qmin, step_max_I, step_max_V, step_max_Q);
%         else
%          lam_max = cpf.cpf.max_lam;
%         end
%                  
% resR(k, :) = (1+lam_max*[cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD)/SBASE;
%    
% end
% close all
% plot(resR(:,1),resR(:,2))
