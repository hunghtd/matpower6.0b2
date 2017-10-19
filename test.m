addpath('/Users/hunghtd/Documents/MATLAB/matpower6.0b2');
              define_constants;
              mpc = loadcase('case5');
              mpc = ext2int(mpc);

              % Until we implement the shunts
              %mpc.branch(:, BR_B) = 0.0;
              %mpc.bus(:, GS) = 0.0;
              %mpc.bus(:, BS) = 0.0;

              %some useful functions
              spd = @(vec) sparse(1:numel(vec), 1:numel(vec), vec);
              spz = @(n, m) sparse([],[],[],n,m);

              %bus type and size
              slack = find(mpc.bus(:, BUS_TYPE) == 3);
              pq = find(mpc.bus(:, BUS_TYPE) == 1);
              pv = find(mpc.bus(:, BUS_TYPE) == 2);
              nslack = size(slack, 1);
              npq = size(pq, 1);
              npv = size(pv, 1);
              pvpq = [pv; pq];

              Y = makeYbus(mpc);
              nb = size(Y,1);
              np = numel(mpc.branch(:,F_BUS));

              sol = runpf(mpc);
              sol = ext2int(sol);
              S_s = makeSbus(sol.baseMVA, sol.bus, sol.gen);

              J_s = makeJac(sol);

              branch = sol.branch;
              nl = size(branch, 1);       %% number of lines

              % Borrowed from makeYbus
              from = branch(:, F_BUS);                           %% list of from buses
              to = branch(:, T_BUS);                           %% list of to buses
              Cf = sparse(1:nl, branch(:, F_BUS), ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
              Ct = sparse(1:nl, branch(:, T_BUS), ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
              stat = branch(:, BR_STATUS);                    %% ones at in-service branches
              Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance

              Ic = Ct' + Cf';
              Ec = Ct' - Cf';

              [Ybus, Yf, Yt] = makeYbus(sol.baseMVA, sol.bus, sol.branch);

              stat = branch(:, BR_STATUS);                    %% ones at in-service branches
              Ybr = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance
              B_br = spd(imag(Ybr));
              G_br = spd(real(Ybr));
              G_bus = spd(Ic*G_br*ones(nl,1));
              B_bus = spd(Ic*B_br*ones(nl,1));

              hSm = [
                   G_bus(pv,:) -Ic(pv,:)*G_br Ec(pv,:)*B_br ;
                   G_bus(pq,:) -Ic(pq,:)*G_br Ec(pq,:)*B_br ;
                  -B_bus(pq,:)  Ic(pq,:)*B_br Ec(pq,:)*G_br
                  ];

              Vm_s = sol.bus(:,VM);
              Va_s = sol.bus(:,VA)/180*pi;

              dVa_s = -Ec'*Va_s;
              W_s = (Ct*Vm_s).*(Cf*Vm_s);

              hphi_s = [Vm_s.*Vm_s; W_s.*cos(dVa_s); W_s.*sin(dVa_s)];

              sc_s = makeSbus(sol.baseMVA, sol.bus, sol.gen);
              s_s = [real(sc_s(pv)); real(sc_s(pq)); imag(sc_s(pq))];
              %assert(norm(hSm*hphi_s - s_s, inf) < 1e-6, 'Base solution not recovered')

              % Switch from cos, sin to cs, sn:
              sp_cos_s = spd(cos(dVa_s));
              sp_sin_s = spd(sin(dVa_s));

              cs_s = ones(size(dVa_s,1),1);
              sn_s = zeros(size(dVa_s,1),1);

              Om = [spd(ones(nb,1)) spz(nb,nl) spz(nb,nl) ;
                    spz(nl, nb)     sp_cos_s  -sp_sin_s;
                    spz(nl, nb)     sp_sin_s   sp_cos_s];
              Sm = hSm*Om;
              phi_s = [Vm_s.*Vm_s; W_s; 0*W_s];
              %assert(norm(Sm*phi_s - s_s, inf) < 1e-6, 'Base solution not recovered')

              % create matrix A
              id = speye(nb);
              E = sparse([(1:np)';(1:np)'], [from; to], [ones(np,1); -ones(np,1)], np, nb)
              A = [0*speye(nb, npv + npq) id(:,pq);
                    E(:,pv) E(:,pq) 0*speye(np, npq)];
                
                C = -A/J_s*Sm