function ef = cpf_qlim_event(cb_data, cx)
%CPF_QLIM_EVENT  Event function to detect the generator reactive power limits
%   EF = CPF_QLIM_EVENT(CB_DATA, CX)
%
%   CPF event function to detect a generator reactive power limits,
%   i.e. Ip  >= Imax.
%
%   Inputs:
%       CB_DATA : struct of data for callback functions
%       CX : struct containing info about current point (continuation soln)
%
%   Outputs:
%       EF : event function value

%   MATPOWER
%   Copyright (c) 2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% event function value is 1 np x 1 vector equal to:
%%      [ Ip - Imax ]

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get updated MPC
d = cb_data;
mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
    d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, cx.V, cx.lam, d.mpopt);
%% calculate branch currents
V = mpc.bus(:, VM) .* exp(1i*mpc.bus(:, VA)/180*pi);
[~,Yf,~] = makeYbus(mpc);
I = Yf*V;
I = abs(I);
%% compute current branch violations for on-line lines
np = size(mpc.branch, 1);
v_Imax = NaN(np, 1);
on = find(mpc.branch(:, BR_STATUS) > 0);  %% which lines are on?
v_Imax(on) = I(on) - mpc.branch(on, end); %% mpc.branch(on, end) is the Imax

%% assemble event function value
ef = v_Imax;
