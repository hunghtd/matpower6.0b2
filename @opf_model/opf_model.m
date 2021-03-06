function om = opf_model(mpc)
%OPF_MODEL  Constructor for OPF model class.
%   OM = OPF_MODEL(MPC)
%
%   This class implements the OPF model object used to encapsulate
%   a given OPF problem formulation. It allows for access to optimization
%   variables, constraints and costs in named blocks, keeping track of the
%   ordering and indexing of the blocks as variables, constraints and costs
%   are added to the problem.
%
%   This class is a sub-class of OPT_MODEL and simply adds the 'mpc'
%   field for storing the MATPOWER case struct used to build the object
%   along with the get_mpc() method.
%
%   The following is the structure of the data in the OPF model object.
%   Each field of .idx or .data is a struct whose field names are the names
%   of the corresponding blocks of vars, constraints or costs (found in
%   order in the corresponding .order field). The description next to these
%   fields gives the meaning of the value for each named sub-field.
%   E.g. om.var.data.v0.Pg contains a vector of initial values for the 'Pg'
%   block of variables.
%
%   om
%       .opt_model  - the corresponding OPT_MODEL object
%       .mpc        - MATPOWER case struct used to create this model object
%           .baseMVA
%           .bus
%           .branch
%           .gen
%           .gencost
%           .A  (if present, must have l, u)
%           .l
%           .u
%           .N  (if present, must have fparm, H, Cw)
%           .fparm
%           .H
%           .Cw
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin == 0
    es = struct();
    s = struct('mpc', es);
    om = opt_model;
    om = class(s, 'opf_model', om);
else
    if isa(mpc,'opf_model') 
        om = mpc;
    else
        if isfield(mpc, 'om')   %% avoid nesting
            s = struct('mpc', rmfield(mpc, 'om'));
        else
            s = struct('mpc', mpc);
        end
        om = opt_model;
        om = class(s, 'opf_model', om);
    end
end
