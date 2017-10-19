function psse2m(rawfile, mfilename)
% psse2m to convert a raw file in psse to m file in MATPOWER
% mfilename with extension .m
% example psse2m('118.raw', '118.m')
%   Detailed explanation goes here
mpc = psse2mpc(rawfile);
savecase(mfilename, mpc)
end