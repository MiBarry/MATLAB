function [phi, tet] = return_phi_tet_2(psi)

phi = angle(psi(2))-angle(psi(1));
tet = 2*acos(abs(psi(1)));
