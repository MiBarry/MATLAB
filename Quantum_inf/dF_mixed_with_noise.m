function [dF, d] = dF_mixed_with_noise(x, po, t, n)

if r == 1
    fprintf('model: r = %d, actual: r = %d', r, rank(po))
    [rad, tet, phi] = return_r_tet_phi_by_po_matrix(po);
    po = build_po_matrix(1, tet, phi);
end

[m, s] = size(x); 

lambda = lambda_for_protocol_noise(x, m, po, T, T1);

C = purification_procedure2(po);
c = [real(C); imag(C)];

Lambda2 = Lambda_for_mixed_noise3(x,r,T,T1);

H = zeros(2*r*s,2*r*s);

Lambda = zeros(2*r*s,2*r*s,m);
for k = 1:m
    Lambda(:,:,k) = [[real(Lambda2(:,:,k)) -imag(Lambda2(:,:,k))]
                     [imag(Lambda2(:,:,k)) real(Lambda2(:,:,k))]];
end

for k=1:m
	if( lambda(k) > 10^(-30))
        H = H + t(k)/lambda(k)*(Lambda(:,:,k)*c)*(Lambda(:,:,k)*c)';
    end
end
H = 2 * H;
if abs(c'*H*c - 2*n) > 1e-3
error("c'*H*c ~= 2n")
end

S = eig(H);
S = sort(S);
nu = (2*s-r)*r;
S = S(length(S)-nu+1:end-1);
d = 1./(2.*S);
d = real(d);

dF = sum(d);
if(isinf(dF))
    dF
end
