function [dF, d] = dF_amplitude_relaxationYYY(X, psi, ver, time)
[N N0] = size(X);
c = [real(psi); imag(psi)];
H = zeros(2*N0,2*N0);
for k=1:N
    Lambda = X(k,:)'*X(k,:);
    Lambda2 = [real(Lambda), -imag(Lambda); imag(Lambda), real(Lambda)];
    if( abs(c'*Lambda2*c) > 10^(-30)) 
         H = H + time/ver(k)*(Lambda2*c*c'*Lambda2); 
    end
end
H = 2 * H;
S = eig(H);
S = sort(S);
d = zeros(length(S)-2, 1);
for k=1:length(d)
   d(k) = 1 / (2 * S(k + 1)); 
end
dF = sum(d);
