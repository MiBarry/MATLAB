function dF = dF_phase_relaxation(X, psi, T, T2, n)
[N N0] = size(X);

c = [real(psi); imag(psi)];
phi = zeros(N,1);
theta = zeros(N,1);
s = zeros(N,1);
for k=1:N
    phi(k) = angle(X(k,2))-angle(X(k,1));
	theta(k) = 2*acos(abs(X(k,1)));
    s(k) = (angle(X(k,1))+angle(X(k,2)))/2;
end

H = zeros(2*N0,2*N0);
% Lambda = zeros(2, 2);
for k=1:N
    buffer = [1+cos(theta(k)),sin(theta(k))*exp(1i*phi(k))*exp(-T/T2);
             sin(theta(k))*exp(-1i*phi(k))*exp(-T/T2),1-cos(theta(k))];
    buffer = buffer*1/N;
%     Lambda = Lambda + buffer;
    buffer2 = [real(buffer), -imag(buffer); imag(buffer), real(buffer)];
    if( abs(c'*buffer2*c) > 10^(-30)) 
         H = H + (buffer2*c*c'*buffer2)/(c'*buffer2*c); 
    end
end
% Lambda
H = 2*4000 * H;

S = eig(H);
S = sort(S);
d = zeros(length(S)-2, 1);
for k=1:length(d)
   d(k) = 1 / (2 * S(k + 1)); 
end

dF = sum(d);
