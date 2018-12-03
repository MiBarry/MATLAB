function dF = dF_phase_flip(X, psi, prob, n)
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

rx = zeros(N,1);
ry = zeros(N,1);
rz = zeros(N,1);
for k=1:N
    rx(k) = sin(theta(k))*cos(phi(k));
    ry(k) = sin(theta(k))*sin(phi(k));
    rz(k) = cos(theta(k));
end

H = zeros(2*N0,2*N0);
% Lambda = zeros(2, 2);
for k=1:N
    buffer = [1+rz(k), (rx(k)+1i*ry(k))*(1-2*prob);
             (rx(k)-1i*ry(k))*(1-2*prob), 1-rz(k)];
    buffer = buffer*1/N;
%     Lambda = Lambda + buffer;
    buffer2 = [real(buffer), -imag(buffer); imag(buffer), real(buffer)];
    if( abs(c'*buffer2*c) > 10^(-30)) 
         H = H + (buffer2*c*c'*buffer2)/(c'*buffer2*c); 
    end
end
% Lambda
H = 2*n * H;
S = eig(H);
S = sort(S);
d = zeros(length(S)-2, 1);
for k=1:length(d)
   d(k) = 1 / (2 * S(k + 1)); 
end
dF = sum(d);
