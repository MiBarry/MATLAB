clear
clc
close all
addpath('Protocols');
state0 = rand;
state1 = rand;
prob1= randi(10);
prob2= randi(10);
prob3= randi(10);
prob4= randi(10);
if (mod(prob1,2)==0)
    state0 = state0 * (-1);
end
if (mod(prob2,2)==0)
    state1 = state1 * (-1);
end
if (mod(prob3,2)==0)
    state0 = state0 * 1i;
end
if (mod(prob4,2)==0)
    state1 = state1 * 1i;
end
C=[state0;state1];
%C = [0; 1];
C=C/norm(C);
time = 1000;
matrix = X_Cube
matrix = matrix * rotation_operator(pi/32,pi/32,pi/32)
ver=zeros(6,1);
for s=1:6
ver(s,1)=(abs(matrix(s,:)*C))^2;
end
graph=zeros(1000,1);
for z=1:1000
R = poissrnd(ver*time);
I=zeros(2);
J=zeros(2);
for s=1:2
    for n=1:2
        res=0;
        for m=1:6
            res=res+conj(matrix(m,s))*matrix(m,n)*time;
        end
        I(s,n)=res;
    end
end

C0 = C;
Ck1=C0;
alf=0.5;
ver1=zeros(6,1);
for t=1:100
    Ck=C0;
    for q=1:6
    ver1(q,1)=(abs(matrix(q,:)*C0))^2;
    end
       for s=1:2
            for n=1:2
                res2=0;
                for m=1:6
                    res2=res2+conj(matrix(m,s))*matrix(m,n)*R(m)/ver1(m);
                end
                J(s,n)=res2;
            end
       end
    C0=(1-alf)*I^(-1)*J*Ck+alf*Ck;
    C0 = C0 / norm(C0); 
end
F=(abs(C'*C0)).^2;
graph(z,1)=1-F;
end
sre=sum(graph)/length(graph); 


fprintf("dF exp = %f\n", sre)
[dF, d] = dF_amplitude_relaxationYYY(matrix, C,ver, time);
fprintf("dF teor = %f\n", dF)
plot_hist_and_generalchi2pdf(graph, d, 1000)
