clear
clc
T=0;
T2=1;
T1=1;
n = 4000;

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
C = [exp(2i*pi/2*rand()); exp(2i*pi/2*rand());exp(2i*pi/2*rand());exp(2i*pi/2*rand())];
C=C/norm(C);
time = 1000;
matrix = [1 0; 0 1; 1/sqrt(2)  1/sqrt(2); 1/sqrt(2)  (-1)/sqrt(2); 1/sqrt(2)  (-sqrt(-1))/sqrt(2); 1/sqrt(2)  (sqrt(-1))/sqrt(2)];
matrix = kron(matrix, matrix)
ampl=zeros(36,1);
ver=zeros(36,1);
for s=1:36
ampl(s,1)=matrix(s,:)*C;
ver(s,1)=(abs(matrix(s,:)*C))^2;
end
graph=zeros(1000,1);
for z=1:1000
R = poissrnd(ver*time);
I=zeros(4);
J=zeros(4);
for s=1:4
    for n=1:4
        res=0;
        for m=1:36
            res=res+conj(matrix(m,s))*matrix(m,n)*time;
        end
        I(s,n)=res;
    end
end
C0=[1/sqrt(5);1/sqrt(5)*2*sqrt(-1)];
C0 = C;
Ck1=C0;
alf=0.5;
ver1=zeros(36,1);
for t=1:100
    Ck=C0;
    for q=1:36
    ver1(q,1)=(abs(matrix(q,:)*C0))^2;
    end
       for s=1:4
            for n=1:4
                res2=0;
                for m=1:36
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
sre
hist(graph,30)
xlabel('�������� 1-F'); 
ylabel('����� �������������'); 
text(1,40,'��������� ���������:');

DDD=[0.000187499988302079;9.37500029237963e-05]
XX=0:0.1:3;
MU =0;
SIGMA =1;
F1=300-1000000*(DDD(1)*normcdf(XX,MU,SIGMA).*normcdf(XX,MU,SIGMA)+DDD(2)*normcdf(XX,MU,SIGMA).*normcdf(XX,MU,SIGMA));
figure
plot(XX,F1)
grid on;
