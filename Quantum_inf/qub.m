clear
clc
close all
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
C=[state0;state1];
C=C/norm(C);
time = 1000;
matrix = [1 0; 0 1; 1/sqrt(2)  1/sqrt(2); 1/sqrt(2)  (-1)/sqrt(2); 1/sqrt(2)  (-sqrt(-1))/sqrt(2); 1/sqrt(2)  (sqrt(-1))/sqrt(2)];
ampl=zeros(6,1);
ver=zeros(6,1);
for s=1:6
ampl(s,1)=matrix(s,:)*C;
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
C0=[1/sqrt(5);1/sqrt(5)*2*sqrt(-1)];
Ck1=C0;
alf=0.5;
ver1=zeros(6,1);
for t=1:100
    Ck=Ck1;
    for q=1:6
    ver1(q,1)=(abs(matrix(q,:)*Ck1))^2;
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
    Ck1=(1-alf)*I^(-1)*J*Ck+alf*Ck;
    Ck1 = Ck1 / norm(Ck1); 
end
F=(abs(C'*Ck1)).^2;
graph(z,1)=1-F;
end
sre=sum(graph)/length(graph); 
sre

%{
XX=0:0.1:3;
MU =0;
SIGMA =1;

DDD=[0.000187499988302079;9.37500029237963e-05]
F1=300-1000000*(DDD(1)*normcdf(XX,MU,SIGMA).*normcdf(XX,MU,SIGMA)+DDD(2)*normcdf(XX,MU,SIGMA).*normcdf(XX,MU,SIGMA));

hist(graph,30)
plot(XX,F1)
xlabel('�������� 1-F'); 
ylabel('����� �������������'); 
text(1,40,'��������� ���������:');
grid on;

%}
%%

[dF, d] = dF_amplitude_relaxationYYY(matrix, C, T, T1, n)
%plot_hist_and_generalchi2pdf(graph, d, 1000)
    dF_set = sort(graph);
    min_x = min(dF_set);
    max_x = max(dF_set);
    d_x = (max_x - min_x)/100;
    x_p = min_x:d_x:max_x;
    p = chi2pdf_general_bogdanov(x_p,d);
    figure
    h = histogram(dF_set,50);
    p = p * h.BinWidth * 4000;
    hold on
    plot(x_p, p, 'r', 'LineWidth', 2)
    hold off
