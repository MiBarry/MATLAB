clear
clc
vert=20;
tau=(1+sqrt(5))/2;
aZ=sqrt(2+3*tau);

dode1=[aZ, aZ, aZ];
dode2=[-aZ, aZ, aZ];
dode3=[aZ, -aZ, aZ];
dode4=[aZ, aZ, -aZ];
dode5=[-aZ, -aZ, aZ];
dode6=[-aZ, aZ, -aZ];
dode7=[aZ, -aZ, -aZ];
dode8=[-aZ, -aZ, -aZ];
dode9=[0, (2*(aZ^2)-1)/3, ((aZ^2)-2)/3];
dode10=[0, -(2*(aZ^2)-1)/3, ((aZ^2)-2)/3];
dode11=[0, (2*(aZ^2)-1)/3, -((aZ^2)-2)/3];
dode12=[0, -(2*(aZ^2)-1)/3, -((aZ^2)-2)/3];
dode13=[(2*(aZ^2)-1)/3, ((aZ^2)-2)/3, 0];
dode14=[-(2*(aZ^2)-1)/3, ((aZ^2)-2)/3, 0];
dode15=[(2*(aZ^2)-1)/3, -((aZ^2)-2)/3, 0];
dode16=[-(2*(aZ^2)-1)/3, -((aZ^2)-2)/3, 0];
dode17=[((aZ^2)-2)/3, 0, (2*(aZ^2)-1)/3];
dode18=[((aZ^2)-2)/3, 0, -(2*(aZ^2)-1)/3];
dode19=[-((aZ^2)-2)/3, 0, (2*(aZ^2)-1)/3];
dode20=[-((aZ^2)-2)/3, 0, -(2*(aZ^2)-1)/3];

dode_1=dode1/sqrt(dode1*dode1');
dode_2=dode2/sqrt(dode2*dode2');
dode_3=dode3/sqrt(dode3*dode3');
dode_4=dode4/sqrt(dode4*dode4');
dode_5=dode5/sqrt(dode5*dode5');
dode_6=dode6/sqrt(dode6*dode6');
dode_7=dode7/sqrt(dode7*dode7');
dode_8=dode8/sqrt(dode8*dode8');
dode_9=dode9/sqrt(dode9*dode9');
dode_10=dode10/sqrt(dode10*dode10');
dode_11=dode11/sqrt(dode11*dode11');
dode_12=dode12/sqrt(dode12*dode12');
dode_13=dode13/sqrt(dode13*dode13');
dode_14=dode14/sqrt(dode14*dode14');
dode_15=dode15/sqrt(dode15*dode15');
dode_16=dode16/sqrt(dode16*dode16');
dode_17=dode17/sqrt(dode17*dode17');
dode_18=dode18/sqrt(dode18*dode18');
dode_19=dode19/sqrt(dode19*dode19');
dode_20=dode20/sqrt(dode20*dode20');

dode=[dode_1' dode_2' dode_3' dode_4' dode_5' dode_6' dode_7' dode_8' dode_9' dode_10' dode_11' dode_12' dode_13' dode_14' dode_15' dode_16' dode_17' dode_18' dode_19' dode_20'];


dode_tet=acos(dode(3, :));
dode_fi=atan(dode(2, :)./dode(1, :));

cos_fi=dode(1, :)./sin(dode_tet);
% sin_fi=dode(2, :)./sin(dode_tet);
for iii=1:vert
    if (cos_fi(iii)<=0) 
        dode_fi(iii)=dode_fi(iii)+pi;
    end
end
      

dode_cos_tet_frak_2=cos(dode_tet/2);
dode_sin_tet_frak_2=sin(dode_tet/2);
dode_cos_fi_frak_2=cos(dode_fi/2);
dode_sin_fi_frak_2=sin(dode_fi/2);
psi_A_1=dode_cos_tet_frak_2.*(dode_cos_fi_frak_2-i*dode_sin_fi_frak_2);
psi_A_2=dode_sin_tet_frak_2.*(dode_cos_fi_frak_2+i*dode_sin_fi_frak_2);
psi_A=[psi_A_1;psi_A_2];

XXX=psi_A';
matrix=XXX;
matrix = matrix*rotation_operator(pi/3, pi/3, pi/3);
%matrix = [1 0; 0 1; 1/sqrt(2)  1/sqrt(2); 1/sqrt(2)  (-1)/sqrt(2); 1/sqrt(2)  (-sqrt(-1))/sqrt(2); 1/sqrt(2)  (sqrt(-1))/sqrt(2)];

[X,Y,Z] = sphere(80);
LL = zeros(81, 81);
n1 = sin(pi/2.5) * cos(pi/2.5);
n2 = sin(pi/2.5) * sin(pi/2.5);
n3 = cos(pi/2.5);
sig_1 = [ 0 1;
          1 0 ];
sig_2 = [ 0 -1i;
          1i 0 ];
sig_3 = [ 1 0;
          0 -1 ];
matrix1 = matrix* expm(-1i * pi/2.5 / 2 * (n1 * sig_1 + n2 * sig_2 + n3 * sig_3));
C1test=zeros(81,1);
C2test=zeros(81,1);
CCN1=zeros(2,81);
CCN=zeros(4,81);
cmap = colormap('jet');
for q0 = 1:81
    for q00 = 1:81
        [phi1,tet1] = cart2sph(X(q0,q00),Y(q0,q00),Z(q0,q00));
        tet1 = pi/2 - tet1;
        C1test(q00,1)=cos(tet1/2);
        C2test(q00,1)=sin(tet1/2)*exp(1j*phi1);
        CCN1(:,q00)=[C1test(q00,1);C2test(q00,1)];
        CCN(:,q00)=[real(C1test(q00,1));real(C2test(q00,1)); imag(C1test(q00,1)); imag(C2test(q00,1))];
        phi = zeros(20,1);
        theta = zeros(20,1);
        s = zeros(20,1);
        for k=1:20
            phi(k) = angle(matrix1(k,2))-angle(matrix1(k,1));
            theta(k) = 2*acos(abs(matrix1(k,1)));
            s(k) = (angle(matrix1(k,1))+angle(matrix(k,2)))/2;
        end
        H = zeros(4,4);
        for k=1:20
            AAA = [1+cos(theta(k)),sin(theta(k))*exp(1i*phi(k));
                     sin(theta(k))*exp(-1i*phi(k)),1-cos(theta(k))];
            AAA = AAA*1/4;
            AAAA = [real(AAA), -imag(AAA); imag(AAA), real(AAA)];
            if( abs(CCN(:,q00)'*AAAA*CCN(:,q00)) > 10^(-30)) 
                 H = H + (AAAA*CCN(:,q00)*CCN(:,q00)'*AAAA )/(CCN(:,q00)'*AAAA*CCN(:,q00)); 
            end
        end
        H = 4000*2* H;
        S = eig(H);
            S = sort(S);
            d = zeros(length(S)-2, 1);
            for k=1:length(d)
               d(k) = 1 / (2 * S(k + 1)); 
            end
        dF = sum(d);
        LL(q0,q00) =4000* dF;
    end
end
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
max1 = max(max(LL));
min1 = min(min(LL));
LL = round((LL - min1) / (max1-min1) * size(cmap,1));
LL(LL==0) = 1;
C= zeros([size(LL) 3]);
for i = 1:size(LL,1)
    for j = 1:size(LL,2)
        color = cmap(LL(i,j),:);
        C(i,j,1) = color(1);
        C(i,j,2) = color(2);
        C(i,j,3) = color(3);
    end
end
h = surf(X,Y,Z,C, 'EdgeAlpha', 0.3);
colorbar;
caxis([min1 max1]);
title(sprintf('Ikosaedr, Lmin = %.3f, Lmax = %.3f', min1, max1));
