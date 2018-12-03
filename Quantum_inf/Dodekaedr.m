clear
clc
vert=12;
tau=(1+sqrt(5))/2;
ikos1=[1, tau, 0];
ikos2=[-1, tau, 0];
ikos3=[1, -tau, 0];
ikos4=[-1, -tau, 0];
ikos5=[0, 1, tau];
ikos6=[0, -1, tau];
ikos7=[0, 1, -tau];
ikos8=[0, -1, -tau];
ikos9=[tau, 0, 1];
ikos10=[tau, 0, -1];
ikos11=[-tau, 0, 1];
ikos12=[-tau, 0, -1];

ikos_1=ikos1/sqrt(ikos1*ikos1');
ikos_2=ikos2/sqrt(ikos2*ikos2');
ikos_3=ikos3/sqrt(ikos3*ikos3');
ikos_4=ikos4/sqrt(ikos4*ikos4');
ikos_5=ikos5/sqrt(ikos5*ikos5');
ikos_6=ikos6/sqrt(ikos6*ikos6');
ikos_7=ikos7/sqrt(ikos7*ikos7');
ikos_8=ikos8/sqrt(ikos8*ikos8');
ikos_9=ikos9/sqrt(ikos9*ikos9');
ikos_10=ikos10/sqrt(ikos10*ikos10');
ikos_11=ikos11/sqrt(ikos11*ikos11');
ikos_12=ikos12/sqrt(ikos12*ikos12');


ikos=[ikos_1' ikos_2' ikos_3' ikos_4' ikos_5' ikos_6' ikos_7' ikos_8' ikos_9' ikos_10' ikos_11' ikos_12'];
ikos_tet=acos(ikos(3, :));
ikos_fi=atan(ikos(2, :)./ikos(1, :));

cos_fi=ikos(1, :)./sin(ikos_tet);
for iii=1:vert
    if (cos_fi(iii)<=0) 
        ikos_fi(iii)=ikos_fi(iii)+pi;
    end
end

ikos_cos_tet_frak_2=cos(ikos_tet/2);
ikos_sin_tet_frak_2=sin(ikos_tet/2);
ikos_cos_fi_frak_2=cos(ikos_fi/2);
ikos_sin_fi_frak_2=sin(ikos_fi/2);
psi_A_1=ikos_cos_tet_frak_2.*(ikos_cos_fi_frak_2-i*ikos_sin_fi_frak_2);
psi_A_2=ikos_sin_tet_frak_2.*(ikos_cos_fi_frak_2+i*ikos_sin_fi_frak_2);
psi_A=[psi_A_1;psi_A_2];

XXX=psi_A';


%matrix = [1 0; 0 1; 1/sqrt(2)  1/sqrt(2); 1/sqrt(2)  (-1)/sqrt(2); 1/sqrt(2)  (-sqrt(-1))/sqrt(2); 1/sqrt(2)  (sqrt(-1))/sqrt(2)];
matrix=XXX;
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
        phi = zeros(12,1);
        theta = zeros(12,1);
        s = zeros(12,1);
        for k=1:12
            phi(k) = angle(matrix1(k,2))-angle(matrix1(k,1));
            theta(k) = 2*acos(abs(matrix1(k,1)));
            s(k) = (angle(matrix1(k,1))+angle(matrix(k,2)))/2;
        end
        H = zeros(4,4);
        for k=1:12
            AAA = [1+cos(theta(k)),sin(theta(k))*exp(1i*phi(k));
                     sin(theta(k))*exp(-1i*phi(k)),1-cos(theta(k))];
            AAA = AAA*1/4;
            AAAA = [real(AAA), -imag(AAA); imag(AAA), real(AAA)];
            if( abs(CCN(:,q00)'*AAAA*CCN(:,q00)) > 10^(-30)) 
                 H = H + (AAAA*CCN(:,q00)*CCN(:,q00)'*AAAA )/(CCN(:,q00)'*AAAA*CCN(:,q00)); 
            end
        end
        H =  2*4000* H;
        S = eig(H);
            S = sort(S);
            d = zeros(length(S)-2, 1);
            for k=1:length(d)
               d(k) = 1 / (2 * S(k + 1)); 
            end
        dF = sum(d);
        LL(q0,q00)= 4000*dF;
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
h = surf(X,Y,Z,C,'EdgeAlpha', 0.3);
colorbar;
caxis([min1 max1]);
title(sprintf('Dodekaedr, Lmin = %.3f, Lmax = %.3f', min1, max1));
