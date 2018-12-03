clear
clc
XXX=[0.908136693119174+0.0i 	0.1510131300963873+0.39049043669347744i; 0.4549936260526747+0.0i 	-0.8603142221076049-0.22987004913825868i; 0.7330605004147114+0.0i 	-0.6654905652320535+0.1405155162921693i ;0.8183040579168933+0.0i 	0.40985938451857046-0.4029810835744083i ;0.9193004590602009+0.0i 	-0.016799455916566616-0.39319771648956886i; 
0.499999999280608+0.0i 	0.617532256756284-0.6071687677943283i ;
0.5000000004673837+0.0i 	-0.6175322560473718+0.6071687675380375i ;
0.7330605005724823+0.0i 	-0.1517574249019325+0.6630173349826955i ;
0.7071067813678041+0.0i 	0.49575122277374883+0.5042129756977674i ;
0.9081366930112436+0.0i 	-0.3878789375650467-0.1575997354059104i ;
0.906930947136282+0.0i 	-0.3003994592373508+0.2953581250218008i ;
0.3427958123053995+0.0i 	0.9188847901015343+0.19529918992648299i; 
0.3427958118243274+0.0i 	-0.17972102583061966-0.9220582325808333i ;
0.45499362662218623+0.0i 	0.2152781197820067+0.864081090451767i ;
0.08862688902268234+0.0i 	-0.7102588413229527+0.6983392104591706i; 
0.707106780984806+0.0i 	-0.49575122300714847-0.5042129760054002i ;
0.6843015026807708+0.0i 	0.014033423985478315-0.7290641373981324i ;
0.9193004591742169+0.0i 	0.3928571138301035+0.023451095393501493i ;
0.6843015029813213+0.0i 	0.7291972195207386-0.0016935349799648031i;
1.0+0.0i 	0.0+0.0i ];

matrix=XXX;
matrix = matrix*rotation_operator(pi/4, pi/4, pi/4);
[X,Y,Z] = sphere(80);
LL = zeros(81, 81);
n1 = sin(pi/4) * cos(pi/4);
n2 = sin(pi/4) * sin(pi/4);
n3 = cos(pi/4);
sig_1 = [ 0 1;
          1 0 ];
sig_2 = [ 0 -1i;
          1i 0 ];
sig_3 = [ 1 0;
          0 -1 ];
matrix1 = matrix* expm(-1i * pi/4 / 2 * (n1 * sig_1 + n2 * sig_2 + n3 * sig_3));
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
        H =  2*4000*H;
        S = eig(H);
            S = sort(S);
            d = zeros(length(S)-2, 1);
            for k=1:length(d)
               d(k) = 1 / (2 * S(k + 1)); 
            end
        dF = sum(d);
        LL(q0,q00) =4000*dF;
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
title(sprintf('Prot_Bel (20 st), Lmin = %.3f, Lmax = %.3f', min1, max1));
