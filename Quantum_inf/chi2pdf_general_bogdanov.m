function p = chi2pdf_general_bogdanov(x, d)
err = 1;
N0=2048*2;
z0=8;
while err>3e-3
    [x_F,p,err]=chi2_like_pdf_new(d,N0,z0);
    N0=2*N0; z0=z0+2;
    if N0>1e7
        break;
    end
end
p = interp1(x_F, p, x);
p(isnan(p)) = 0;
    
end

function [x_new,p_new,err]=chi2_like_pdf_new(d,N,z)

d=1./sort(1./d);
s=prod(size(d));
mju=sum(d);
sigma=sqrt(2*sum(d.^2));
x_max=mju+z*sigma;
dx=x_max/N;
du=2*pi/x_max;
j=(1:N)';
k=j;
u=(j-1)*du;
x=(k-1)*dx;
fi=1./sqrt(1-2*i*u*d(1));
for j1=2:s
fi=fi.*(1./sqrt(1-2*i*u*d(j1)));   
end
p=(du/(1*pi))*real(fft(fi));
x=x(1:N-20);
p=p(1:N-20)-min(p);

Norm=sum(p.*dx);
x_mean=sum(x.*p.*dx);
x_2_mean=sum(x.^2.*p.*dx);
sig=sqrt(x_2_mean-x_mean^2);
x_mean_teor=sum(d);
sig_teor=sqrt(2*sum(d.^2));
err1=abs(1-Norm);
err2=abs(x_mean_teor-x_mean)/x_mean_teor;
err3=abs(sig_teor-sig)/sig_teor;

err=max([err1 err2 err3]);
skew=sum((x-x_mean).^3.*p.*dx)/sig^3;
skew_teor=8*sum(d.^3)/sig_teor^3;
exc=sum((x-x_mean).^4.*p.*dx)/sig^4-3;
exc_teor=48*sum(d.^4)/sig_teor^4;

x_new=x;
p_new=p;

end
