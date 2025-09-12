#n = 1,2...
#l = 0,1...n-1
#m = -l...l

n=4;
l=3;
m=0;

function tpart = THETA(theta,m)
  if (m==0)
    tpart = 1/sqrt(2*pi);
  elseif (m>0)
    tpart = cos(m*theta)./sqrt(2*pi);
  else
    tpart = sin(m*theta)./sqrt(2*pi);
  endif;
 endfunction;

lmod = mod(l, 2);

al=1;
if(l>1)
    if (lmod==0)
        for k = 1:l/2
            j=2*k;
            al=[al,al(k)*((j-2)*(j-1)-l*(l+1))/(j*(j-1))];
        endfor;
    else
        for k = 1:(l-1)/2
        j=2*k+1;
        al=[al,al(k)*((j-2)*(j-1)-l*(l+1))/(j*(j-1))];
        endfor;
    endif;
endif;

sigma=@(x,f,a,b)sum(arrayfun(f,x*ones(1,b-a+1),a:b));

derivfacms=[];
for q=0:floor(l/2)
  derivfacms = [derivfacms,prod((2*q + lmod - abs(m)+1):(2*q + lmod))];
endfor

AL=@(x,k)al(k+1)*derivfacms(k+1)*x^(2*k + lmod - abs(m));
slegendre=@(w)sigma(w,AL,0,floor(l/2));
Legendre=@(w)arrayfun(slegendre,w);

Plm=@(w)(1-w.^2).^(abs(m)/2);

PHI=@(phi)Plm(cos(phi)).*Legendre(cos(phi));

u=9.1*10^-31;
Z=1;
e0=1.6*10^-19;
h=1.1*10^-34;
eps0=8.9*10^-12;

E=-u*Z^2*e0^4/(2*(4*pi*eps0*h*n)^2);

af = 1;
if (l<n-1)
    for k=1:n-1-l
        af=[af,af(k)*(l+k-n)/(k*(k+2*l+1))];
    endfor;
endif;

prob=@(x)abs(x).^2;

AF=@(x,k)af(k+1)*x^k;
sf=@(x)sigma(x,AF,0,n-1-l);
F=@(rho)arrayfun(sf,rho);

rho=@(x)x*sqrt(8*u*abs(E))/h;
U=@(x)x.^(l+1).*e.^(-x/2).*F(x);

R=@(r)U(rho(r))./r;

W=@(r,t,p)R(r).*THETA(t,m).*PHI(p);
W1=@(r,t,p)prob(W(r,t,p));

cart2rho=@(x,y,z)sqrt(x.^2+y.^2+z.^2);
cart2thet=@(x,y,z)atan2(y,x);
cart2phi=@(x,y,z)atan2(sqrt(x.^2+y.^2),z);

W1cart=@(x,y,z)W1(cart2rho(x,y,z),cart2thet(x,y,z),cart2phi(x,y,z));

width =  1.3*10^-9;
distr = 25;

tx=ty=tz=linspace(-width,width,distr);
[xx,yy,zz]=meshgrid(tx,ty,tz);

Ws1 = W1cart(xx,yy,zz);
approxscale = max(max(max(Ws1)));
Ws2 = max(Ws1./approxscale,0.001);
Ws2s = Ws2.*20;

scatter3(xx(:),yy(:),zz(:), Ws2s(:), Ws2(:), 'filled');
colormap hot
colorbar
axis equal tight

#n=1, width = 1.5*10^-10;
#n=2, width = 8*10^-10;
#n=3, width = 1.3*10^-9;

#{
integralarea=triplequad(@(r,t,p)W1(r,t,p).*r.^2*sin(p),0,10^-9,0,2*pi,0,pi,10^-11)
scalefactor = sqrt(integralarea)

Wscaled=@(r,t,p)W(r,t,p)./scalefactor;
#}

#{
x=0:10^-10:10^-8;
plot(x,R(x))
axis("equal","auto");
#}


