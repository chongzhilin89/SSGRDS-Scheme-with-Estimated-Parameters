//------------------------------------------------------------------------------
function [EANOS1]=rlSSGRDSe(m,n,n1,n2,L1,L2,L3,L,dmin,dmax)
//------------------------------------------------------------------------------
[argout,argin]=argn()
if (argin~=10)
  error("incorrect number of arguments")
end
if (m<=0)|(m~=floor(m))
  error("argument ''m'' must be an integer >= 1")
end
if (n<=0)|(n~=floor(n))
  error("argument ''n'' must be an integer >= 1")
end
if (n1<=0)|(n1~=floor(n1))
  error("argument ''n1'' must be an integer >= 1")
end
if (n2<=0)|(n2~=floor(n2))
  error("argument ''n2'' must be an integer >= 1")
end
if L1<=0
  error("argument ''L1'' must be > 0")
end
if L2<L1
  error("argument ''L2'' must be >= L1")
end
if L2<=0
  error("argument ''L2'' must be > 0")
end
if dmin<0
  error("argument ''deltamin'' must be >= 0")
end
if dmax<0
  error("argument ''deltamax'' must be >= 0")
end
h=15
[u,wu]=quadhermite(h)
a=m*(n-1)/2
b=2/(m*(n-1))
sv=sqrt((a-exp(2*(gammaln(a+0.5)-gammaln(a))))*b)
vmin=max(0,1-5*sv)
vmax=1+5*sv
[v,wv]=quadlegendre(h,vmin,vmax)
U=ones(h,1)*u'
V=v*ones(1,h)
W=(ones(h,1)*wu').*(wv*ones(1,h))
FU=pdfnormal(U)
FV=2*V.*pdfgamma(V.^2,m*(n-1)/2,2/(m*(n-1)))
sn1=sqrt(n1)
sn2=sqrt(n2)
sn1mn=sqrt(n1/(m*n))
sn2mn=sqrt(n2/(m*n))
sn12=sqrt(n1+n2)
[xi,we]=lgwt_table(9)
EANOS1sum=0
for id=1:9
    delta=((dmax-dmin)/2)*xi(id)+((dmax+dmin)/2)
P1=cdfnormal(U*sn1mn+V*L1-delta*sn1)-cdfnormal(U*sn1mn-V*L1-delta*sn1)
P2=zeros(P1)
Pal=zeros(P1)
[za,wza]=quadlegendre(h,-L2,-L1)
[zb,wzb]=quadlegendre(h,L1,L2)
for iu=1:h
  uiu=u(iu)
  for iv=1:h
    viv=v(iv)
    fza=viv*pdfnormal(uiu*sn1mn+viv*za-delta*sn1)
    fzb=viv*pdfnormal(uiu*sn1mn+viv*zb-delta*sn1)
    p21=sum(wza.*(cdfnormal(uiu*sn2mn+viv*(L3*sn12-za*sn1)/sn2-delta*sn2)-..
                  cdfnormal(uiu*sn2mn-viv*(L3*sn12+za*sn1)/sn2-delta*sn2)).*fza)
    p22=sum(wzb.*(cdfnormal(uiu*sn2mn+viv*(L3*sn12-zb*sn1)/sn2-delta*sn2)-..
                  cdfnormal(uiu*sn2mn-viv*(L3*sn12+zb*sn1)/sn2-delta*sn2)).*fzb)
    Pal(iv,iu)=sum(wzb.*(cdfnormal(uiu*sn2mn+viv*(L3*sn12-zb*sn1)/sn2-delta*sn2)).*fzb)
    P2(iv,iu)=p21+p22
  end
end
P=1-P1-P2
A=1-(1-P).^L
A2=A.^2
al=(1-cdfnormal(U*sn1mn+V*L1-delta*sn1)-Pal)./P
ARLP=((1-al.*(1-al).*A2)./((1+al.*(1-al).*(A-2)).*A2))./P
ARL=sum(W.*FU.*FV.*ARLP)
P3=cdfnormal(U*sn1mn-V*L1-delta*sn1)-cdfnormal(U*sn1mn-V*L2-delta*sn1)+..
   cdfnormal(U*sn1mn+V*L2-delta*sn1)-cdfnormal(U*sn1mn+V*L1-delta*sn1)
n=n1+n2*P3
A=n.*ARLP
ANOS1=sum(W.*FU.*FV.*A)
EANOS1=we(id)*(1/(dmax-dmin))*ANOS1
EANOS1sum=EANOS1sum+EANOS1
end
EANOS1=((dmax-dmin)/2)*EANOS1sum

//------------------------------------------------------------------------------
function [ANOS]=rlSSGRDS0e(m,n,n1,n2,L1,L2,L3,L,delta)
//------------------------------------------------------------------------------
[argout,argin]=argn()
if (argin<8)|(argin>9)
  error("incorrect number of arguments")
end
if (m<=0)|(m~=floor(m))
  error("argument ''m'' must be an integer >= 1")
end
if (n<=0)|(n~=floor(n))
  error("argument ''n'' must be an integer >= 1")
end
if (n1<=0)|(n1~=floor(n1))
  error("argument ''n1'' must be an integer >= 1")
end
if (n2<=0)|(n2~=floor(n2))
  error("argument ''n2'' must be an integer >= 1")
end
if L1<=0
  error("argument ''L1'' must be > 0")
end
if L2<L1
  error("argument ''L2'' must be >= L1")
end
if L2<=0
  error("argument ''L2'' must be > 0")
end
if argin==8
  delta=0
end
if delta<0
  error("argument ''delta'' must be >= 0")
end
h=15
[u,wu]=quadhermite(h)
a=m*(n-1)/2
b=2/(m*(n-1))
sv=sqrt((a-exp(2*(gammaln(a+0.5)-gammaln(a))))*b)
vmin=max(0,1-5*sv)
vmax=1+5*sv
[v,wv]=quadlegendre(h,vmin,vmax)
U=ones(h,1)*u'
V=v*ones(1,h)
W=(ones(h,1)*wu').*(wv*ones(1,h))
FU=pdfnormal(U)
FV=2*V.*pdfgamma(V.^2,m*(n-1)/2,2/(m*(n-1)))
sn1=sqrt(n1)
sn2=sqrt(n2)
sn1mn=sqrt(n1/(m*n))
sn2mn=sqrt(n2/(m*n))
sn12=sqrt(n1+n2)
P1=cdfnormal(U*sn1mn+V*L1-delta*sn1)-cdfnormal(U*sn1mn-V*L1-delta*sn1)
P2=zeros(P1)
Pal=zeros(P1)
[za,wza]=quadlegendre(h,-L2,-L1)
[zb,wzb]=quadlegendre(h,L1,L2)
for iu=1:h
  uiu=u(iu)
  for iv=1:h
    viv=v(iv)
    fza=viv*pdfnormal(uiu*sn1mn+viv*za-delta*sn1)
    fzb=viv*pdfnormal(uiu*sn1mn+viv*zb-delta*sn1)
    p21=sum(wza.*(cdfnormal(uiu*sn2mn+viv*(L3*sn12-za*sn1)/sn2-delta*sn2)-..
                  cdfnormal(uiu*sn2mn-viv*(L3*sn12+za*sn1)/sn2-delta*sn2)).*fza)
    p22=sum(wzb.*(cdfnormal(uiu*sn2mn+viv*(L3*sn12-zb*sn1)/sn2-delta*sn2)-..
                  cdfnormal(uiu*sn2mn-viv*(L3*sn12+zb*sn1)/sn2-delta*sn2)).*fzb)
    Pal(iv,iu)=sum(wzb.*(cdfnormal(uiu*sn2mn+viv*(L3*sn12-zb*sn1)/sn2-delta*sn2)).*fzb)
    P2(iv,iu)=p21+p22
  end
end
P=1-P1-P2
A=1-(1-P).^L
A2=A.^2
al=(1-cdfnormal(U*sn1mn+V*L1-delta*sn1)-Pal)./P
ARLP=((1-al.*(1-al).*A2)./((1+al.*(1-al).*(A-2)).*A2))./P
ARL=sum(W.*FU.*FV.*ARLP)
P3=cdfnormal(U*sn1mn-V*L1-delta*sn1)-cdfnormal(U*sn1mn-V*L2-delta*sn1)+..
   cdfnormal(U*sn1mn+V*L2-delta*sn1)-cdfnormal(U*sn1mn+V*L1-delta*sn1)
n=n1+n2*P3
A=n.*ARLP
ANOS=sum(W.*FU.*FV.*A)
//------------------------------------------------------------------------------
function ASS=asSSGRDSe(m,n,n1,n2,L1,L2,delta)
//------------------------------------------------------------------------------
[argout,argin]=argn()
if (argin<6)|(argin>7)
  error("incorrect number of arguments")
end
if (m<=0)|(m~=floor(m))
  error("argument ''m'' must be an integer >= 1")
end
if (n<=0)|(n~=floor(n))
  error("argument ''n'' must be an integer >= 1")
end
if (n1<=0)|(n1~=floor(n1))
  error("argument ''n1'' must be an integer >= 1")
end
if (n2<=0)|(n2~=floor(n2))
  error("argument ''n2'' must be an integer >= 1")
end
if L1<=0
  error("argument ''L1'' must be > 0")
end
if L2<L1
  error("argument ''L2'' must be >= L1")
end
if argin==6
  delta=0
end
if delta<0
  error("argument ''delta'' must be >= 0")
end
h=15
[u,wu]=quadhermite(h)
a=m*(n-1)/2
b=2/(m*(n-1))
sv=sqrt((a-exp(2*(gammaln(a+0.5)-gammaln(a))))*b)
vmin=max(0,1-5*sv);
vmax=1+5*sv
[v,wv]=quadlegendre(h,vmin,vmax)
U=ones(h,1)*u'
V=v*ones(1,h)
W=(ones(h,1)*wu').*(wv*ones(1,h))
FU=pdfnormal(U)
FV=2*V.*pdfgamma(V.^2,m*(n-1)/2,2/(m*(n-1)))
sn1=sqrt(n1)
sn1mn=sqrt(n1/(m*n))
sn2mn=sqrt(n2/(m*n))
P3=cdfnormal(U*sn1mn-V*L1-delta*sn1)-cdfnormal(U*sn1mn-V*L2-delta*sn1)+..
   cdfnormal(U*sn1mn+V*L2-delta*sn1)-cdfnormal(U*sn1mn+V*L1-delta*sn1)
n=n1+n2*P3
ASS=sum(W.*FU.*FV.*n)
//------------------------------------------------------------------------------
function dif=optdsxbareL2(L2,m,n,n1,n2,L1)
//------------------------------------------------------------------------------
if L2<L1
  dif=%inf
else
  ASS=asSSGRDSe(m,n,n1,n2,L1,L2)
  dif=n-ASS
end

//------------------------------------------------------------------------------
function dif=optdsxbareL3(L3,m,n,n1,n2,L1,L2,L,anos0)
//------------------------------------------------------------------------------
if L3<=0
  dif=%inf
else
  ANOS=rlSSGRDS0e(m,n,n1,n2,L1,L2,L3,L)
  dif=anos0-ANOS
end

//------------------------------------------------------------------------------
function dif=optdsxbareL1min(L1,m,n,n1,n2,L1max,anos0)
//------------------------------------------------------------------------------
if (L1<=0)|(L1>=L1max)
  dif=%inf
else
  L2=simplexolve(2*L1,optdsxbareL2,list(m,n,n1,n2,L1),tol=1e-6)
  ANOS=rlSSGRDS0e(m,n,n1,n2,L1,L2,1000,1000)
  dif=anos0-ANOS
end

//------------------------------------------------------------------------------
function dif=optdsxbareL1max(L1,m,n,n1,n2)
//------------------------------------------------------------------------------
if L1<=0
  dif=%inf
else
  ASS=asSSGRDSe(m,n,n1,n2,L1,1000)
  dif=n-ASS
end

//------------------------------------------------------------------------------
function EANOS1=optdsxbareL1(L1,m,n,n1,n2,dmin,dmax,L1min,L1max,anos0)
//------------------------------------------------------------------------------
if (L1<=L1min)|(L1>=L1max)
  EANOS1=%inf
else
  L2=simplexolve(2*L1,optdsxbareL2,list(m,n,n1,n2,L1),tol=1e-6)
  L3=simplexolve(1,optdsxbareL3,list(m,n,n1,n2,L1,L2,anos0),tol=1e-6)
  EANOS1=rlSSGRDSe(m,n,n1,n2,L1,L2,L3,1000,dmin,dmax)
end
//------------------------------------------------------------------------------
function [L1,L2,L3,L]=optdsxbare(m,n,n1,n2,anos0,dmin,dmax)
//------------------------------------------------------------------------------
[argout,argin]=argn()
if (argin~=7)
  error("incorrect number of arguments")
end
if (m<=0)|(m~=floor(m))
  error("argument ''m'' must be an integer >= 1")
end
if (n<=0)|(n~=floor(n))
  error("argument ''n'' must be an integer >= 1")
end
if (n1<=0)|(n1~=floor(n1))
  error("argument ''n1'' must be an integer >= 1")
end
if (n2<=0)|(n2~=floor(n2))
  error("argument ''n2'' must be an integer >= 1")
end
if (n<=n1)|(n>=n1+n2)
  error("arguments ''n,n1,n2'' must verify n1<n<n1+n2")
end
if anos0<=0
  error("argument ''anos0'' must be >= 1")
end
if dmin<0
  error("argument ''deltamin'' must be >= 0")
end
if dmax<0
  error("argument ''deltamax'' must be >= 0")
end
L1max=simplexolve(1,optdsxbareL1max,list(m,n,n1,n2),tol=1e-6)
L1min=simplexolve(L1max/2,optdsxbareL1min,list(m,n,n1,n2,L1max,anos0),tol=1e-6)
L1=neldermead((L1min+L1max)/2,optdsxbareL1,list(m,n,n1,n2,dmin,dmax,L1min,L1max,anos0),tol=1e-6,opt="min")
L2=simplexolve(2*L1,optdsxbareL2,list(m,n,n1,n2,L1),tol=1e-6)
L=1
EANOSo=%inf
while %t
L3=simplexolve(1,optdsxbareL3,list(m,n,n1,n2,L1,L2,L,anos0),tol=1e-6)
EANOS1=rlSSGRDSe(m,n,n1,n2,L1,L2,L3,L,dmin,dmax)
if EANOS1>=EANOSo
    L=L-1
    break
  else
    L=L
    L3o=L3
    L30=L3
    L=L+1
    EANOSo=EANOS1
  end
end
L=L
L3=L3o
//------------------------------------------------------------------------------
function [n1,n2,L1,L2,L3,L,EANOSo]=opt_est(m,n,anos0,dmin,dmax)
//------------------------------------------------------------------------------
fil=sprintf("n%d_dmin%02d_dmax%02d_m%2d.txt",n,dmin*10,dmax*10,m)
f=mopen(fil,"w")
EANOSo=%inf
  for n1=1:n-1
    for n2=n-n1+1:2*n
      [L1,L2,L3,L]=optdsxbare(m,n,n1,n2,anos0,dmin,dmax)
      EANOS1=rlSSGRDSe(m,n,n1,n2,L1,L2,L3,L,dmin,dmax)
      mprintf("%2d %2d %9.6f %9.6f %9.6f %11.6f %9.6f \n",[n1,n2,L1,L2,L3,L,EANOS1])
      mfprintf(f,"%2d %2d %9.6f %9.6f %9.6f %11.6f %9.6f\n",[n1,n2,L1,L2,L3,L,EANOS1])
      if EANOS1<EANOSo
        EANOSo=EANOS1
        sol=[n1,n2,L1,L2,L3,L,EANOSo]
      end
    end
  end

  n1=sol(1)
  n2=sol(2)
  L1=sol(3)
  L2=sol(4)
  L3=sol(5)
  L=sol(6)
  mprintf("%2d %2d %6.5f %6.5f %6.5f %6.1f %6.5f \n",[n1,n2,L1,L2,L3,L,EANOSo])
  mfprintf(f,"%2d %2d %6.5f %6.5f %6.5f %6.1f %6.5f \n",[n1,n2,L1,L2,L3,L,EANOSo])

mclose(f)


//------------------------------------------------------------------------------
function [x,w]=lgwt_table(n)
//------------------------------------------------------------------------------

if n==9      
 x=[-0.9681602395;
    -0.8360311073;
    -0.6133714327;
    -0.3242534234;
     0.0;
     0.3242534234;
     0.6133714327;
     0.8360311073;
     0.9681602395]

 w=[0.0812743884;
    0.1806481607;
    0.2606106964;
    0.3123470770
    0.3302393550;
    0.3123470770;
    0.2606106964;
    0.1806481607;
    0.0812743884]
end
