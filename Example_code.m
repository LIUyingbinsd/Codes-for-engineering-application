%% Example code used for heave 1 induced by grouting and torsion induced by uneven heave 1
clear;clc;close all;
%% parameter
E=30e9;vs=0.2;
E1=E/(1-vs^2);G=E/2/(1+vs);
Es = 7e6;v=0.3;% soil parameter
H = 6.7;% Tunnel depth
R = 3.5/2;
t = 0.25;
Ac = pi*(R^2-(R-t)^2);
I=pi/64*(R^4-(R-t)^4)*2^4;
L=200; % calculation length
l=0.5; % element length
m=L/l; 
C=cell(m);
xint=zeros(m+1,1);
j=0;
fb = 1/7;fs = 1/7;ft = 0.04;% effective factor for tunnel stiffness
kn = 1.3*Es/(1+1/(1.7*(H)/2/R))/2/R/(1-v^2)*(Es*(2*R)^4/(fb*E*I))^(1/12);% soil coefficient
kt = kn/3;
q = 1.2e6;% air pressure
d = 2.45 + R/2;% tunnel center-grouting center distance
Gs = 2.692e6;% fill
G2 = 22.69e6;% CDG
y1 = -1.5/2;y2 = -y1;Rg = y2;% column diameter
%% Define soil reaction matrix and tunnel stiffness matrix
for i=1:m   
    C{1,i}=Node_Stiffness(E1,G,R,t,l,kn,kt,fb,fs,ft);%soil reaction matrix
    xint(i)=j*l;
    j=j+1;
end
xint(m+1)=L;
KK=zeros(3*(m+1),3*(m+1));
for i=1:m
    KK=Node_Assembly(KK,C{1,i},i,i+1);
end
%% Define force matrix
D=cell(m);
h1 = 10.4;% fill depth
for i=1:m
    y0 = i*l-L/2;
    if -21/2 <= y0 && y0<=21/2  % limited force range imposed by permanent wall
    dh1 = -7.2;% distance between tunnel and grouting center
    k1 = 0.6;% reduction factor
    z2 = 22.7;z1= 1.9; % Grouting range
    z2t = (z2 - h1)*(G2/Gs)^1/3 + h1;% effective depth
    D{1,i} = Node_eqforce(y0,z1,h1,v,Gs,k1*q,d,dh1,H,R,l,kn,kt,Rg)...
         + Node_eqforce(y0,h1,z2t,v,Gs,k1*q*(z2-h1)/(z2t-h1),d,dh1,H,R,l,kn,kt,Rg);%
    else
        D{1,i} = [0;0;0;0;0;0];
    end
end
fq=zeros(3*(m+1),1);
for i=1:m
    fq=Node_FAssembly(fq,D{1,i},i,i+1);
end
%% Calculation
n=3*(m+1);
k=KK(4:n-3,4:n-3);
p=fq(4:n-3);
u=k\p;
u=[0;0;0;u;0;0;0;];% Fixed ends
F=KK*u;
b=u(1:3:end);v=u(2:3:end);o=u(3:3:end);% calculated deformation,bending,shearing,torsion
Q=[];M=[];T=[];
for i=1% Force computation
    ki=C{1,i};ui=u(3*i-2:3*i+3);P=ki*ui;
    M(i)=P(1);Q(i)=P(2);T(i)=P(3);
end
for i=2:m
    ki=C{1,i};ui=u(3*i-2:3*i+3);Fm=D{1,i};P=ki*ui-Fm;
    M(i)=P(1);Q(i)=P(2);T(i)=P(3);
end
for i=m%
    ki=C{1,i};ui=u(3*i-2:3*i+3);P=-ki*ui;
    M(i+1)=P(4);Q(i+1)=P(5);T(i+1)=P(6);
end

%% matrix formation
function k=Node_Stiffness(E,G,R,t,l,kn,kt,fb,fs,ft)
koo1=fs*pi*G*R*t;koo2=pi*E*R^3*t*fb;kvv=koo1;kov=fs*pi*G*R*t;
kpp=2*G*pi*R^3*t*ft;
kvv1=2*kn*R*pi/4;kvv3=2*kt*R*pi/4;
koo=kt*pi*R^3;
kp=2*kt*pi*R^3;
k=[koo1*l/3+koo2/l  -1/2*kov   0   koo1*l/6-koo2/l   1/2*kov  0;
    -1/2*kov   kvv/l  0      -1/2*kov    -kvv/l  0;
    0        0       kpp/l     0        0       -kpp/l;
    koo1*l/6-koo2/l  -1/2*kov   0  koo1*l/3+koo2/l  1/2*kov 0 ;
    1/2*kov   -kvv/l  0    1/2*kov    kvv/l    0;
    0        0       -kpp/l     0        0       kpp/l;];
ks=[2*koo  0  0   koo   0    0;
    0    2*(kvv1+kvv3)  0   0  (kvv1+kvv3)   0;
    0    0     2*kp   0    0     kp;
    koo  0   0  2*koo  0 0;
    0   (kvv1+kvv3)   0       0    2*(kvv1+kvv3)  0;
    0    0     kp   0    0     2*kp;];
k=k+l/6*ks;
end

function Z=Node_Assembly(KK,k,i,j)
DOF(1)=3*i-2;DOF(2)=3*i-1;DOF(3)=3*i;DOF(4)=3*j-2;DOF(5)=3*j-1;DOF(6)=3*j;
for n1=1:6
    for n2=1:6
        KK(DOF(n1),DOF(n2))=KK(DOF(n1),DOF(n2))+k(n1,n2);
    end
end
Z=KK;
end
% Considering the complexity and difficulty of multiple integrations, to simplify the calculation, the load is regarded as the vertical load determined by three points at the bottom.
function fe=Node_eqforce(y0,z1,z2,v,G,q,D,dh,H,R,l,kn,kt,Rg)
W= (JGv(y0,z1,z2,v,G,q,D,dh,H,Rg) + (JGv(y0,z1,z2,v,G,q,D-R,dh,H,Rg) ...
    + JGv(y0,z1,z2,v,G,q,D+R,dh,H,Rg))/2)*kn*R;
Fv1=@(z)(1-z/l).*W;
Fv2=@(z)z/l.*W;
f1=integral(Fv1,0,l);f4=integral(Fv2,0,l);
q2 = JGv(y0,z1,z2,v,G,q,D+R,dh,H,Rg)*kn;
q1 = JGv(y0,z1,z2,v,G,q,D-R,dh,H,Rg)*kn;
m = -(q2-q1)*R^2/3;
Fo1=@(z)(1-z/l)*m;Fo2=@(z)z/l*m;
f3=integral(Fo1,0,l);f6=integral(Fo2,0,l);
fe=[0;f1;f3;0;f4;f6];
end

% vertical soil displacement calculated by Mindlin
function ov = JGv(y0,z1,z2,v,G,q,D,dh,H,R)
c = @(z)z;
r1 = @(o,z)((D - R*sin(o)).^2 + (dh + R*cos(o) - y0).^2 + (H - c(z)).^2).^0.5;
r2 = @(o,z)((D - R*sin(o)).^2 + (dh + R*cos(o) - y0).^2 + (H + c(z)).^2).^0.5;
fun1 = @(o,z) q/16/G/pi/(1-v) * R *(D - R*sin(o)).*( (H-c(z))./(r1(o,z)).^3 + (3-4*v)*(H-c(z))./(r2(o,z)).^3 - ...
    6*H*c(z).*(H+c(z))./(r2(o,z)).^5 + 4*(1-v)*(1-2*v)./(r2(o,z).*(r2(o,z)+H+c(z))));
ov = -integral2(fun1,0,2*pi,z1,z2);
end

function qq=Node_FAssembly(fq,fe,i,j)
DOF(1)=3*i-2;DOF(2)=3*i-1;DOF(3)=3*i;DOF(4)=3*j-2;DOF(5)=3*j-1;DOF(6)=3*j;
for n1=1:6
    fq(DOF(n1))=fq(DOF(n1))+fe(n1);
end
qq=fq;
end