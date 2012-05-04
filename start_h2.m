clc
clear
close all

A = [-0.1778 0 0 0 0 -9.7807 -9.7807 0 0 0 0;
      0 -0.3104 0 0 9.7807 0 0 9.7807 0 0 0;
      -0.3326 -0.5353 0 0 0 0 75.7640 343.8600 0 0 0;
      -0.1903 -0.2490 0 0 0 0 172.6200 -59.9580 0 0 0;
       0 0 1.0000 0 0 0 0 0 0 0 0;
       0 0 0 1.0000 0 0 0 0 0 0 0;
       0 0 0 -1.0000 0 0 -8.1222 4.6535 0 0 0;
       0 0 -1.0000 0 0 0 -0.0921 -8.1222 0 0 0;
       0 0 0 0 0 0 0 0 -0.6821 -0.0535 0;
       0 0 0 0 0 0 0 0 -0.2892 -5.5561 -36.6740;
       0 0 0 0 0 0 0 0 0 2.7492 -11.1120];
B = [zeros(6,4);
     0.0496 2.6224 0 0;
     2.4928 0.1741 0 0;
     0 0 7.8246 0;
     0 0 1.6349 -58.4053;
     0 0 0 0];
C = [eye(6,11); zeros(2,8) [1 0 0; 0 1 0]];
D = zeros(8,4);
E =[A(:,1),A(:,2),A(:,9)];

x0 = [1 0 0 0 0 -0.1 0 0 0 0 0]';

epsilon = 1;

t =0:0.1:40;
wind= max(cos(2*pi/40*(t-20)),0);
uwind=10*wind;
vwind=10*wind;
wwind=3*wind;
len=length(t);
wind_dis=[t',uwind',vwind',wwind'];
z=zeros(len,1);
meanoise=[t',0.01*randn(len,4),0.001*randn(len,3),z,0.01*randn(len,2),0.001*randn(len,1)];

Ecap = horzcat(E,[epsilon*eye(7);zeros(4,7)],zeros(11,1));
C1 = C;
D1cap = [zeros(8,2) zeros(8,1) epsilon*eye(8)];
C2cap = [zeros(1,11); epsilon*eye(11); zeros(4,11)];   
D2cap = [zeros(1,4); zeros(11,4); epsilon*eye(4)];   

P = h2care(A,B,C2cap,D2cap);
Q = h2care(A',C1',Ecap',D1cap');

eigp = eig(P)
eigq = eig(Q)

C2 = C2cap;
D1 = D1cap;
D2 = D2cap;

F = -inv((D2'*D2))*(D2'*C2+B'*P);
K = -(Q*C1'+Ecap*D1')/(D1*D1');

Acmp = A+B*F+K*C1;
Bcmp = -K;
Ccmp = F;

gam2star = sqrt(trace(Ecap'*P*Ecap)+trace((A'*P+P*A+C2'*C2)*Q))