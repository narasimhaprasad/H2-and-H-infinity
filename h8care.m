function [P,err,flag]=h8care(A,B,C,D,E,gamma,tol)

%H8CARE  Solution to H-infinity Continuous-time Algebraic Riccati Equation
%
%        P = h8care(A,B,C,D,E,gamma)
%
%        returns a positive semi-definite solution, if existent, for
%        the following algebraic Riccati equation for continuous-time
%        H-infinity control:
%
%                               [ B'P+D'C ]'        [ B'P+D'C ]
%          0 = PA + A'P + C'C - |         |  G^{-1} |         |
%                               [   E'P   ]         [   E'P   ]
%
%        where
%                       [ D'D      0      ]
%                   G = |                 |
%                       [  0   -gamma^2 I ]
%
%        This CARE is related to H-infinity control for the following
%        continuous-time system:
%               .
%               x = A x + B u + E w
%               h = C x + D u
%
%        Note that a positive semi-definite stabilizing solution is
%        existent if and only if the quadruple (A,B,C,D) has no
%        invariant zeros on the jw axis, D is of full column rank,
%        and gamma > gamma_\infty^*, which can be determined using
%        GM8STAR.
%
%        See also H2CARE, H8STATE, H8DARE and GM8STAR.
%

if nargin==6
   tol=1e-8;
end

[n,q]=size(E);
t1=[B';E'];t2=[D'*C;zeros(q,n)];

if rank(D)<size(D,2)
   P=[];err=inf;flag=0;
   return
end
invG=blkdiag(inv(D'*D),-1/gamma/gamma*eye(q));

At=A-t1'*invG*t2;
Bt=t1'*invG*t1;
Ct=C'*C-t2'*invG*t2;

P=are(At,Bt,Ct);
flag=(min(real(eig(P)))>-tol);

%verify
err=norm(P*A+A'*P+C'*C-[B'*P+D'*C;E'*P]'*invG*[B'*P+D'*C;E'*P]);

%if flag==0
%   P=[];
%end