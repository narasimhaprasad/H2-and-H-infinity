function [P,err_h2care]=h2care(A,B,C,D,tol)

%H2CARE  Solution to H2 Continuous-time Algebraic Riccati Equation
%
%        P = h2care(A,B,C,D)
%
%        returns a positive semi-definite solution, if existent, for
%        the following algebraic Riccati equation for continuous-time
%        H2 optimal control:
%
%            0 = PA + A'P + C'C - (PB+C'D)(D'D)^{-1}(PB+C'D)'
%
%        Note that a positive semi-definite stabilizing solution is
%        existent if and only if the quadruple (A,B,C,D) has no
%        invariant zeros on the jw axis and D is of full column rank.
%
%        See also H2STATE, H8CARE and H2DARE.

if nargin==4
   tol=1e-8;
end

if rank(D)~=size(D,2)
   P=[];err=Inf;
   disp('D is not full column rank...')
   return
end

R=inv(D'*D);
At=A-B*R*D'*C;
Bt=B*R*B';
Bt=Bt/2+Bt'/2;
Ct=C'*C-C'*D*R*D'*C;
Ct=Ct/2+Ct'/2;

P=are(At,Bt,Ct);

%verify
err_h2care=norm(P*A+A'*P+C'*C-(P*B+C'*D)*inv(D'*D)*(P*B+C'*D)');
