function [ x,error ] = GMRES_Projection( A,b,x_0,m,tol )
%This function does 1 projection with constant Krylov subspace dimension 'm'
%With monitoring of the inner residual

[n,mA]=size(A);
assert(n==mA,'Input matrix must be square');

%Initialzations
V=zeros(n,m+1);
W=zeros(n,m);
H=zeros(m+1,m);


r_0=b-A*x_0;
betta=norm(r_0);
V(:,1)= r_0/betta;
inner_residual=zeros(m,1);

for j=1:m
    
    W(:,j)=A*V(:,j);
    
    for i=1:j
        %Calculate the dot product between the vectors (angle)
        H(i,j)=dot(W(:,j),V(:,i));
        %Calculate the projection vector (angle * vector)
        %Subtract the projection from the vector we want to orthonormalize
        W(:,j)=W(:,j)-H(i,j)*V(:,i);
    end
    H(j+1,j)=norm(W(:,j));
    
    %Tracking the inner residual
    e_1=zeros(j+1,1);
    e_1(1)=1;

    [Q,~]=qr(H(1:j+1,1:j));
    g=ctranspose(Q)*e_1*betta;
    
    inner_residual(j)=abs(g(j+1));
    display(['Inner residual at m = ' num2str(j) ' is ' num2str(inner_residual(j))]);
    
    if (inner_residual(j)<=tol)
        break
    end

    if (H(j+1,j)>1e-10)
        V(:,j+1)=W(:,j)/H(j+1,j);
    else
        break
    end
end
e_1=zeros(j+1,1);
e_1(1)=1;

[Q,R]=qr(H(1:j+1,1:j));
g=ctranspose(Q)*e_1*betta;
y=R\g;
x=x_0+V(:,1:j)*y;
error=abs(g(j+1));

end