function [ x,error,n,j ] = GMRES_Iterative( A,b,x_0,m,r_tol,N )

%Inpputs:
%           A,b: System matrix and RHS such that A*x=b
%           x_0: initial guess
%           m: Krylov subspace dimension
%           r_tol: Accepted Relative Residual error
%           N: maximum number of projections (outer iterations)

%Outputs:
%           x: Final solution obtained
%           error: Final relative residual error
%           n: number of projections (outer iterations)
%           j: dimesion of last krylov subspace (inner iterations)




%Get the size of the matrix
[nA,mA]=size(A);
assert(nA==mA,'Input matrix must be square');

%Display 
display([ 'RUNNING GMRES for (' num2str(nA) 'x' num2str(nA) ') Matrix'])
display(['r_tol = ' num2str(r_tol)]);
display(['N = ' num2str(N)]);
display(['m = ' num2str(m)]);

%Calculate the tol for the residuals based on r_tol
tol=r_tol*norm(b);

%Initialzations
V=zeros(nA,m+1);
W=zeros(nA,m);
H=zeros(m+1,m);

%Vector of the outer residuals in case we want to plot it
outer_residual=zeros(N,1);

%Set x to the initial guess
x=x_0;

%PERFORM N PROJECTIONS (OUTER ITERATIONS)
for n=1:N
    
    %Calculate the residual and normalize it to get the first vector in
    %Krylov subspace
    r=b-A*x;
    betta=norm(r);
    V(:,1)= r/betta;
    
    %Store the outer residual of iteration n
    outer_residual(n)=betta;
    
    %Check if the outer residual less than tol. 
    %DONT DO THE ITERATION (BREAK)
    if (outer_residual(n)<=tol)
        break
    end
    
    %IF WE DIDN'T BREAK DO THE ITERATION
    display(['Outer residual at Outer Iteration = ' num2str(n) ' is ' num2str(outer_residual(n))]);
    
    %Vector for the inner residual in case we want to plot it
    inner_residual=zeros(m,1);
    
    %The next vectors in Krylov subspace
    for j=1:m
        
        %Calculate the next vector as based on Krylo subspace definition
        W(:,j)=A*V(:,j);
        
        %Orthonormalize it
        for i=1:j
            %Calculate the dot product between the vectors (angle)
            H(i,j)=dot(W(:,j),V(:,i));
            %Calculate the projection vector (angle * vector)
            %Subtract the projection from the vector we want to orthonormalize
            W(:,j)=W(:,j)-H(i,j)*V(:,i);
        end
        H(j+1,j)=norm(W(:,j));
        %Check if the norm approximatly equal zero 
        %stop becaues V(:,j) and A*V(:,j) are linearly dependent
        if (H(j+1,j)>1e-10)
            V(:,j+1)=W(:,j)/H(j+1,j);
        else
            %don`t include the last vector 
            %use only first j vectors of V
            j=j-1;
            display(['Break at m = ' num2str(j) ' Exact solution ']);
            break
        end
    
        %TRACKING THE INNER RESIDUAL
        
        %Calculate g as defined
        e_1=zeros(j+1,1);
        e_1(1)=1;
        %[Q,~]=qr(H(1:j+1,1:j));
        %Calculate Full QR decomposition using Givens Rotations
        [Q,~]=Givens_QR_Hessn(H(1:j+1,1:j));
        g=ctranspose(Q)*e_1*betta;
        %inner residual last element in g
        inner_residual(j)=abs(g(j+1));
        display(['Inner residual at m = ' num2str(j) ' is ' num2str(inner_residual(j))]);
        
        %Check if the inner residual less than tol. 
        %DONT DO ADD MORE VECTORS FROM KRYLOV SUBSPACE (BREAK)
        if (inner_residual(j)<=tol)
            display(['Break at m = ' num2str(j) ' Inner residual < tol ']);
            break
        end

    end
    
    %UPDATE 'x'
    
    %Calculate y as defined
    e_1=zeros(j+1,1);
    e_1(1)=1;
    %[Q,R]=qr(H(1:j+1,1:j));
    %Calculate the full QR decomposition using Givens rotation
    [Q,R]=Givens_QR_Hessn(H(1:j+1,1:j));
    g=ctranspose(Q)*e_1*betta;
    %Solve for y using backward substitution 
    y=R\g;
    %Calculate new 'x'
    x=x+V(:,1:j)*y;

end

%If j~=m it means that we breakend in an inner iteration and we will not do
%the next outer iteration, therefore decrement the number of outer
%iterations
if(j~=m)
    n=n-1;
end

%If n=0
%Initial Guess is good w.r.t tol (we didn't do any iterations)
if(n~=0)
    %calculate the relative residual error
    error=abs(g(j+1))/norm(b);
else
    %Calculate the relative residual error
    error=betta/norm(b);
    j=0;
end

end