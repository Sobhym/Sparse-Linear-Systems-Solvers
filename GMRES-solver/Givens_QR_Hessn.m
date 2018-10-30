function [ Q,R ] = Givens_QR_Hessn( H )
%This function gets the Full QR decomposition of a Real Hessenberg matrix H.
%Inputs:    H: Hessenberg matrix of size (n*m) 

%Outputs:   Q*R=H
%           Q: Orthogonal Martrix of size (n*n)
%           R: Upper triangular matrix of size (n*m)

H2=triu(H,-1);
E1=H2-H;
assert(norm(E1)<eps,'Matrix is not Hessenberg');
[n,m]=size(H);
assert(m<=n,'Columns are linearly dependent (m>n)');

R=H;
Q=eye(n,n);

%Only the first lower subdiagonal to be eleminated 
for j=1:m
    
    %row index (1st lower subdiagonal)
    i=j+1;
    if(i>n)
        break
    end

    %Calculate c,s to do the transformation
    if R(i,j)==0
        c=1;
        s=0;
    elseif R(j,j)==0
        c=0;
        s=R(i,j)/abs(R(i,j));
    else
        c=abs(R(j,j))/sqrt((R(j,j)^2) + (R(i,j)^2));
        s=R(j,j)/abs(R(j,j));
        s=s*R(i,j)/sqrt((R(j,j)^2) +(R(i,j)^2));    
    end
    
    %Update R and Q by transforming them based on the computed values of c,s
    % Efficient calculation of R=O*R where O is the 
    % transformation matrix 
    row_j=c*R(j,:)+s*R(i,:);
    row_i=-s*R(j,:)+c*R(i,:);
    R(j,:)=row_j;
    R(i,:)=row_i;
    % Efficient calculation of Q=Q*transpose(O) where O is the
    % transformation matrix
    col_j=c*Q(:,j)+s*Q(:,i);
    col_i=-s*Q(:,j)+c*Q(:,i);
    Q(:,j)=col_j;
    Q(:,i)=col_i;
    
end


end



