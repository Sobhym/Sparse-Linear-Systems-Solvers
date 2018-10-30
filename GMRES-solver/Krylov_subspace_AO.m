function [ V,H ] = Krylov_subspace_AO( A,r,m )
%This function Gets the Krylov subspace based on Modified Gram-Schmidt

%Inputs     A: system matrix  (n*n)
%           r: vector (n*1) to construct krylov subspace
%           m: number of orthonormal vectors to be constructed other than r
%               total vectors m+1
%Outputs:
%           V: m+1 orthnonormal vectors (n * m+1)
%           H: Hessenberg matrix  (m+1 * m)

[n,~]=size(A);

%Initialzations
V=zeros(n,m+1);
W=zeros(n,m);
H=zeros(m+1,m);

%First vector is r normalized
betta=norm(r); 
V(:,1)= r/betta;

%Construct the next m vectors
for j=1:m
    
    %Vector based on Krylo subspace definition
    W(:,j)=A*V(:,j);
    for i=1:j
        %Calculate the dot product between the vectors (angle)
        H(i,j)=dot(V(:,i),W(:,j));
        %Calculate the projection vector (angle * vector)
        %Subtract the projection from the vector we want to orthonormalize
        W(:,j)=W(:,j)-H(i,j)*V(:,i);
    end
    %Normalize the j-th vector
    H(j+1,j)=norm(W(:,j));
    
    %Check it the norm ~=0 stop becaues V(:,j) and A*V(:,j) are linearly
    %dependent
    if (H(j+1,j)>1e-10)
        V(:,j+1)=W(:,j)/H(j+1,j);
    else
        %don`t include the last vector 
        %use only first j vectors of V
        %dimesion of Krylov subspace = j
        V=V(:,1:j);
        H=H(1:j,1:j-1);
        break
    end
end


end

