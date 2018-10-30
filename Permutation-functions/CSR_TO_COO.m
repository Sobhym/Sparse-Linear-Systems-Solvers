function [ AS ] = CSR_TO_COO( IA,JA,VA,m )
% Convert a Matrix stored in Compressed Sparse Row (CSR) format to Coordinate List (COO) format

%Inputs
% IA, JA: Indices
% IV: values
% m: number of columns 
%Outputs
% AS: Matrix stored in COO format


cols=JA;
values=VA;
n=length(IA)-1;

rows=zeros(length(values),1);



for i=1:n
    for k=IA(i):IA(i+1)-1
    
        assert(rows(k)==0)
        rows(k)=i;
    end
end

AS=sparse(rows,cols,values,n,m);
end

