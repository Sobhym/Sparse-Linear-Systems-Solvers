function [ IA,JA,VA,m ] = TO_CSR( A )
% Store a Matrix in Compressed Sparse Row (CSR) format

%Inputs
% A: Input Matrix

%Outputs
% IA, JA: Indices
% IV: values
% m: number of columns 

[n,m]=size(A);
%Get the values and indeces (transpose to search by rows)
[cols, rows,values]=find(transpose(A));
JA=cols;
VA=values;

%Count the the number of elements in each row
all_rows=1:n;
[row_count,~]=hist(rows,all_rows);
IA=zeros(n+1,1);

index=1;
for i=1:n
   
    IA(i)=index;
    index=index+row_count(i);
end

IA(n+1)=index;

end

