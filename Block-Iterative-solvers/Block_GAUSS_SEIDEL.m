function [ x,F_error,iter_count ] = Block_GAUSS_SEIDEL( A,b,x0,tol,n_iter,n_b )
%Solve "A x = b" for x using block Gauss Seidel iterative method
%Inputs
%x0: Initial Guess to be used in the iterative method
%tol: Maximum error accepted in the soln 
%     (such that the error at each iteration = norm(Ax-b)/norm(b))
%n_iter: Maximum number of iterations
%n_b: Number of rows/columns per block (square blocks)

%Outputs
%x: Final solution
%F_error: Final error
%iter_count: Number of Iterations done to get the solution 

[n,~]=size(A);
x=x0;

%Calculate the error defined as "norm(Ax-b)/norm(b)"
iter_count=1;
error(iter_count)= norm(A*x-b)/norm(b);

%Get the number of blocks based on the block size and the number of rows
assert(n>n_b);
blocks_count=round(n/n_b);

while(error(iter_count)>tol && iter_count<n_iter)
    
    iter_count=iter_count+1;
    %Update x
    x_new=zeros(n,1);
    for block_i=0:blocks_count-1        
        %Get the row Range
        [row_start,row_end]=Get_range(n,blocks_count,block_i);
        %Set part of x to the corresponding part of b
        x_new(row_start:row_end)=b(row_start:row_end);
            
        for block_j=0:block_i-1
            %Get the column range
            [col_start,col_end]=Get_range(n,blocks_count,block_j);
            %Update x by subtracting the matrix vector products with new x
            x_new(row_start:row_end)=x_new(row_start:row_end)-A(row_start:row_end,col_start:col_end)*x_new(col_start:col_end); 
        end
        
        for block_j=block_i+1:blocks_count-1
            %Get the column range
            [col_start,col_end]=Get_range(n,blocks_count,block_j);
            %Update x by subtracting the matrix vector products with old x
            x_new(row_start:row_end)=x_new(row_start:row_end)-A(row_start:row_end,col_start:col_end)*x(col_start:col_end); 
        end
        %Multiply x by the inverse of the diagonal
        x_new(row_start:row_end)=A(row_start:row_end,row_start:row_end)\x_new(row_start:row_end);
    end
    
    x=x_new;
    %Calculate the error defined as "norm(Ax-b)/norm(b
    error(iter_count)= norm(A*x-b)/norm(b);

end

%Get the Final error
F_error=error(iter_count);
%Plot Convergence
plot(0:iter_count-1,log(error)/log(10),'--','color',rand(1,3))
xlabel('Number of iterations');
ylabel('Log_1_0(||Ax-b|| / ||b||)');
title(['Block Gauss Seidel convergence for matrix A (' num2str(n) ' x ' num2str(n) ')']);
end

