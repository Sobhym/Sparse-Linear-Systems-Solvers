function [ inner_product, error] = MOA_resampling( A,b,h,k1,k2,N1,N2,method )
% This function calculates the inner product of the solution with a vector
% 'h', by running a large number of random walks on the states of the
% system to estimate the inner product, <h,x_k>, where: 
% - x_k: is the exact iterative solution after (k=k1+k2) iterations  
% - h: is the vector to calculate the inner product with
% Inputs:
%           A,b: System matrix and RHS such that A*x=b
%           h:   Vector we estimate the inner product with.
%           k1:  Number of steps before resampling
%           k2:  Number of steps after resampling
%           N1:  Number of random walks before resampling
%           N2:  Number of random walks after resampling
%           method: (choice depends on the sparsity of 'h' and 'b')
%               - B_F: Backward sampling then Forward sampling
%               - B_B: Backward sampling then Backward sampling
%               - F_B: Forward sampling then Backward sampling
%               - F_F: Forward sampling then Forward sampling
% Outputs:
%           Inner_product: Estimate of the inner product <h,x_k>
%           Error:         Estimate of the probabilistic error bound.

switch upper(method)
    case 'B_F'
        [inner_product_1,error_1,hBk]=Backward_sampling(A,b,h,k1,N1);   
        [inner_product_2,error_2]=Forward_sampling(A,b,hBk,k2,N2);
    case 'B_B'
        [inner_product_1,error_1,hBk]=Backward_sampling(A,b,h,k1,N1);   
        [inner_product_2,error_2]=Backward_sampling(A,b,hBk,k2,N2);
    case 'F_B'
        [inner_product_1,error_1,bBk]=Forward_sampling(A,b,h,k1,N1);   
        [inner_product_2,error_2]=Backward_sampling(A,bBk,h,k2,N2);
    case 'F_F'
        [inner_product_1,error_1,bBk]=Forward_sampling(A,b,h,k1,N1);   
        [inner_product_2,error_2]=Forward_sampling(A,bBk,h,k2,N2);
end
        
inner_product=inner_product_1+inner_product_2;
error=error_1+error_2;
        
end

