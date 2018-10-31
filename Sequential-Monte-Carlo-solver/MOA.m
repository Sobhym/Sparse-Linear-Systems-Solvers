function [ inner_product, error] = MOA( A,b,h,k,N )
% This function calculates the inner product of the solution with a vector
% 'h', by running a large number of random walks on the states of the
% system to estimate the inner product, <h,x_k>, where: 
% - x_k: is the exact iterative solution after k iterations  
% - h: is the vector to calculate the inner product with
% Inputs:
%           A,b: System matrix and RHS such that A*x=b
%           h:   Vector we estimate the inner product with.
%           k:   Number of steps in a random walk (depends on system convergence)
%           N:   Number of random walks
% Outputs:
%           Inner_product: Estimate of the inner product <h,x_k>
%           Error:         Estimate of the probabilistic error bound.


%Calculate "B = I - A" where the system is written "x = B x + b"
[n,m]=size(A);
assert(n==m,'Input matrix must be square');
B=eye(n)-A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Setting for the random walks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create p (Initial probability vector) based on h
p=transpose(h);
p=abs(p);
p=p/sum(p);
%p_cdf is the cumulative distribution function (cdf) of 'p'
%used to sample from 'p' pdf
p_cdf=cumsum([0,p]);

%Create P (Transition matrix) based on B
P=abs(B);
%P_cdf: holds the cdf of every cdf (every row) in P
P_cdf=zeros(n,n+1);
for i=1:n
   total=sum(P(i,:));
   P(i,:)=P(i,:)/total;
   P_cdf(i,:)=cumsum([0,P(i,:)]);
end

%Z_k the random variable to be calculated at the end of every walk
%Z_k(1) value of Z_k calculated at the first run
%Z_k(N) vaule of Z_k calculated at last run
Z_k=zeros(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running N random walks %
%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:N
    
    %%%%%%%%%%%%%%%%%
    % A random walk %
    %%%%%%%%%%%%%%%%%
    
    %Get the initail state based on the 'p' pdf
    %Get a random number [0,1]
    r=rand;
    %Check for r >= 'p_cdf' (the cdf of 'p')
    %The random state is the number of times this condition is true
    S_0=sum(r >= p_cdf);
    
    %Initialize the weight function W_0 =1 
    W_0=1;
    %intialize the summation Z = W_0 * b_state_0
    Z=W_0*b(S_0);
    %Store the current state and weight as the old ones
    S_old=S_0;
    W_old=W_0;
    
    %Walk k steps 
    for i=1:k
    
        %P(S_old,:) is the pdf of the transition from state 'S_old' to any other state
        %The cdf of P(S_old,:) to sample the new state is pre computed in P_cdf(S_old,:)

        %Sample the current state
        r=rand;
        S=sum(r >= P_cdf(S_old,:));
        %Calculate the weight `W` of the current state based on the previous weight `W_old
        W=W_old* B(S_old,S)/P(S_old,S);
        %Update the summation Z
        Z=Z+W*b(S);
        %Store the current state and weight as the old ones
        W_old=W;
        S_old=S;
        
    end
    %This condition will not fail as we know from the begining that the
    %matrix is converging
    assert(Z<1e9,'The answer is not converging, spectral radius of the system must be smaller than one');
    %Calculate the random variable `Z_k` at the end of walk 
    Z_k(j)=h(S_0)/p(S_0)*Z;
    
end

%Calculate <h,x_k+1> defined as the mean of the random variable 'Z_k'
inner_product=mean(Z_k);
error=sqrt(var(Z_k)/N);
end

