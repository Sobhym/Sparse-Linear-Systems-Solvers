function [ X_full] = MCMC_full_soln( A,b,k,N )
%This function calculates 'x_k+1' the full solution vector
%Run every state as the initial state N walks to calculate X_full(j) 
%defined as the mean of the random variable 'Z_k(j)' j=1,2,..n
%where     -one sample of the random varible 'Z_k(j)' is calculated at the end 
%            of walk ‘i’ after walking k steps on the states of the matrix i=1,2,..N            
%          -'x_k+1' is iterative solution at iteration k+1

%A,b are the system and RHS respectively such that "A x = b"
%k is the number of steps in a random walk 
%N is the number of random walks

%Calculate "B = I - A" where the system is written "x = B x + b"
[n,m]=size(A);
assert(n==m,'Input matrix must be square');
B=eye(n)-A;
%Calculate the eigen values to check if the sequence will converge
%In a practical impelenation the calculation of the eigen values will not
%be included but i added it for illustration. 
Eig_B=abs(eig(B));
Spectral_radius=max(Eig_B);
sprintf([ 'Spectral radius is ' num2str(Spectral_radius) ]) 
assert(Spectral_radius < 1,'Spectral radius bigger than 1 (non-converging)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Setting for the random walks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

%The full solution
X_full=zeros(n,1);

%Run every state as the initial state N times to calculate X_full(j) 
%defined as the mean of the random variable 'Z_k(j)' j=1,2,..n
for S_0=1:n
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running N random walks %
%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:N
    
    %%%%%%%%%%%%%%%%%
    % A random walk %
    %%%%%%%%%%%%%%%%%
    
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
    assert(Z<1e9,'The answer is not converging');
    %Calculate the random variable `Z_k` at the end of walk 
    Z_k(j)=Z;
    
end

%Calculate <h,x_k+1> defined as the mean of the random variable 'Z_k'
X_full(S_0)=mean(Z_k);
end

end

