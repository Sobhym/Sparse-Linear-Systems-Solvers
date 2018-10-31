function [ inner_product,error,bBk ] = Forward_sampling( A,b,h,k,N )

%Calculate "B = I - A" where the system is written "x = B x + b"
[n,m]=size(A);
assert(n==m,'Input matrix must be square');
B=eye(n)-A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Setting for the random walks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create q (Initial probability vector) based on b
q=transpose(b);
q=abs(q);
q=q/sum(q);
%p_cdf is the cumulative distribution function (cdf) of 'p'
%used to sample from 'p' pdf
q_cdf=cumsum([0,q]);

%Create P (Transition matrix) based on B
Q=abs(B);
Q=transpose(Q);
%P_cdf: holds the cdf of every cdf (every row) in P
Q_cdf=zeros(n,n+1);
for i=1:n
   total=sum(Q(i,:));
   Q(i,:)=Q(i,:)/total;
   Q_cdf(i,:)=cumsum([0,Q(i,:)]);
end

%The distribution after k steps
bBk=zeros(n,1);
%Z_k the random variable to be calculated at the end of every walk
Z_k=zeros(N,1);
for j=1:N
    
    %%%%%%%%%%%%%%%%%
    % A random walk %
    %%%%%%%%%%%%%%%%%
    
    %Get the initail state based on the 'q' pdf
    %Get a random number [0,1]
    r=rand;
    %Check for r >= 'q_cdf' (the cdf of 'q')
    %The random state is the number of times this condition is true
    S_0=sum(r >= q_cdf);
    
    %Initialize the weight function W_0 =1 
    W_0=b(S_0)/q(S_0);
    %intialize the summation Z = W_0 * b_state_0
    Z=W_0*h(S_0);
    %Store the current state and weight as the old ones  
    S_old=S_0;
    W_old=W_0;
    
    %Walk k steps from h using Q Transition Matrix
    for i=1:k

        %Sample the current state
        r=rand;
        S=sum(r >= Q_cdf(S_old,:));
        
        %Calculate the weight `W` of the current state based on the previous weight `W_old
        W=W_old* B(S,S_old)/Q(S_old,S);
        %Update the summation Z
        Z=Z+W*h(S);
        %Store the current state and weight as the old ones
        W_old=W;
        S_old=S;
       
    end
    %This condition will not fail as we know from the begining that the
    %matrix is converging
    assert(Z<1e9,'The answer is not converging, spectral radius of the system must be smaller than one');
    %Calculate the random variable `Z_k` at the end of walk 
    Z_k(j)=Z;
    %walk on more step to calc Bkb
    %Sample the current state
    r=rand;
    S=sum(r >= Q_cdf(S_old,:));
    %Calculate the weight `W` of the current state based on the previous weight `W_old
    W=W_old* B(S,S_old)/Q(S_old,S);
    bBk(S)=bBk(S)+W;
end
%normalize Bkb
bBk=bBk/N;
%Calculate <h,x_k+1> defined as the mean of the random variable 'Z_k'
inner_product=mean(Z_k);
error=sqrt(var(Z_k)/N);

end