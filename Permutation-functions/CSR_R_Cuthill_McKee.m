function [ PI ] = CSR_R_Cuthill_McKee( IA,JA )
% Create a permutation set based on Reverse Cuthill-McKee method for matrices stored in Compressed Sparse Row format (CSR)

%Inputs
% IA, JA: Input matrix in CSR format
%Output
% PI: Permutation set


%Number of rows (nodes)
n=length(IA)-1;
%To track the visted nodes (rows)
marked=zeros(n,1);
%PI Reordered index set
PI=zeros(n,1);
%number of elements in a row = connections = degree of the node
n_row=zeros(n,1);
for i=1:n
    n_row(i)=IA(i+1)-IA(i);
end
%Highest degree over all nodes
max_n_row=max(n_row);

%Find initial node (row) with least connections (degree)
[~,i]=min(n_row);
%Update the node degree to more than the max 
%(not to showup again when getting the min as this node is already visited)
n_row(i)=max_n_row+1;

%Set the initial level Set  (only contain one node)
S=zeros(0,0);
S(1)=i;
%Mark the node as visted
marked(i)=1;
%Set the number of visted nodes (only one node visted yet)
n_visited=1;

%add initial node to the set PI
PI(1)=i;

%Check if there are nodes not visited yet 
while n_visited<n
   
    %Initialize the new set (S_new)
    S_new=zeros(0,0);
    %Intialize the degrees of the new sets
    %S_n_row to keep the track of the degrees of the nodes in S_new
    S_n_row=zeros(0,0);
    
    %For all the nodes in the previous set (S)
    for i=1:length(S)
        
        %Check the nodes connected to every node in the set (S(i))
        for j=IA(S(i)):IA(S(i)+1)-1
            
            %If conected to unmarked node
            if (marked(JA(j))==0) 
                %Mark the node
                marked(JA(j))=1;
                %Add it to the new set
                S_new(end+1)=JA(j);
                %Get the degree of the node
                S_n_row(end+1)=n_row(JA(j));
                %Set the degree to be more than max
                assert(n_row(JA(j))~= (max_n_row+1))
                n_row(JA(j))=max_n_row+1;

            end
        end
    end

    %Sort S_new ascendingly based on the degree
    [~,order]=sort(S_n_row);
    S_new=S_new(order);

    
    %Add Ordered S_new to PI
    for i=1:length(S_new)
        %increase the number of nodes visited
        n_visited=n_visited+1;
        %Add it to the set PI
        PI(n_visited)=S_new(i);
    end
                
    %If the new set is empty
    if(isempty(S_new))
        
        %Get a new node with least degree
        [~,i]=min(n_row);
        %Set the degree to be more than max
        assert(n_row(i)~= (max_n_row+1));
        n_row(i)= (max_n_row+1);
        assert(marked(i)==0,num2str(i))

        %Mark the node
        marked(i)=1;
        %Add it to S_new set
        S_new(end+1)=i;
        %increase the number of nodes visited
        n_visited=n_visited+1;
        %Add it to the set PI
        PI(n_visited)=i;
    end
    
    S=S_new;
end
PI=flipud(PI);
end