function [ PI ] = CSR_BFS( IA,JA,i )
% Create a permutation set based on Breadth First Search (BFS) method for matrices stored in Compressed Sparse Row format (CSR)

%Inputs
% IA, JA: Input matrix in CSR format
% i: initial node
%Output
% PI: Permutation set


%Number of rows (nodes)
n=length(IA)-1;
%To track the visted nodes (rows)
marked=zeros(n,1);
%PI Reordered index set
PI=zeros(n,1);
%Set the initial level Set  (only contain one node)
S(1)=i;
%add initial node to the set PI
PI(1)=i;
%Mark the node as visted
marked(i)=1;
%Set the number of visted nodes (only one node visted yet)
n_visited=1;

%Check if there are nodes not visited yet 
while n_visited<n
   
    %Initialize the new set (S_new)
    S_new=zeros(0,0); 
    
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
                %increase the number of nodes visited
                n_visited=n_visited+1;
                %Add it to the set PI
                PI(n_visited)=JA(j);
                
            end
        end
    end
    
    %If the new set is empty
    if(isempty(S_new))
        %find unmarked node
        i=find(marked==0,1);
        %mark the node
        marked(i)=1;
        %increase number of nodes visted
        n_visited=n_visited+1;
        %add it to PI
        PI(n_visited)=i;
        %Set S_new
        S_new(end+1)=i;
    end
    S=S_new;
end

end

