function [ start,stop ] = Get_range( n,blocks_count,i )


divisor=floor(n/blocks_count);
remainder=rem(n,blocks_count);

if (i<remainder)
    offset=i;
else
    offset=remainder;
end

start=offset+i*divisor+1;
stop=start+divisor-1;
if(remainder>i)
    stop=stop+1;
end

end

