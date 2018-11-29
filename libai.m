function objval = libai(chrom)

global sequence
n = length(sequence);
dim = 2*n - 5;  
array = zeros(1,n);

array(1,2:(n+1)) = chrom(1,1:n);
beta(1,3:(n-1)) = chrom(1,(n-1):(2*n-5));


v1 = 0;
v2 = 0;
for x = 2:(n-1)
    v1 = v1 + (1 - cosd(array(x)))./4;
end

position(1,1:3) = [0,0,0];
position(2,1:3) = [0,1,0];
position(3,1:3) = [cosd(array(1,2)),1 + sind(array(1,2)),0];
for j = 4:n
    ax = position((j-1),1);
    ay = position((j-1),2);
    az = position((j-1),3);
    position(j,1:3) = [ax + cosd(array(j-1)).* cosd(beta(j-1)), ay + sind(array(j-1)).* cosd(beta(j-1)), az + sind(beta(j-1))];
end


for i=1:(n-2)
    for j=(i+2):n
        aai = sequence(i);
        aaj = sequence(j);
        if ((aai==aaj))
            if (aai==1)
                ccc=1;
            else
                ccc=0.5;
            end
        else
            ccc=-0.5;
        end
        
        
        % º∆À„r

        r = norm(position(j,1:3) - position(i,1:3));
        %
        v2 = v2 + 4*((r^(-12)) - (r^(-6))*ccc);
        
    end
end
objval =v1+v2;


