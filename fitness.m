function objval = fitness(chrom)
global sequence
dim=size(chrom,2);  
array=zeros(1,dim+2);
for x=1:dim
    array(1,x+1)=chrom(x);
end
n=dim+2;
v1=0;
v2=0;
for x=2:(n-1)
    v1=v1+(1-cosd(array(x)))*0.25;
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
        
        zhengxian=0;
        yuxian=0;
        for k=(i+1):(j-1)
            angle=0;
            for l=(i+1):k
                angle=angle + array(l);
            end
            yuxian = yuxian + cosd(angle);
            zhengxian = zhengxian + sind(angle);
        end
        r = ((1 + yuxian)^2 + (zhengxian)^2)^0.5;
        v2 = v2 + 4*(r^(-12)-(r^(-6))*ccc);
        
    end
end
objval =v1+v2;