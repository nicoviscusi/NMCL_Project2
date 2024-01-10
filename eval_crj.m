function crj = eval_crj(k)

crj = zeros(k+1,k);

for r=-1:k-1  
    for j=0:k-1
        for m=j+1:k
            p = 1.0;
            for l=0:k
                if(l~=m)
                    p = p*(m-l);
                end
            end
            s = 0.0;
            for l=0:k
                p1 = 1.0;
                for q=0:k
                    if(q~=m && q~=l)
                        p1 = p1*(r-q+1.0);
                    end
                end
                if(l~=m)
                    s = s + p1;
                end
            end
            crj(r+2,j+1) = crj(r+2,j+1) + s/p;
        end
    end
end
            