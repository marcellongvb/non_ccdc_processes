function [Vout] = RemoveEquivalentVectors(V)


cont=0; 

number_vertices=size(V,2);
for i=1:number_vertices
    for j=1:number_vertices
        if i==j
        elseif norm(V(:,i)-V(:,j))<10^(-8)
            V(:,i)=[];
        end
    end
end

end