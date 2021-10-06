function [V,points] = sphere_approximation(n)


points = zeros(3,2*n^2);

total_vertices=0;
for i1=1:n
    for i2=1:n
        total_vertices=total_vertices+1;
        points(:,total_vertices)=[cos(i1*pi/n)*cos(i2*pi/n) sin(i1*pi/n)*cos(i2*pi/n)  sin(i2*pi/n)]; 
        
        total_vertices=total_vertices+1;
        points(:,total_vertices)=[-cos(i1*pi/n)*cos(i2*pi/n) -sin(i1*pi/n)*cos(i2*pi/n)  -sin(i2*pi/n)]; 
    end
end

decimal_precision=7;
V=fix(points*10^decimal_precision)/10^decimal_precision;
V=unique(V','rows')';
end