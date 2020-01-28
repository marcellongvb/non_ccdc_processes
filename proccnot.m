D = 0;

CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
% phi = 1/2 * (ketbra(0,0,0,0) + ketbra(1,1,1,1) + ketbra(0,0,1,1) + ketbra(1,1,0,0));

% Removing aux
for i=0:1
    for j=0:1
        for k=0:1
            for l=0:1
             D = D + kron2(ketbra(i,j,k,l),PartialTrace(CNOT*ketbra(i,j,k,l)*CNOT,2,[2 2]));    
            end
        end
    end
end

% for i=0:1
%     for j=0:1
%         for k=0:1
%             for l=0:1
%              D = D + kron2(ketbra(i,j,k,l),CNOT*ketbra(i,j,k,l)*CNOT);    
%             end
%         end
%     end
% end

clear i j k l CNOT
