function [qubits, etaShrinkAN, etaShrinkNUM ] = uniform_qubit_states(n)
    %Creates a set of 2n^2 states
    
    vectors=sphere_approximation(n);
        total_vectors=size(vectors,2);
    for i=1:total_vectors
        qubits(:,:,i)=1/2*[1+vectors(3,i), vectors(1,i)-sqrt(-1)*vectors(2,i); vectors(1,i)+sqrt(-1)*vectors(2,i), 1-vectors(3,i)];
    end
    
    etaShrinkAN=cos(pi/(2*n))^2;
    [A,b,Aeq,beq]=vert2lcon(vectors');
    etaShrinkNUM=max(b);
    n_Vertices=size(qubits);

end