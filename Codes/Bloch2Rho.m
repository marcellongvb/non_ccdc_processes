function [Rho] = Bloch2Rho(BlochVector)

    Rho =1/2*[1+BlochVector(3), BlochVector(1)-sqrt(-1)*BlochVector(2); BlochVector(1)+sqrt(-1)*BlochVector(2), 1-BlochVector(3)];
end