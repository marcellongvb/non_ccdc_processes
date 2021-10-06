function [Bloch] = Rho2Bloch(Rho)

Bloch=NaN(3,1);
    X=[0 1;1 0];
Y=[0 -sqrt(-1);sqrt(-1) 0];
Z=[1 0;0 -1];
    Bloch(1)= trace(X*Rho);
    Bloch(2)= trace(Y*Rho);
    Bloch(3)= trace(Z*Rho);
end