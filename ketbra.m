function ketbraout = ketbra(i,j,k,l)

ketbraout = kron2(ket(i),ket(j))*kron2(bra(k),bra(l));

end

