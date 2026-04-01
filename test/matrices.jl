
function MBO(Ol,pos,l)
    # many body operator
    # pos and Ol should be sorted by index previously
    I = [1 0; 0 1]
    O = I
    k = 1
    if in(1,pos)
        O = Ol[1]
        k = k+1
    end
    
    for i in 2:l
        if in(i,pos)
            O = kron(O,Ol[k])
            k = k+1
        else
            O = kron(O,I)
        end
    end
    return O
end
