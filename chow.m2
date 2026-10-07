needsPackage "Resultants"

-- Extract the Stiefel coordinates (row span) of L
rowspan = L -> (
    K = coefficientRing ring L;
    N = sub(jacobian gens trim L, K);
    transpose gens ker transpose N
);

-- Function to compute the reciprocal linear space
reciprocalSpace = L -> (
    R = ring L;
    n = numgens R;
    P = product gens R;
    invImages = apply(toList(0..(n-1)), i -> product drop(gens R, {i,i}));
    invMap = map(R, R, invImages);
    saturate(invMap(L), P)
);

-- Function to extract the Plucker coordinates as a list
pluckerListFromStiefel = M -> (
    d = min(numgens source M, numgens target M);
    pluckerMatrix = exteriorPower(d, M);
    flatten entries pluckerMatrix 
);

-- Function to calculate the LTDR for a reciprocal linear space
ltdrOfReciprocal = L -> (
    R = ring L;
    K = coefficientRing R;
    n = dim R;
    M = rowspan(L);
    d = min(numgens source M, numgens target M);
    k = binomial(n-1, d-1);
    V = K^k;
    
    -- Generate the coordinate ring of the Grassmannian
    G = Grass(d-1, n-1, K);
    pLList = pluckerListFromStiefel M;

    -- Index lists to convert Plucker indeces to the iterator
    pIndeces = subsets(toList(0..(n-1)), d);
    vIndeces = subsets(toList(0..(n-2)), d-1);

    ltdr = matrix mutableMatrix(G, k, k);
    
    -- Check if the uniform case applies (no zero Plucker coordinates)
    if not member(0_K, pLList) then (
        for j from 0 to binomial(n,d) -1 do (
            -- Check if n is in I (Macaulay2 is 0-indexed, so we check for n-1)
            if member(n-1, pIndeces#j) then (
                vI = V_(j - binomial(n-1, d));
            )
            else (
                -- Case where n is not in I
                I = pIndeces#j;
                vI = 0_V;
                for l in I do(
                    Iwithoutl = delete(l, I);
                    -- Find the position of I \setminus \{l\} in the index list
                    IListIndex = position(vIndeces, x -> x==Iwithoutl);   

                    -- Add the alternating sum components                                     
                    vI = vI + (-1)^(position(I, x -> x==l)) * V_IListIndex;                                  
                );
            );
            ltdr = ltdr + (1 / pLList#j) * (gens G)#j * matrix vI * transpose matrix vI;
        );
    );
    ltdr
);

-- Main function to compute the Chow form using the LTDR
chowFormOfReciprocal = L -> (
    R = ring L;
    K = coefficientRing R;
    n = dim R;
    M = rowspan(L);
    d = min(numgens source M, numgens target M);
    k = binomial(n-1, d-1);
    
    G = Grass(d-1, n-1, K);
    
    pIndeces = subsets(toList(0..(n-1)), d);
    pIndecesDual = subsets(toList(0..(n-1)), n-d);

    ltdr = ltdrOfReciprocal(L);
    GDual = Grass(n-d-1, n-1, K);

    imagesList = {};
    -- Iterate through all Plucker coordinates to construct the wedge product
    for j from 0 to binomial(n,d)-1 do (
        jDual = toList(0..(n-1)) - set pIndeces#j;
        sI = (-1)^((1/2)*d*(d-1) + sum(pIndeces#j));
        betaI = (gens GDual)#(position(pIndecesDual, x -> x==jDual));
        imagesList = append(imagesList, sI * betaI);
    );
    wedgeMap = map(GDual, G, imagesList);
    chowMatrix = wedgeMap(ltdr);
    det(chowMatrix)
);

K = toField(QQ[i]/ideal(i^2+1));
R = K[x_0..x_3];

R = K[x_0, x_1, x_2];
L1 = ideal(x_0-2*x_1+x_2);

M1 = matrix(R, {{1,0}, {1,1}, {1,2}, {1,3}});
L2 = image M1;

chowForm(reciprocalSpace(L1))
chowFormOfReciprocal(L1)