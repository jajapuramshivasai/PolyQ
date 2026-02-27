"""
̈-Theta Reduction sub Routine

A(O) = sigma { Œ^P(I,M,O) }

I:input space
M:intermediate space
O:output space

Œ: exp(i.2pi/n)
n: base of ring polynomials

Typically I = 0
P(I,M,O) = f(I,M,O) + g(M)

Compute g(M) and cache it as hash map G: M -> g(M) #  Ø(K.2^|M|) // parallalizable

when required calculate f(I,M',O) + G(M) # Ø(K.2^|M'+O|)) + Ø(1) // parallalizable

net complexity reduction: Ø(K.2^|M|) + Ø(K.2^|M'+O|)) + Ø(1) << Ø(K.2^|I+M+O|))

"""