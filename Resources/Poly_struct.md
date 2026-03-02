# methods in circuit phase polynomial evaluations

serial routine

F(inputs,outputs,x) -> apply inputs -> seperate output poly f(outputs,x) and g(x) -> write g(x) as g1(x1,x) + g2(x2,x) +... for pivot ie we can evaluate cache G, G(x) + gi(xi=0,x) ,G(x) + gi(xi=1,x) -> g(x) = {Z8}g(x) + (Z4)g(x) + {Z2}g(x) send to clifford engine and t engine seperatly -> bilinear form and S reduction -> evaluate at specific output O -> pre define out vec and iter gray coe/pivot distributed along threads -> collect