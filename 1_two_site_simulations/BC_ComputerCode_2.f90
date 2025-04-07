
    
!   Program for Properties of Equilibrium Multi-site Gene Model with no epistasis
!   This program uses the approximations in the last section of the Appendix of the LD paper
!   to calculate the properties of the gene model with multiple sites.

    program geneq
!   Solution for equilibrium q and population genetic statistics for the multi-site gene model

    CHARACTER*20 FOUT
    CHARACTER*10 CONT
    
    write(*,*) 'Output file?'
    read(*,*)  FOUT
    OPEN (1,FILE=FOUT)

    write(1,*) 'Solution for equilibrium q for gene model'
    write(1,*) 'Calculates population genetic statistics using solution'
    write(*,*) 'Number of selected sites?'
    read(*,*) G
    write(1,*) 'Number of selected sites= ',G
    write(1,*) ''
    write(*,*) 'Initial q?'
    read(*,*) q0
    write(1,*) 'Initial q= ',q0
    write(*,*) 'Maximum number of iterations? '
    read(*,*) NMAX
    write(1,*) 'Maximum number of iterations= ',NMAX
    write(*,*) 'Selection coefficient?'
    read(*,*) s
    write(1,*) 'Selection coefficient= ',s
    write(*,*) 'Mutation rate?'
    read(*,*) u
    write(1,*) 'Mutation rate= ',u

    write(*,*) 'Continue? Type 0 if yes'
    read (*,*) CONT
!  allows program to be cancelled if desired
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
10  write(*,*) 'Dominance coefficient?'
    read(*,*) h
    if(h.gt.1) go to 50
    write(1,*) ''
    write(1,*) 'Dominance coefficient= ',h

    write(1,*) 'Newton-Raphson iteration of q'
    x=G*q0
    a=2*G*u/s

    do 20 i=1,NMAX
    ex=exp(-x)
    fx=a-x*(1+(x-1)*(1-2*h)*ex)
    dfx=0.0-(1+(x-1)*(1-2*h)*ex)-x*(1+(2-x)*ex*(1-2*h))
    x1=x-fx/dfx
    write(*,*) 'i=',i,' x1=',x1,'fx=',fx
    xe=(x1-x)/x
    if(abs(xe).le.0.0001) go to 30
    if(i.eq.NMAX) go to 40
    x=x1
20  continue

40  write(1,*) 'Iteration does not converge'

30  write(1,*) 'Iteration= ',i
    write(1,*) 'fx= ',fx
    q=x/G
    write(1,*) 'Equilibrium q= ',q
    write(1,*) 'Mean number of mutations per haplotype= ',x
    write(1,*) ''
    write(1,*) 'LD statistics'
    R=q*x*ex*(0.5-h)
    R=R/(1+(x-1)*ex*(1-2*h))
    D=R*q
    write(1,*) 'R= ',R,' D= ',D
    write(1,*) ''
    AL=(x-(1-2*h)*x*ex)*s
    B=x*s-AL
    V=0.5*x*(s*(1+(x-1)*ex*(1-2*h)))**2
    write(1,*) 'Load statistics'
    write(1,*) 'Genetic load= ',AL
    write(1,*) 'Inbreeding load= ',B
    write(1,*) 'Fitness Variance= ',V
    write(1,*) ''
    go to 10

50  end program geneq
