!  The programs are written in GFortran for the MacIntosh 
!    (https://gnuc.org.wiki/GFortranBinariesMacOS)
!    Written by Brian Charlesworth. Contact: brian.charlesworth@ed.ac.uk for questions.
!
!    Program for Properties of Deterministic Equilibrium Populations
!    This program uses iteration to calculate the equilibrium haplotype frequencies and load 
!    statistics for the two-locus gene and sites models with no recombination, allowing for 
!     epistasis.

     program dom1b
!   program for 2-locus dominance model
    double precision :: h,u,x,y,p,q,ax,bx,cx,discx,k1,Lg0,Ls0,qr
    double precision :: ay,by,cy,discy,fx,fy,D,R,z,z1,k2,Lg,Ls,Bg,Bs,Vg,Vs
    integer :: it,nmax

    CHARACTER*20 FINP
    CHARACTER*20 FOUT

    WRITE (*,*) 'Input file?'
    READ (*,*) FINP
    OPEN (2,FILE=FINP)
    WRITE (*,*) 'Output file?'
    READ (*,*) FOUT
    OPEN (1,FILE=FOUT)
    write (1,*) 'Equilibria for two-locus dominance models'
    write (1,*) 'No recombination; epistasis allowed'
    write (1,*) ''
    READ (2,*) u
    READ (2,*) h
    READ (2,*) nmax
    READ (2,*) neset
    write (1,*) ''
    write (1,*) 'Ratio of mutation rate to selection coefficient= ',u
    write (1,*) 'Dominance coefficient= ',h
    write (1,*) 'Maximum number of iteration= ',nmax
    write (1,*) 'Number of e values= ',neset
    write (1,*) ''
    write (*,*) 'Ratio of mutation rate to selection coefficient= ',u
    write (*,*) 'Dominance coefficient',h
    write (*,*) 'Maximum number of iteration= ',nmax
    write (*,*) 'Number of e values= ',neset
    
    do 100 i=1,neset
    READ (2,*) e
    write (1,*) ''
    write (1,*) 'Epistasis coefficient= ',e
    write (1,*) ''
    write (1,*) 'Gene model results'
    write (1,*) ''
    write (*,*) 'Gene model results'
    write (*,*) ''
    y=0
    it=0
10  ax=2*(1-h)
    bx=h+3*u+y*(1.5*(1+e)-h)
    cx=0-u*(1-y)
    write(*,*) 'ax=',ax,' bx= ',bx,' cx= ',cx
    discx=(bx**2)-4*ax*cx
    discx=dsqrt(discx)
    x=0.0-bx+discx
    x=x/(2*ax)
    fx=ax*(x**2)+bx*x+cx
    
    ay=(1+e)*(1-h)
    by=(1+e)*(h+x*(1.5-2*h))
    cy=0-u*x
    write(*,*) 'ay=',ay,' by= ',by,' cy= ',cy
    discy=(by**2)-4*ay*cy
    discy=dsqrt(discy)
    y=0-by+discy
    y=y/(2*ay)
    fy=ay*(y**2)+by*y+cy
    it=it+1
    q=x+y
    p=1-q
    D=p*q-x
    R=D/(p*q)

    write (*,*) 'Iteration= ',it
    write (*,*) 'fx= ',fx,' fy= ',fy
    write (*,*) 'x= ',x,' y= ',y
    write (*,*) 'q= ',q,' D= ',D,' Correlation= ',R
    
    if(it.gt.nmax) then
    write (1,*) ''
    write (1,*) 'Exact results'
    write (1,*) 'Iteration= ',it
    write (1,*) 'fx= ',fx,' fy= ',fy
    write (1,*) ''
    write (1,*) 'x= ',x,' y= ',y
    write (1,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (1,*) ''
    qr=q/(u/h)
    write (1,*) 'Relative q= ',qr
    write (1,*) ''

    z=1-2*x-y
!   frequency of wild-type haplotype
    Lg0=4*z*h*(x+y*(1+e))+4*x*(x+1.5*y*(1+e))+2*(1+e)*y**2
    Lg=Lg0/(4*u)
!   gene model load relative to 4 x mutation rate (additive value)
    Bg=(x+y*(1+e))/(2*u)-Lg
    if(h.le.001) go to 50
    Bg=Bg/(1.0/(2*h)-1)
!   inbreeding load relative to additive value
50  write (1,*) 'Relative Lg= ',Lg,' Relative Bg= ',Bg
    write (1,*) ''
    Vg=4*z*x*h**2+4*z*y*(h*(1+e))**2+4*x**2+4*x*y*(1.5*(1+e))**2
    Vg=Vg+(y*2*(1+e))**2-Lg0**2
    Vg=Vg/(4*u)
    write (1,*) 'Relative Vg x h = ',Vg
    write (1,*) ''
    write (1,*) 'Approximate results 1'
    z=h+3*u
    x=u/z
    x=x-2*(1-h)*(u**2)/z**3
    y=0.0-cy/by
    q=x+y
    p=1-q
    D=p*q-x
    R=D/(p*q)
    write (1,*) 'x= ',x,' y= ',y
    write (1,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (1,*) ''
    write (1,*) 'Approximate results 2'
    z1=u/h
    x=z1*(1-2*(1-h)*z1/h)
    k1=(1.5/h)-2
    y=x*z1*(1-k1*x)/(1+e)
    q=x+y
    D=p*q-x
    R=D/(q*(1-q))
    write (1,*) 'x= ',x,' y= ',y
    write (1,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (1,*) ''
    write (1,*) 'Approximate results 3'
    D=(z1**3)*((1.0/(2*h))-e/z1)
    R=D/z1
    write (1,*) 'D= ',D,' Correlation= ',R
    write (1,*) ''
    go to  30
    end if

    go to 10

30  write (1,*) ''
    write (1,*) 'Site model results'
    write (1,*) ''
    write (*,*) 'Site model results'
    write (*,*) ''
    y=0
    it=0
40  ax=1+2*h*e
    bx=h+3*u+y*(1+(1+h)*e)
    cx=0-u*(1-y)
    write(*,*) 'ax=',ax,' bx= ',bx,' cx= ',cx
    discx=(bx**2)-4*ax*cx
    discx=dsqrt(discx)
    x=0-bx+discx
    x=x/(2*ax)
    fx=ax*(x**2)+bx*x+cx
    
    ay=(1+e)*(1-h)
    by=(h+x*(1-h))*(1+e)
    cy=0-u*x
    write(*,*) 'ay=',ay,' by= ',by,' cy= ',cy
    discy=(by**2)-4*ay*cy
    discy=dsqrt(discy)
    y=0-by+discy
    y=y/(2*ay)
    fy=ay*(y**2)+by*y+cy
    it=it+1
    q=x+y
    p=1-q
    D=p*q-x
    R=D/(p*q)
    
    write (*,*) 'Iteration= ',it
    write (*,*) 'fx= ',fx,' fy= ',fy
    write (*,*) 'x= ',x,' y= ',y
    write (*,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (*,*) ''

    if(it.gt.nmax) then
    write (1,*) ''
    write (1,*) 'Iteration= ',it
    write (1,*) 'fx= ',fx,' fy= ',fy
    write (1,*) ''
    write (1,*) 'Exact results'
    write (1,*) 'x= ',x,' y= ',y
    write (1,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (1,*) ''
    qr=q/(u/h)
    write (1,*) 'Relative q= ',qr
    write (1,*) ''

    z=1-2*x-y
!   frequency of wild-type haplotype
    Ls0=4*z*(x*h+y*h*(1+e))+4*h*(1+e)*x**2+2*x**2+4*x*y*(1+h)*(1+e)+2*(1+e)*y**2
    Ls=Ls0/(4*u)
    Bs=(x+y*(1+e))/(2*u)-Ls
    if(h.le.001) go to 55
    Bs=Bs/(1.0/(2*h)-1)
!   inbreeding load relative to additive value
55  write (1,*) 'Relative Ls= ',Ls,' Relative Bs= ',Bs
    write (1,*) ''
    Vs=4*z*x*(h**2)+4*z*y*(h*(1+e))**2+8*(x*h*(1+e))**2+2*x**2
    Vs=Vs+4*x*y*((1+h)*(1+e))**2+(y*2*(1+e))**2-Ls0**2
    Vs=Vs/(4*u)
    write (1,*) 'Relative Vs x h = ',Vs
    write (1,*) ''
    write (1,*) 'Approximate results 1'
    z2=h+3*u
    x=u/z2
    x=x-(u**2)*(1+2*h*e)/z2**3
    y=0.0-cy/by
    q=x+y
    p=1-q
    D=p*q-x
    R=D/(p*q)
    write (1,*) 'x= ',x,' y= ',y
    write (1,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (1,*) ''
    write (1,*) 'Approximate results 2'
    x=z1-(z1**2)*(1+2*h*e)/h
    k2=(1/h)-1
    y=x*z1*(1-k2*x)/(1+e)
    q=x+y
    p=1-q
    D=p*q-x
    R=D/(q*(1-q))
    write (1,*) 'x= ',x,' y= ',y
    write (1,*) 'q= ',q,' D= ',D,' Correlation= ',R
    write (1,*) ''
    write (1,*) 'Approximate results 3'
    D=(z1**3)*(1-e/z1)
    R=D/z1
    write (1,*) 'D= ',D,' Correlation= ',R
    write (1,*) ''
    
    go to  100
    end if

    go to 40

100 continue
    end program dom1b
    
