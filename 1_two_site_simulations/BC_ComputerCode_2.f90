
!     The programs are written in GFortran for the MacIntosh 
!    (https://gnuc.org.wiki/GFortranBinariesMacOS)

    
    program dom5
!   recursions for 2-locus mutation model with multiplicative fitnesse
    double precision :: h,u,x,y,z,p,q,q0,wx,wy,wbar,wz,s,x0,y0,disc,wprod
    double precision :: D,R,x1,y1,z1,u1,Lg,Lg0,Ls,Ls0,Bg,Bs,Vg,Vs,qr,hs
    double precision :: Dapp,Rapp
    integer :: it,nmax,ngen,ind

    CHARACTER*20 FINP
    CHARACTER*20 FOUT

    WRITE (*,*) 'Input file?'
    READ (*,*) FINP
    OPEN (2,FILE=FINP)
    WRITE (*,*) 'Output file?'
    READ (*,*) FOUT
    OPEN (1,FILE=FOUT)
    write (1,*) 'Exact recursions for two-locus mutation model'
    write (1,*) 'Multiplicative fitnesses without epistasis'
    write (1,*) 'No recombination; D=0 initially'
    write (1,*) ''
    READ (2,*) u
    READ (2,*) s
    READ (2,*) ngen

    write (*,*) 'Sites or genes (type 0 or 1)'
    READ (*,*) ind
    
    if(ind.eq.0) then
    write (1,*) 'Sites model'
    else
    write (1,*) 'Gene model'
    end if

    write (1,*) ''
    write (1,*) 'Mutation rate= ',u
    write (1,*) 'Selection coefficient= ',s
    write (1,*) 'Maximum number of generations= ',ngen
    write (1,*) ''
    write (*,*) 'Mutation rate= ',u
    write (*,*) 'Selection coefficient',s
    write (*,*) 'Maximum number of generations= ',ngen
    write (1,*) ''
   
    u1=u/s
    write(*,*) 'Continue?'
    read (*,*) CONT
    if(CONT.eq.1) go to 100

10  write (*,*) 'Dominance coefficient'
    read (*,*) h
    if(h.ge.1) go to 100
    write (1,*) ''
    write (1,*) 'Dominance coefficient',h
    write (*,*) 'Dominance coefficient',h
    hs=h*s
    
    if(h.le.0.001) then
    q0=dsqrt(u1)
    else
    disc=((hs+u)**2)+4*u*s*(1-3*h)
    disc=sqrt(disc)
    q0=(-(hs+u)+disc)/(2*s*(1-3*h))
    if(h.ge.0.499) then
    q0=u/hs
    end if
    end if
!   approximate equilbrium allele frequency with D=0
    p=1-q0
    x0=p*q0
    y0=q0**2
!   starts iteration from approximate D=0 equilibrium
    write (1,*) 'Initial state (approximate equilibrium equation; no LD)'
    write (1,*) 'x =',x0,' y= ',y0,' q= ',q0
    write (1,*) ''

    y=y0
    x=x0
    z=1-2*x-y
    q=x+y
    write (1,*) ''
    write (1,*) 'Iteration of haplotype frequencies'
    igen=0
 50 if(ind.eq.0) then
    wx=1-(h+q-x*s*(h**2)-0.5*y*h*s)*s
    wy=1-2*(h+(1-h)*q+s*(x+0.5*y)*(h**2)-x*s*h-0.5*y*s)*s
    wz=1-(2*q-y*h*s)*h*s
    else
    wx=1-(h+2*(1-h)*x+(1.5-h)*y-0.5*y*s)*s
    wy=1-2*(h+(1.5-2*h)*x+(1-h)*y+(x+0.5*y)*s*(h**2)-x*s-0.5*y*s)*s
    wz=1-(2*h*q-y*s*h**2)*s
    end if

    wbar=z*wz+2*x*wx+y*wy
    wprod=wz*wy-wx**2

    x1=(x*wx)/wbar
    y1=(y*wy)/wbar
    z1=1-2*x1-y1
    y1=y1+2*x1*u+(z1*u**2)
    x1=x1*(1-u)+z1*u*(1-u)
    igen=igen+1
    q=x1+y1
    write (*,*) 'Gen= ',igen,' x= ',x1,'y= ',y1,' q= ',q
    write (*,*) 'wx= ',wx,' wy= ',wy,' wz= ',wz
    write (*,*) 'wy*wz-wx**2= ',wprod
    write (*,*) 'wbar= ',wbar
    D=q*(1-q)-x1
    R=D/(q*(1-q))
    write(*,*) 'D= ',D,' Correlation= ',R
    x=x1
    y=y1
    z=1-2*x-y
    if(igen.lt.ngen) go to 50
    if(ind.eq.0) then
!   population genetic statistics for sites model
    write (1,*) ''
    write (1,*) 'Population genetic statistics for sites model'
    write (1,*) 'Gen= ',igen,' x= ',x1,'y= ',y1,' q= ',q
    write (1,*) ''
    write (1,*) 'wx= ',wx,' wy= ',wy,' wz= ',wz
    write (1,*) 'wy*wz-wx**2= ',wprod
    write (1,*) 'wbar= ',wbar
    write (1,*) ''
    write (1,*) 'D= ',D,' Correlation= ',R
    Dapp=0
    Rapp=0
    write (1,*) 'Approx. D= ',Dapp,' Approx. Corr.= ',Rapp
    write (1,*) ''
    qr=q/q0
    write (1,*) 'q0= ',q0

    if(h.le.0.001) then
    qr=q/dsqrt(u1/s)
    end if
    write (1,*) 'q relative to approx. equilib. value= ',qr
    
    write (1,*) ''
    go to 10
    end if

    Rapp=(q0**2)*(1-q0)*(0.5-h+(0.5-h)*s)
    Rapp=Rapp/(h+2*(1-3*h)*q0)
    Dapp=Rapp*q0*(1-q0)
    write (1,*) ''
    write (1,*) 'Population genetic statistics for gene model'
    write (1,*) 'Gen= ',igen,' x= ',x1,'y= ',y1,' q= ',q
    write (1,*) 'wx= ',wx,' wy= ',wy,' wz= ',wz
    write (1,*) 'wy*wz-wx**2= ',wprod
    write (1,*) 'wbar= ',wbar
    write (1,*) ''
    D=q*(1-q)-x1
    R=D/(q*(1-q))
    write (1,*) 'D= ',D,' Correlation= ',R
    write (1,*) 'Approx. D= ',Dapp,' Approx. Corr.= ',Rapp
    write (1,*) ''
    write (1,*) 'q0= ',q0
    qr=q/q0
    if(h.le.0.001) then
    qr=q/dsqrt(u1)
    end if
    write (1,*) 'q relative to approx. equilib. value= ',qr
    write (1,*) ''
    write (1,*) 'Load statistics for gene model'
    Lg0=4*z*h*(x+y*(1+e))+4*x*(x+1.5*y*(1+e))+2*(1+e)*y**2
    Lg=Lg0/(4*u1)
!   gene model load relative to 4 x mutation rate (additive value)
    Bg=2*x+y*(1+e)-Lg0
    z2=4*u1*(1.0/(2*h)-1)

    if(h.le.0.001) then
    Bg=Bg/(2*(dsqrt(u1)-u1))
    go to 120
    end if
        if(h.le.0.4999) then
        Bg=Bg/z2
        else
        Bg=1-Bg
        end if
!   inbreeding load relative to additive value
120 write (1,*) ''
    write (1,*) 'Relative L= ',Lg,' Relative B= ',Bg
    write (1,*) ''
    Vg=4*z*x*h**2+4*z*y*(h*(1+e))**2+4*x**2+4*x*y*(1.5*(1+e))**2
    Vg=Vg+(y*2*(1+e))**2-Lg0**2
    Vg=Vg/(4*u1)
    write (1,*) 'Relative V x h = ',Vg
    write (1,*) ''
    go to 10

100 end program dom5
    
