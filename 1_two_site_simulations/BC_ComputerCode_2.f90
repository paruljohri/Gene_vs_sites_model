!The programs are written in GFortran for the MacIntosh 
!    (https://gnuc.org.wiki/GFortranBinariesMacOS)
!
!   Written by Brian Charlesworth. Contact: brian.charlesworth@ed.ac.uk for questions.
!
!Program for Simulating Two Segregating Loci
!   This program simulates a haploid population with two segregating sites with no   
!   recombination, allowing weights to be applied towards low allele frequencies by the 
!   methods of Garcia & Lohmueller (2021) or Good (2022)

    program dsim11
!   program for 2-locus D (no recombination)
!   finite sample size
    CHARACTER*20 FOUT
    CHARACTER*20 FIN
    integer :: N,N1,nrep,ngen,ns,isamp,tsamp1,tsamp2,tsamp
    real :: AN,f(0:1000),g(1000),r,q1b1s,q1b2s,q2b1s,q2b2s,D1bs,D2bs,pr1bs,pr2bs
    real :: y1,y3,y4

    write(*,*) 'Input file?'
    read(*,*)  FIN
    write (*,*) 'Output file?'
    read  (*,*)  FOUT
    OPEN (1,FILE=FOUT)
    OPEN (2,FILE=FIN)

    read(2,*) s
    read(2,*) eps
    read(2,*) N
    read(2,*) nrep
    write(*,*) 'Allele frequency weight?'
    read(*,*) f0
!   this is the reciprocal of the Good (2022) weight
!   f0=0 if no weight is applied
    write(*,*) 'Sample size?'
    read(*,*) ns
    write(*,*) 'Is sampling required (enter 0 for next item if not)?'
    write(*,*) 'Sample frequency selection?'
    read(*,*) isamp

    AN=N
    anrep=nrep
    ANs=2*AN*s
!   scaled selection coefficent

    write(1,*) 'Program for 2-locus D with no recombination'
    write(1,*) 'Finite sample size'
    write(1,*) 'Size of haploid population= ',N
    write(1,*) 'Sample size= ',ns
    write(1,*) 'Sample frequency selection= ',isamp
    write(1,*) 'Selection coefficient= ',s
    write(1,*) 'Scaled selection coefficient= ',ANs
    write(1,*) 'Epistasis coefficient= ',eps
    write(1,*) 'Number of replications= ',nrep
    write(1,*) 'Allele frequency weight= ',f0
    write(1,*) ''

    write(*,*) 'Program for 2-locus neutral D with no recombination'
    write(*,*) 'Finite sample size'
    write(*,*) 'Size of haploid population= ',N
    write(*,*) 'Sample size= ',ns
    write(*,*) 'Sample frequency selection= ',isamp
    write(*,*) 'Selection coefficient= ',s
    write(*,*) 'Scaled selection coefficient= ',ANs
    write(*,*) 'Epistasis coefficient= ',eps
    write(*,*) 'Number of replications= ',nrep
    write(*,*) 'Allele frequency weight= ',f0
    write(*,*) ''

    write(*,*) 'Continue?'
    read (*,*) CONT
    if(CONT.eq.1) go to 400

    N1=N-1
    f(0)=0
    
    do 10 i=1,N1
    ai=i
    if(ANs.le.0.0001) then
    f(i)=f(i-1)+1/ai
!   cumulative frequency distribution with no selection
    else
    qi=ai/AN
    pi=1-qi
    g(i)=(exp(ANs*pi)-1)/(exp(ANs)-1)
    g(i)=g(i)/(pi*qi)
    f(i)=f(i-1)+g(i)
!   cumulative frequency distribution with selection
    end if
10  continue

    do 20 i=1,N1
    f(i)=f(i)/f(N1)
!   write(*,*) 'i= ',i,' f(i)= ',f(i)
20  continue
!   normalized initial cumulative frequency distribution

    q1bar1=0
    q1bar2=0
    q2bar1=0
    q2bar2=0
    D1bar=0
    D2bar=0
    D1bar2=0
    D2bar2=0
    pr1bar=0
    pr2bar=0
    q1b1s=0
    q1b2s=0
    q2b1s=0
    q2b2s=0
    D1bs=0
    D2bs=0
    D1b2s=0
    D2b2s=0
    pr1bs=0
    pr2bs=0
    
!   initializes sums and sums of squares of frequency and D statistics for the 2 initial states
!   pr is the crossproduct of allele frequencies
    ns1=0
    ns2=0
    tns1=0
    tns2=0
    tsamp1=0
    tsamp2=0
!   initializes numbers of 2 initial states and sums of their sojourn times

    w2=1-s
    w3=1-s
    w4=1-2*s*(1+eps)
!   haplotype fitnesses


    call init_random_seed()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   starts simulations
    do 100 irep=1,nrep
!   write(*,*) ''
!   write(*,*) 'rep= ',irep
    call random_number(r)
    z0=r
!   write(*,*) 'rn= ',z0
    do 30 i=1,N1
    if(z0.le.f(i)) then
    ai=i
!   write(*,*) 'i= ',i
    q0=ai/AN
    go to 31
    end if
    
30  continue
!   sets initial allele frequency for 1st mutation
31  p0=1-q0

!   write(*,*) 'q0=',q0
    call random_number(r)
    z1=r
!   write(*,*) 'rn= ',z1
    if(z1.le.q0) then
    i1=1
    ns1=ns1+1
    else
    i1=2
    ns2=ns2+1
    end if
!   decides if 2nd mutation is in coupling (1) or repulsion (2) with 1st one

    if(i1.eq.1) then
    x10=p0
    x20=0
    x40=1.0/AN
    x30=q0-x40
    else
    x20=1.0/AN
    x10=p0-x20
    x30=q0
    x40=0
    end if
!   sets up initial haplotype frequencies for each initial state

!   write(*,*) 'Initial state=',i1
    x1=x10
    x2=x20
    x3=x30
    x4=x40
    j1=0
!   write(*,*) 'x1,x2,x3,x4 =',x1,x2,x3,x4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   starts next generation
35  j1=j1+1
    
!   write(*,*) 'Generation',j1
    wbar=x1+(x2+x3)*w2+x4*w4
    x1=x1/wbar
    x2=x2*(1-s)/wbar
    x3=x3*(1-s)/wbar
    x4=1-x1-x2-x3
!   post-selection haplotype frequencies

    y1=x1
    if(i1.eq.1) then
    y2=x3
    y3=x4
    else
    y2=x2
    y3=x3
    end if
!   relabels haplotypes to match the 2 different initial states
    y0=y1+y2
    y5=y1/y0
    n1=0
    n2=0
    n3=0
!   write(*,*) 'y1,y2,y3',y1,y2,y3
    do  40 j2=1,N
    call random_number(r)
    z1=r
!   write(*,*) 'random no 1 = ',z1
    if(z1.gt.y0) then
    n3=n3+1
    else
    call random_number(r)
    z2=r
!   write(*,*) 'random no 2 = ',z2
     if(z2.le.y5) then
     n1=n1+1
     end if
    end if
40  continue
    n2=N-n1-n3
!   counts numbers of haplotypes after multinomial sampling

!   write(*,*) 'n1,n2,n3',n1,n2,n3
    if(n3.eq.0) then
      if(i1.eq.1) then
      q2n=0
      tn1=j1
      else
      q1n=0
      tn2=j1
      end if
    go to 50
    end if
!   loss of 1st or 2nd mutation

    if(n3.eq.N) then
     if(i1.eq.1) then
     q2n=1
     tn1=j1
     else
     q1n=1
     tn2=j1
     end if
    go to 50
    end if
!   fixation of 1st or 2nd mutation

    if(i1.eq.2) then
      if(n2.eq.0) then
      q2n=0
      tn2=j1
      go to 50
      end if
!   loss of 2nd mutation
      if(n2.eq.N) then
      q1n=1
      tn2=j1
      go to 50
!   fixation of 2nd mutation
      end if
    else
    n4=n2+n3
      if(n4.eq.0) then
      q1n=0
      tn1=j1
!   loss of 1st mutation
      go to 50
      end if
      if(n4.eq.N) then
      q1n=1
      tn1=j1
!   fixation of 1st mutation
      go to 50
      end if
    end if
!   fixations or losses of mutations; run terminated; times tn1 or tn2 recorded

    y11=n1/AN
    y31=n3/AN
    y21=1-y11-y31
!   new frequencies from multinomial sampling of haplotypes

    if(isamp.gt.0) then
    call sample(ns,isamp,y11,y21,y31,i1,q1b1s,q2b1s,q1b2s,q2b2s,D1bs,D2bs,D1b2s,D2b2s,pr1bs,pr2bs,tsamp1,tsamp2)
!   estimates sums over generations and replicates of sample statistics for chosen sample allele frequencies
    end if

    x1=y11
    if(i1.eq.1) then
    x2=0
    x3=y21
    x4=y31
    q1bar1=q1bar1+x3+x4
    q2bar1=q2bar1+x4
    D1=x1*x4*exp(0.0-f0*(x3+2*x4))
    D1bar=D1bar+D1
    D1bar2=D1bar2+D1**2
    pr1bar=pr1bar+(x3+x4)*(1-x3-x4)*x4*(1-x4)*exp(0.0-f0*(x3+2*x4))
    else
    x2=y21
    x3=y31
    x4=0
    q1bar2=q1bar2+x3
    q2bar2=q2bar1+x2
    D2=0-x2*x3*exp(0.0-f0*(x2+x3))
    D2bar=D2bar+D2
    D2bar2=D2bar2+D2**2
    pr2bar=pr2bar+x2*(1-x2)*x3*(1-x3)*exp(0.0-f0*(x2+x3))
    end if
!   D statistics with Good's (2022) weights
!   write(*,*) 'x1,x2,x3,x4',x1,x2,x3,x4
    
!   reverses haplotype labels and accumulates sums of q's and D statistics over all generations
    go to 35

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   replicate ended due to fixation or loss at either locus
50  if(i1.eq.1) then
    tns1=tns1+tn1
    
    else
    tns2=tns2+tn2
    end if
!   accumulates sums of sojourn times
    
100 continue
    write(*,*) 'Simulation completed'

    tns=tns1+tns2
!   sum of the two sums of sojourn times
    Dbar=(D1bar+D2bar)/tns
    Dbar2=(D1bar2+D2bar2)/tns
    D1bar=D1bar/tns1
    D1bar2=D1bar2/tns1
    D2bar=D2bar/tns2
    D2bar2=D2bar2/tns2
    prbar=(pr1bar+pr2bar)/tns
    pr1bar=pr1bar/tns
    pr2bar=pr2bar/tns
!   sums of D statistics normalised by sums of sojourn times

    write(1,*) 'Summary statistics for whole population'
    write(1,*) ''
    ans1=ns1
    ans2=ns2
    ps1=ans1/anrep
    ps2=ans2/anrep
    q1bar=(q1bar1+q1bar2)/tns
    q1bar1=q1bar1/tns1
    q1bar2=q1bar2/tns2
    q2bar=(q2bar2+q2bar2)/tns
    q2bar1=q2bar1/tns1
    q2bar2=q2bar2/tns2
    write(1,*) 'ns1= ',ns1,' ps1= ',ps1
    write(1,*) 'ns2= ',ns2,' ps2= ',ps2
    write(1,*) ''
    write(1,*) 'Mean allele frequencies at 1st locus, conditioned on segregation'
    write(1,*) 'q1bar1= ',q1bar1,' q1bar2= ',q1bar2,' q1bar= ',q1bar
    write(1,*) 'Mean allele frequencies at 2nd locus, conditioned on segregation'
    write(1,*) 'q2bar1= ',q2bar1,' q2bar2= ',q2bar2,' q2bar= ',q2bar
    write(1,*) ''
    tbar1=tns1/ans1
    tbar2=tns2/ans2
    tbar=(tns1+tns2)/anrep
    write(1,*) 'Mean sojourn times'
    write(1,*) 'tbar1= ',tbar1,' tbar2= ',tbar2,' tbar= ',tbar
    write(1,*) ''
    write(1,*) 'D statistics'
    write(1,*) 'Mean D1= ',D1bar,'Mean D2= ',D2bar,' Mean D= ',Dbar
    write(1,*) 'Mean D1^2= ',D1bar2,'Mean D2^2= ',D2bar2,' Mean D^2= ',Dbar2
    write(1,*) ''
    write(1,*) 'Mean pr1= ',pr1bar,' Mean pr2= ',pr2bar,' Mean pr= ',prbar
    sig1d1=D1bar/pr1bar
    sig1d2=D2bar/pr2bar
    sig1d=Dbar/prbar
    write(1,*) 'sig1d1= ',sig1d1,' sig1d2= ',sig1d2,' sig1d= ',sig1d
    sqtpr1=sqrt(pr1bar)
    sqtpr2=sqrt(pr2bar)
    sqtpr=sqrt(prbar)
    sigd1=D1bar/sqtpr1
    sigd2=D2bar/sqtpr2
    sigd=Dbar/sqtpr
    write(1,*) 'sigd1= ',sigd1,' sigd2= ',sigd2,' sigd= ',sigd
    sig2d1=D1bar2/pr1bar
    sig2d2=D2bar2/pr2bar
    sig2d=Dbar2/prbar
    write(1,*) 'sig2d1= ',sig2d1,' sig2d2= ',sig2d2,' sig2d= ',sig2d
    sig2d=0.454545
    write(1,*) 'Theoretical value of unweighted sigma2d= ',sig2d
    write(1,*) ''
    write(1,*) ''

    if(isamp.eq.0) go to 400
    write(1,*) 'Summary statistics for sample'
    write(1,*) ''
    atsamp1=tsamp1
    atsamp2=tsamp2
    tsamp=atsamp1+atsamp2
    atsamp=tsamp
    write(1,*) 'Sums of generation times for samples which match chosen allele frequency'
    write(1,*) 'tsamp1= ',tsamp1,' tsamp2= ',tsamp2,' tsamp= ',tsamp

    q1bs=(q1b1s+q1b2s)/atsamp
    q1b1s=q1b1s/atsamp1
    q1b2s=q1b2s/atsamp2
    q2bs=(q2b1s+q2b2s)/atsamp
    q2b1s=q2b1s/atsamp1
    q2b2s=q2b2s/atsamp2
    write(1,*) 'Mean allele frequencies at 1st locus, conditioned on segregation'
    write(1,*) 'q1bar1= ',q1b1s,' q1bar2= ',q1b2s,' q1bar= ',q1bs
    write(1,*) 'Mean allele frequencies at 2nd locus, conditioned on segregation'
    write(1,*) 'q2bar1= ',q2b1s,' q2bar2= ',q2b2s,' q2bar= ',q2bs
    write(1,*) ''

    Dbs=(D1bs+D2bs)/atsamp
    D1bs=D1bs/atsamp1
    D2bs=D2bs/atsamp2
    Db2s=(D1b2s+D2b2s)/atsamp
    D1b2s=D1b2s/atsamp1
    D2b2s=D2b2s/atsamp1
    prbs=(pr1bs+pr2bs)/atsamp
    pr1bs=pr1bs/atsamp1
    pr2bs=pr2bs/atsamp2
!   sums of D statistics normalised by sums of generation times
    write(1,*) 'D statistics'
    write(1,*) 'Mean D1= ',D1bs,'Mean D2= ',D2bs,' Mean D= ',Dbs
    write(1,*) 'Mean D1^2= ',D1b2s,'Mean D2^2= ',D2b2s,' Mean D^2= ',Db2s
    write(1,*) ''
    write(1,*) 'Mean pr1= ',pr1bs,' Mean pr2= ',pr2bs,' Mean pr= ',prbs

    sig1d1=D1bs/pr1bs
    sig1d2=D2bs/pr2bs
    sig1d=Dbs/prbs
    write(1,*) 'sig1d1= ',sig1d1,' sig1d2= ',sig1d2,' sig1d= ',sig1d
    sqtpr1=sqrt(pr1bs)
    sqtpr2=sqrt(pr2bs)
    sqtpr=sqrt(prbs)
    sigd1=D1bs/sqtpr1
    sigd2=D2bs/sqtpr2
    sigd=Dbs/sqtpr
    write(1,*) 'sigd1= ',sigd1,' sigd2= ',sigd2,' sigd= ',sigd
    sig2d1=D1b2s/pr1bs
    sig2d2=D2b2s/pr2bs
    sig2d=Db2s/prbs
    write(1,*) 'sig2d1= ',sig2d1,' sig2d2= ',sig2d2,' sig2d= ',sig2d
    
400 end program dsim11


    subroutine sample(ns,isamp,y1,y2,y3,i1,q1b1s,q2b1s,q1b2s,q2b2s,D1bs,D2bs,D1b2s,D2b2s,pr1bs,pr2bs,tsamp1,tsamp2)
!   sampling of haplotypes from population
    integer :: ns,isamp,i1,tsamp1,tsamp2,nsamp1,nsamp2
    real :: y1,y2,y3,q1b1s,q2b1s,q1b2s,q2b2s,D1bs,D2bs,D1b2s,D2b2s,pr1bs,pr2bs

!   write(*,*) 'y1,y2,y3 ',y1,y2,y3
!   write(*,*) 'i1=',i1
    y0=y1+y2
    y5=y1/y0
    n1=0
    n2=0
    n3=0
    ans=ns
    do  10 j=1,ns
    call random_number(r)
    z1=r
!   write(*,*) 'random no 1 = ',z1
    if(z1.gt.y0) then
    n3=n3+1
    else
    call random_number(r)
    z2=r
!   write(*,*) 'random no 2 = ',z2
      if(z2.le.y5) then
      n1=n1+1
      end if
    end if
10  continue
    n2=ns-n1-n3
!   counts numbers of haplotypes after multinomial sampling into sample of size ns

    n11=n1
    x1=n11/ans
    if(i1.eq.1) then
    n12=0
    n13=n2
    n14=n3
    na1=n13+n14
    na2=n14
    x2=0
    x3=n13/ans
    x4=1-x1-x3
!   write(*,*) 'na1=',na1,' na2=',na2
        if(na1.eq.isamp) then
        if(na2.eq.isamp) go to 20
        end if
        go to 100
!   exits subroutine and does not collect statistics if neither locus has right state
    
20      q1b1s=q1b1s+x3+x4
        q2b1s=q2b1s+x4
        D1s=x1*x4
        D1bs=D1bs+D1s
        D1b2s=D1b2s+D1s**2
        pr1bs=pr1bs+(x3+x4)*(1-x3-x4)*x4*(1-x4)
        tsamp1=tsamp1+1
!       write(*,*) 'tsamp1=',tsamp1
        go to 100
    else
    n12=n2
    n13=n3
    n14=0
    na1=n13
    na2=n12
    x4=0
    x2=n12/ans
    x3=1-x1-x2
!   write(*,*) 'na1=',na1,' na2=',na2
    if(na1.eq.isamp) then
    if(na2.eq.isamp) go to 40
    end if
    go to 100
!   exits subroutine and does not collect statistics

 40     q1b2s=q1b2s+x3
        q2b2s=q2b2s+x2
        D2s=0-x2*x3
        D2bs=D2bs+D2s
        D2b2s=D2b2s+D2s**2
        pr2bs=pr2bs+x2*(1-x2)*x3*(1-x3)
        tsamp2=tsamp2+1
!   write(*,*) 'tsamp2=',tsamp2
    end if
!   reverses haplotype labels and accumulates sums of q,t and D statistics


100 end subroutine sample


    subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
     form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
   read(un) seed
   close(un)
    else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
   call system_clock(t)
   if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
    end if
    call random_seed(put=seed)
    contains
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
    function lcg(s)
  integer :: lcg
  integer(int64) :: s
  if (s == 0) then
     s = 104729
  else
     s = mod(s, 4294967296_int64)
  end if
  s = mod(s * 279470273_int64, 4294967291_int64)
  lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
    end subroutine init_random_seed
    
