
!The programs are written in GFortran for the MacIntosh 
!(https://gnuc.org.wiki/GFortranBinariesMacOS)



!Program for Simulating Samples from a Population with Two Segregating Loci
!This program simulates repeatied sampling from a population with two loci with the same allele frequency at each locus 
!and a specified value of  D. The full set of D values for a fixed number of replicates can be stored.
	
    program Dsampa
!   program for 2-locus D with finite sample size

    CHARACTER*20 FOUT
    CHARACTER*20 FIN
    integer :: ns,nrep,nind,nseg,n1,n2,n3,n4
    real :: A
    real :: y1,y2,y3,y4,x1,x2,x3,x4,DH(10000)

    write(*,*) 'Input file?'
    read(*,*)  FIN
    write (*,*) 'Output file?'
    read  (*,*)  FOUT
    OPEN (1,FILE=FOUT)
    OPEN (2,FILE=FIN)

    read(2,*) q
    read(2,*) D
    read(2,*) ns
    read(2,*) nrep
    read(2,*) nhist

    AN=N
    anrep=nrep
    y2=q*(1-q)-D
    y3=y2
    y4=D+q**2
    y1=1-y2-y3-y4
!   Population haplotype frequencie

    write(1,*) 'Program for 2-locus D in finite sample'
    write(1,*) 'Allele frequency= ',q
    write(1,*) 'Population D= ',D
    write(1,*) 'y1= ',y1,' y2= ',y2,' y3= ',y3,' y4= ',y4
    write(1,*) ''
    write(1,*) 'Sample size= ',ns
    write(1,*) 'Number of replicates= ',nrep
    write(1,*) 'Number of values for histogram= ',nhist
    write(1,*) ''

    write(*,*) 'Program for 2-locus D in finite sample'
    write(*,*) 'Allele frequency= ',q
    write(*,*) 'Population D= ',D
    
    write(*,*) 'Sample size= ',ns
    write(*,*) 'Number of replicates= ',nrep
    write(*,*) ''
   
    write(*,*) 'Continue?'
    read (*,*) CONT
    if(CONT.eq.1) go to 400

    
    q1s=0
    q1ss=0
    q2s=0
    q2ss=0
    Ds=0
    Dnegs=0
    Dpos=0
    Dss=0
    nseg=0
    nx4=0
    ndneg=0
    ndpos=0
    x2s=0
    x3s=0
    x4s=0

    do 10 i=1,nrep
    call sample(ns,y1,y2,y3,y4,x1,x2,x3,x4,n1,n2,n3,n4,nind)
    
    if(nind.eq.0) go to 10
    nseg=nseg+nind
!   only increments if nind = 1
!   this happens if both loci segregate
    x2s=x2s+x2
    x3s=x3s+x3
    x4s=x4s+x4

    q1=x3+x4
    q2=x2+x4

    if(n4.eq.0) then
    D=0.0-x2*x3
    else
    nx4=nx4+1
    D=x1*x4-x2*x3
    end if

    if(nseg.le.nhist) then
    DH(nseg)=D
    end if

    iD=n1*n4-n2*n3
    if(iD.lt.0) then
    ndneg=ndneg+1
    Dnegs=Dnegs+D
    else
    ndpos=ndpos+1
    Dpos=Dpos+D
    end if

    q1s=q1s+q1
    q1ss=q1ss+q1**2
    q2s=q2s+q2
    q2ss=q2ss+q1**2
    Ds=Ds+D
    Dss=Dss+D**2
10 continue

    anseg=nseg
    andneg=ndneg
    andpos=ndpos
    x2m=x2s/anseg
    x3m=x3s/anseg
    x4m=x4s/anseg
    px4=nx4/anseg
    pdneg=ndneg/anseg
    q1m=q1s/anseg
    q1v=(q1ss-(q1s**2)/anseg)/(anseg-1)
    q1se=sqrt(q1v/anseg)
    q2m=q2s/anseg
    q2v=(q2ss-(q2s**2)/anseg)/(anseg-1)
    q2se=sqrt(q2v/anseg)

    Dm=Ds/anseg
    D2m=Dss/anseg
    Dv=(Dss-(Ds**2)/anseg)/(anseg-1)
    Dse=sqrt(DV/anseg)
    rm=Dm/(sqrt(q1m*q2m))
    r2m=D2m/(q1m*q2m)

    Dnegm=Dnegs/andneg
    Dposm=Dpos/andpos
    
    write(1,*) 'Sample statistics'
    write(1,*) 'Mean x2,x3,x4= ',x2m,x2m,x4m
    write(1,*) 'Number of samples with 2 segregating sites= ',nseg
    write(1,*) 'No. of samples with A2B2 present= ',nx4,' Prop.= ',px4
    write(1,*) 'No. of samples with D<0 ',ndneg,' Prop.= ',pdneg
    write(1,*) 'Mean negative D= ',Dnegm,' Mean positive D= ',Dposm
    write(1,*) 'Mean q1= ',q1m,' Variance= ',q1v,' s.e.= ',q1se
    write(1,*) 'Mean q2= ',q2m,' Variance= ',q2v,' s.e.= ',q2se
    write(1,*) 'Mean D= ',Dm,' Variance= ',Dv,' s.e. = ',Dse
    write(1,*) 'Mean D^2= ',D2m
    write(1,*) 'Mean sigma-d= ',rm,' Mean sigma-d^2= ',r2m
    write(1,*) ''
    write(1,*) 'Distribution of D'
    do 20 j=1,nhist
    write(1,*) DH(j)
20  continue

400 end program Dsampa


    subroutine sample(ns,y1,y2,y3,y4,x1,x2,x3,x4,n1,n2,n3,n4,nind)
!   sampling of haplotypes from population
    integer :: ns,n1,n2,n3,n4,nind
    real :: y1,y2,y3,y4,x1,x2,x3,x4
!   y's are the population frequencies of the haplotypes
!   write(*,*) 'y1,y2,y3,y4 ',y1,y2,y3,y4

    n1=0
    n2=0
    n3=0
    n4=0
    ans=ns
    do  10 j=1,ns
    call random_number(r)
    z1=r
!   write(*,*) 'random no 1 = ',z1
    if(z1.le.y1) then
!   A1B1 is chosen
    n1=n1+1
    go to 10
    end if

    x1=y2/(y2+y3+y4)
    call random_number(r)
    z2=r
!   write(*,*) 'random no 2 = ',z2
      if(z2.le.x1) then
      n2=n2+1
!   A1B2 is chosen
      go to 10
      end if

        call random_number(r)
        z3=r
        x2=y3/(y3+y4)
          if(z3.le.x2) then
          n3=n3+1
!   A2B1 is chosen
          else
          n4=n4+1
!   A2B2  is chosen
          end if
10  continue
    
    n11=n1+n2
    n12=n1+n3
    if(n11.eq.0) go to 20
    if(n11.eq.ns) go to 20
    if(n12.eq.0) go to 20
    if(n12.eq.ns) go to 20
    go to 30

20  nind=0
    go to 100
!   discards sample that does not segregate at both sites
30  nind=1
    x2=n2/ans
    x3=n3/ans
    x4=n4/ans
    x1=1-x2-x3-x4
    
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




    
