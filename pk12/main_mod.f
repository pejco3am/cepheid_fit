
	program main

        implicit real*8 (a-h,o-z)
        include 'bindat01.com'

        character*30 ifile1,ifile2,ifile3,ifile4
        character*30 ofile1,ofile2,ofile3,ofile4,ofile5,ofile6

        character*30 text(20)
        character*9 tname
        character*2 cgal
        character*5 cnum

c temporary storage for magitude and error output lines
        real mout(NCMAX),eout(NCMAX)
c plotting variables
        real*8 tphi(360),temp(360),rad(360),vrad(360)
c covariances variables
        real*8 egnv1, egnv2, discriminant, nx1, ny1, nx2, ny2, denomdir
c vector to be minimized and its derivative
        real*8 p(NVMAX),dp0(NVMAX),dp1(NVMAX)
c number of measurements available in color i
        integer nmcol(NCMAX)
c temporary storage of magnitude in color i
        real*8 tv(NCMAX)
c for reference i, refvnum(i)=number of velocity measurements, refvres(i)=rms residual
c for reference i, refpnum(i,j)=number of magnitude measurements in color j, refvres(i,j)=rms residual
        integer refvnum(NREFMAX)
        real*8  refvres(NREFMAX)
        real*8  refvbar(NREFMAX)
        integer refpnum(NREFMAX,NCMAX)
        real*8  refpres(NREFMAX,NCMAX)
        real*8  refpbar(NREFMAX,NCMAX)
        real   tmin,tmax
c total number of measurements for each cepheid
        integer nmeastot(NCMAX)
c period derivative sensitivity
        real*8 sense(NFMAX)

        character*2 squote
        character*1 quote

c temporary variables to fix templates
        real*8 phase(1000),amp(1000)

c this puts in a single quote
        quote  = "'"
        squote = " '"

        ln10 = log(10.0D0)
        pi   = 4.0D0*atan(1.0D0)
c velocity expansion p-factor
        pexp = 1.36
c units are 10 solar radii, km/s and days, 100 cm/s^2, solar masses
c 80.55 = 10 Solar radii per day
        vcon  = -80.55*ln10/pexp
        acon  =  93.23/100.0
        mcon  =  0.340


c clear the velocity/magnitude zero point values, set default
c uncertainties for them -- clear various statistical vectors
c for number of data points and typical ucertainties
        do i=1,NREFMAX
          refvnum(i)       = 0
          refvres(i)       = 0.0
          refvbar(i)       = 0.0
          evrms(i)         = 0.3
          evbar(i)         = 0.0
          izero(i,NCMAX+1) = 0
          mzero(i,NCMAX+1) = 0.0D0
          mzeroe(i,NCMAX+1)= 1.0D0
          do j=1,NCMAX
            refpnum(i,j)   = 0
            refpres(i,j)   = 0.0
            refpbar(i,j)   = 0.0
            emrms(i,j)     = 0.02
            embar(i,j)     = 0.00
            izero(i,j)     = 0
            mzero(i,j)     = 0.0
            mzeroe(i,j)    = 0.01
            enddo
          enddo



c READ IN THE DATA FOR THE GALAXIES THAT WILL BE FIT
        print*,'fit how many galaxies '
        read*,ngal
        nc = 1
        nd = 0

        do kk=1,ngal
          print*,'enter galaxy ID# '
          read*,idgalsave
          if (idgalsave.gt.NGMAX) then
            print*,'exceeded allowed number of galaxies '
            stop
            endif
c igal0 provides a pointer to the first Cepheid in the galaxy = k-igal0(idgal)
c where k is the Cepheid number in the order read in
          igal0(idgalsave) = nc-1
c ignum/ignum0 allow you to match the true galaxy ID with its array number and 
c vice versa
          ignum(kk)        = idgalsave
          ignum0(idgalsave)= kk
c ncl will count the number of Cepheids read for the current galaxy
          ncl              = 1
          print*,'enter input file 1 for galaxy ',idgalsave
          read(*,'(a30)')ifile1
  	  open(unit=13,file=ifile1,form='unformatted',status='old')
c now read in the Cepheids -- in the Galaxy and clouds there is
c position information (xx/yy) otherwise set to zero
100       if (idgalsave.le.3) then 
            read(13,end=110)npt(nc),name(nc),period(nc),xx(nc),yy(nc)
          else
            read(13,end=110)npt(nc),name(nc),period(nc)
            xx(nc) = 0.0
            yy(nc) = 0.0
            endif
c this is the literature period of the Cepheid -- save its log for use in 
c scaling relations
            lper0(nc)  = log10(period(nc)/10.0)
c the data for Cepehid nc will run from i0(nc) to i1(nc) 
            i0(nc)     = nd + 1
            i1(nc)     = nd + npt(nc)
c store the galaxy ID of the cepheid, and its number in the current galaxy
            idsys0(nc) = idgalsave
            idnum0(nc) = ncl
c now read the data for the Cepheid
            do i=1,npt(nc)
c increment the total number of data points read, associate data point with the
c current Cepheid number
             nd     = nd + 1
             ic(nd) = nc
             if (nd.gt.NPMAX) then
               print*,'TOO MUCH PHOTOMATRIC DATA REDIMENSION ARRAYS ',nd,NPMAX
               stop
               endif
             read(13)time(nd),vr(nd),evr(nd),ndat(nd),iref(nd)
c if there is a radial velocity measurement, turn on the associated velocity
c zero point for that reference
             if (vr(nd).gt.-500.0) then
               izero(iref(nd),NCMAX+1) = 1
               endif
             if (ndat(nd).gt.5) then
               print*,'TOO MANY COLORS IN ONE LINE ',ndat(nd),iref(nd)
               stop
               endif
c read the photometric data, turn on magnitude zero points for the reference
             if (ndat(nd).gt.0) then
               read(13)(icdat(nd,jj),jj=1,ndat(nd))
               read(13)(mag  (nd,jj),jj=1,ndat(nd))
               read(13)(emag (nd,jj),jj=1,ndat(nd))
               do jj=1,ndat(nd)
                 izero(iref(nd),icdat(nd,jj)) = 1
                 enddo
               endif
             enddo
c increment global Cepheid count
           nc = nc + 1
           if (nc.gt.NFMAX) then
             print*,'Cepheid number limit ',NFMAX,' exceeded '
             print*,'REDIMENSION ARRAYS '
             stop
             endif
c increment local Cepheid count
           ncl= ncl + 1
           go to 100
110      close(unit=13)
         enddo

        nc = nc - 1


c write out a summary of these mappings for doing checks
c   do i=1,nc
c      write(85,'(i4,i4,i4,1x,a10,1x,i5,1x,f8.3)')i,idsys0(i),idnum0(i),name(i),i1(i)-i0(i),period(i)
c   enddo

c set the number of colors to be used

      open(unit=400,file='ven.dat',form='formatted',status='unknown') 

        ncol   = 36
        do k=1,nc
          ii0 = i0(k)
          ii1 = i1(k)
            do i=ii0,ii1
               if (vr(i).gt.-500.0) then
               write(400,'(i2,1x,i3,1x,f10.3,2(1x,f7.3),1x,i3,1x,a9)')0,iref(i),
     1           time(i),vr(i),evr(i),k,name(k)
               endif
              do j=1,ndat(i)
                icol = icdat(i,j)
                if (mag(i,j).gt.-5.0) then
                  write(400,'(i2,1x,i3,1x,f10.3,2(1x,f7.3),1x,i3,1x,a9)')icol,iref(i),
     1               time(i),mag(i,j),emag(i,j),k,name(k)
                  endif
               enddo
            enddo
        enddo
      close(unit=400)
      end
