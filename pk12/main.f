
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
        do i=1,nc
          write(85,'(i4,i4,i4,1x,a10,1x,i5,1x,f8.3)')i,idsys0(i),idnum0(i),name(i),
     1       i1(i)-i0(i),period(i)
          enddo

c set the number of colors to be used
        ncol   = 36

c read in the templates: ncos = Fourier order, ntempl = number of period terms
c tlp = defining periods for the templates
        open(unit=13,file='template.dat',form='formatted',status='old')
        read(13,*)ncos(1),ntempl(1)
        read(13,*)(tlp(1,i),i=1,ntempl(1))
c read in templates for fundamental mode
        ntem(1) = 2*ncos(1)
        do j=1,ncos(1)
          read(13,*)jj,(ctt(1,i,j),stt(1,i,j),i=1,ntempl(1))
          enddo
        do j=1,ncos(1)
          read(13,*)jj,(crt(1,i,j),srt(1,i,j),i=1,ntempl(1))
          enddo
c read in templates for overtone mode
        print*,'getting rid of overtones '
        ncos(2)   = 0
        ntempl(2) = 0
        close(unit=13)


c this is just a grid of phases for plotting light curves
        print*,'setting phase grid '
        nplt = 360
        dphi = 360.0/float(nplt)
        do i=1,nplt
          tphi(i)      = dphi*float(i-1)
          enddo

      print*,' setting which templated in period '
      call settemplate(1,1.35)
      print*,' building old template '
      do i=1,nplt
        ang     = pi*tphi(i)/180.0
        call gettemplate(0,1,1,-1,ang)
        rtin = rt(1)
        ttin = tt(1)
        eps  = 0.001
        call gettemplate(0,1,1,-1,ang+eps)
        gradr = (rt(1)-rtin)/eps
        gradt = (tt(1)-ttin)/eps
        write(13,'(10(1x,f13.6))')ang,rt(1),tt(1),rt(2),tt(2),gradr,gradt
        enddo


      iconvert = 0
      if (iconvert.eq.1) then
        print*,' renormalizing templates '
c convert the phases
c CONVERSION CODE TO RENOMALIZE TO NEW TEMPLATE FORM
c USE AMP/PHASE TO RENORMALIZE AMPLITUDES and PHASES
c in order to make the conversion, one seems to require
c    b0 ->  -b0/2.5
c if rescale radial template amplitude by F
c    r  ->   F r               = r/amp
c    t  ->   F t ln10          = ln 10 (t/amp)
c    A  ->  sqrt(A/F ln 10)    = sqrt (A amp / ln10)
c where the new A is also now squared to get the actual amplitude
        do i=1,ntempl(1)
          amp(i)    = sqrt(crt(1,i,1)**2+srt(1,i,1)**2) 
          phase(i)  = atan2(srt(1,i,1),crt(1,i,1))
          if (phase(i).lt.0.0) phase(i) = phase(i) + 2.0*pi
          print*,' for term ',i,' amp/phase ',amp(i),phase(i)
          do j=1,ncos(1)
            rj    = float(j)
            tempc =  crt(1,i,j)*cos(rj*phase(i))+srt(1,i,j)*sin(rj*phase(i))
            temps = -crt(1,i,j)*sin(rj*phase(i))+srt(1,i,j)*cos(rj*phase(i))
            crt(1,i,j) = tempc/amp(i)
            srt(1,i,j) = temps/amp(i)
            tempc =  ctt(1,i,j)*cos(rj*phase(i))+stt(1,i,j)*sin(rj*phase(i))
            temps = -ctt(1,i,j)*sin(rj*phase(i))+stt(1,i,j)*cos(rj*phase(i))
c there is no sign flip here 
            ctt(1,i,j) = (tempc/amp(i))*ln10
            stt(1,i,j) = (temps/amp(i))*ln10
            if (i.eq.1) then
              write(*,'(a4,1x,10(1x,f13.6))')'new ',crt(1,i,j),srt(1,i,j),ctt(1,i,j),stt(1,i,j)
              endif
            enddo
          enddo
      print*,'setting template '
      juse = 3
      call settemplate(1,tlp(1,juse))
      print*,' building new template '
      do i=1,nplt
        ang     = pi*tphi(i)/180.0
        call gettemplate(0,1,1,-1,ang)
        write(14,'(3(1x,f13.6))')ang+phase(juse),rt(1)*amp(juse),tt(1)*amp(juse)
        enddo
      endif


c read in the zero point shifts (mzero) and their uncertainties (mzeroe)
c can optimize mzero constrained by mzeroe
c the velocity zerio points are stored in entry NCMAX+1
c if there is no zeros.dat file, the defaults are used
        open(unit=13,file='zeros.dat',form='formatted',status='old')
131     read(13,*,end=137)iiref,iidat,zval,eval
           if (iidat.eq.0) then
             mzero(iiref,NCMAX+1)  = zval
             mzeroe(iiref,NCMAX+1) = eval
           else 
             mzero(iiref,iidat)  = zval
             mzeroe(iiref,iidat) = eval
             endif
           go to 131
137        close(unit=13)

c these are how well data is fit by reference -- these are used as
c the error bars in the chisq if the file exists -- otherwise it
c uses the defaults from earlier
        open(unit=13,file='errors.dat',form='formatted',status='old')
130     read(13,*,end=136)iiref,iidat,eval
           if (iidat.eq.0) then
             evrms(iiref) = eval
           else
             emrms(iiref,iidat) = eval
             endif
           go to 130
136        close(unit=13)

1222    continue
c read in the Cepheid model
        nfit = 0
        ifund= 0
        iover= 0
        do kk=1,ngal
          print*,'enter input file for galaxy ',kk
          read(*,'(a30)')ifile1
          open(unit=13,file=ifile1,form='formatted',status='old')
          print*,'how many cepheids to fit '
          read(13,*)nfit0
          print*,'enter galaxy ID #, mean distance '
          read(13,*)igalt,aa,extbar1(igalt),extsig(igalt),iextgal(igalt)
          ncingal(igalt) = nfit0
          print*,'SETTING NCINGAL ',ncingal(igalt),igalt,nfit0
          nfit = nfit + nfit0
          do j=1,nfit0
            i = nfit - nfit0 + j
c debugging
	    print*,ifile1
            print*,'enter id #, dist, extinction, temp of cepheid i '
            read(13,*)iceph(i),per1(i),mass(i),dbar(i),ebar(i),tbar(i),vbar(i),rbar(i),
     1         lper(i),ramp(i),dpdt1(i),zbar0(i),zbar1(i),time0(i),imode(i),igal(i),name(i)
c experiment by increasing rbar and zero point by ln(2)
c            rbar(i) = rbar(i) + log(2.0D0)
c experiment by increasing extinction by 0.06
c            ebar(i) = ebar(i)+0.037



            write(86,'(i4,i4,i4,1x,a10)')i,igal(i),iceph(i),name(i)
c count Cepheids in different modes
            if (imode(i).eq.1) ifund = ifund + 1
            if (imode(i).eq.2) iover = iover + 1
c connect Cepheid model properties to data -- the slow way
            ifound   = 0
            do k=1,nc
              if ((igal(i).eq.idsys0(k)).and.(iceph(i).eq.idnum0(k)).and.(ifound.eq.0)) then
                iceph(i) = k
                ifound   = 1
                write(*,'(a7,i4,1x,i4,1x,a20,2(1x,f13.6))')'found: ',igal(i),iceph(i),name(iceph(i)),xx(k),yy(k)
                if ((abs(yy(k)).gt.400.0).and.(igal(i).le.3)) then
                  print*,'no coord ',igal(i),iceph(i),' ',name(iceph(i)),xx(k),yy(k)
                  if ((igal(i).eq.2).or.(igal(i).eq.3)) then
                    xx(k) = 0.0
                    yy(k) = 0.0
                    endif
                  endif
                endif
              enddo
            if (ifound.eq.0) then 
              print*,'failed to find ',igalt,igal(i),i
              stop
              endif
            per0(i)       = period(iceph(i))
c if period entry is negative reset to intrinsic period
            if (per1(i).lt.0) then
              per1(i) = per0(i)
              lper(i) = log10(per1(i)/10.0)
              endif
            lper(i)=log10(per1(i)/10.0)
c CONVERSION CODE TO RENOMALIZE TO NEW TEMPLATE FORM
            if (iconvert.eq.1) then
              call settemplate(1,lper(i))
              dphase = tw(itemp0)*phase(itemp0)+tw(itemp1)*phase(itemp1)
              print*,'shift: ',itemp0,itemp1,tw(itemp0),tw(itemp1)
              print*,'shift: ',phase(itemp0),phase(itemp1)
              print*,'shifting phase by ',dphase
              time0(i) = time0(i) + per1(i)*dphase/(2.0*pi)
              damp = tw(itemp0)*amp(itemp0)+tw(itemp1)*amp(itemp1)
              print*,'shifting amplitude by ',damp
              ramp(i) = ramp(i)*damp/ln10
              endif
c END CONVERSION CODE
c if phase reference time is negative, set to photometric minimum
            if (time0(i).lt.0) then
              ict   = iceph(i)
              apmin = 0.0 
              do ii=i0(ict),i1(ict)
                do jj=1,ndat(ii)
                  if (mag(ii,jj).gt.apmin) then
                    apmin    = mag(ii,jj)
                    time0(i) = time(ii)
                    endif
                  enddo
                enddo
              endif
c as a numerical variable us the square root of the amplitude
            ramp(i)       = sqrt(abs(ramp(i)))
            enddo
          close(unit=13)
          enddo

c set the temperature and extinction priors
        print*,'found ',nc,' Cepheids in data base '
        print*,'found ',nfit,' Cepheids to fit '

c read in the galaxy distances etc
c mean distance is fixed if idfix  = 1
c galaxy is tilted       if itcode = 1
c Cepheids in galaxy have fixed distance if idcode = 0    LMC or other
c                   individual distances           = 1    galaxy
c  individual distances constrained by range       = 2    SMC
        open(unit=13,file='distance.dat',form='formatted',status='old')
        read(13,*)ngal0
        do i=1,ngal0
          read(13,*)distbar1(i),tiltx(i),tilty(i),distsig(i),idfix(i),idcode(i),itcode(i),it,galname(i)
          enddo
        close(unit=13)

c read in the Cepheid dependence vectors and PLC priors
        print*,'reading Cepheid vectors '
        open(unit=13,file='vector.dat',form='formatted',status='old')
c Prior on period-radius relation coefficients
        read(13,*)rhoa,rhob
	read(13,*)taua,taub
c Period-radius relation
        read(13,*)(rvec(i),i=1,6),text(1)
        rvec(3) = 0.0
        rvec(4) = 0.0
        rvec(5) = 0.0
c period-temperature relation
        read(13,*)(rvec(i),i=7,11),text(2)
        rvec(9) = 0.0
        rvec(10)= 0.0
c amplitude period relation for fundamentals
        read(13,*)(rvec(i),i=12,15),text(3)
        rvec(12) = 0.0
        rvec(13) = 0.0
        rvec(14) = 0.0
c overtone/fundamental period offset (in log space)
        read(13,*)rvec(16),text(4)
        rvec(16) = 0.0
c amplitude period relation for overtones
        read(13,*)(rvec(i),i=17,20),text(5)
        rvec(17) = 0.0
        rvec(18) = 0.0
        rvec(19) = 0.0
        rvec(20) = 0.0
        print*,'warning -- setting higher order terms to zero '
        do i=1,ncol
          print*,'debug: reading vector.dat entry ',i,' of ',ncol
          read(13,*)v0(i),b0(i),r0(i),flam(i),fdlam(i),fname(i)
c experiment by increasing rbar and zero point by ln(2)
c            v0(i) = v0(i) + 5.0*log(2.0D0)
          enddo
        close(unit=13)
c CONVERSION CODE TO RENOMALIZE TO NEW TEMPLATE FORM
        if (iconvert.eq.1) then
          do i=1,ncol
            b0(i) = -b0(i)/2.5
            enddo
          endif
c END CONVERSION CODE
c set the V extinction component to be the V prior value
c this also reads the priors for the temperature vector
        print*,'setting extinction priors and force V band entry to value of prior'
        call rvprior
        r0(3) = r00(3)
	r0(2) = r00(3)+1.
c set the V component of the temperature dependence 
        b0(3) = b00(3)

c clear the flags for presence of companions, period derivatives, binaries, parallaxes
        do i=1,nfit
          ipara(i) = 0
          enddo

c convert period data into angular frequency data
        do i=1,nfit
          w0(i) = 2.0D0*pi/per1(i)
          enddo

c read in the parallax file
        print*,'reading parallax list '
        open(unit=13,file='parallax.dat',form='formatted',status='old')
        read(13,*)npara
        do i=1,npara
          read(13,*)igt,ict,tvpara,tepara
          do j=1,nfit
            if ((igt.eq.igal(j)).and.(ict.eq.idnum0(iceph(j)))) then
              print*,name(iceph(j)),' has a parallax '
              ipara(j) = 1
              vpara(j) = tvpara
              epara(j) = tepara
              endif
            enddo
          enddo
        close(unit=13)

c set up the correspondance between physical variables and the 1D numerical vector
        print*,'setting correspondance tables ncol/ncos ',ncol,ncos
        call eyeset

c set all variables to be optimized, initialize the save vector, clear the error bars
        do i=1,nvar
          pserr(i) = -10.0
          do j=1,nvar
             pmaterr(i,j) = -10.0
          enddo
          iuse(i)  = 1
          enddo
        print*,'converting physical variables to numerical vector '
        call ytop


c do not vary V value of temperature (b0) or the extinction (r0) vectors 
        icv = 3
c vary V value of temperature
        iuse(ib0+icv)  = 0
        iuse(ie0+icv)  = 0
	iuse(ie0+icv-1) = 0
c do not vary the higher order dependences of the radius on period
        iuse(irp0+3) = 0.0
        iuse(irp0+4) = 0.0
        iuse(irp0+5) = 0.0
c do not vary the higher order dependences of the temperature on period
        iuse(irp0+9) = 0.0
        iuse(irp0+10)= 0.0
c do not vary the radius  amplitude period relation
        iuse(irp0+12) = 0.0
        iuse(irp0+13) = 0.0
        iuse(irp0+14) = 0.0
c do not vary the relation of fund/overtone periods
        iuse(irp0+16) = 0.0
c do not vary the temperature amplitude period relation
        iuse(irp0+17) = 0.0
        iuse(irp0+18) = 0.0
        iuse(irp0+19) = 0.0
        iuse(irp0+20) = 0.0
c shut off galaxies not in use
        do k=1,NGMAX
          iuse(idgal0+k)  = 0
          iuse(iwgal0+k)  = 0
          iuse(itxgal0+k) = 0
          iuse(itygal0+k) = 0
          enddo
 
c leave one template coefficient in each template is degenerate with the
c amplitudes -- turn off its variation
c one template phase is degenerate with a global shift of the phase reference times
c turn off its variation
        j=1
        do km=1,2
          do ii=1,ntempl(km)
            iir = irt0(km)+(ii-1)*ntem(km)+2*j
            iuse(iir-1)   = 0
            iuse(iir  )   = 0
            enddo
          enddo

c if there are no overtone pulsators shut off associated variables
        print*,'FOUND ',ifund,' FUNDAMENTALS '
        print*,'FOUND ',iover,' OVERTONES '
        if (iover.eq.0) then
           print*,'NO OVERTONES '
           iuse(irp0+16) = 0
           iuse(irp0+17) = 0
           iuse(irp0+18) = 0
           iuse(irp0+19) = 0
           iuse(irp0+20) = 0
           do k=1,ntempl(2)
             do j=1,ncos(2)
               iit = itt0(2)+(k-1)*ntem(2)+2*j
               iir = irt0(2)+(k-1)*ntem(2)+2*j
               iuse(iit-1) = 0
               iuse(iit  ) = 0
               iuse(iir-1) = 0
               iuse(iir  ) = 0
               enddo
             enddo
           endif
c if there are no fundamental pulsators shut off associated variables
        if (ifund.eq.0) then
           print*,'NO FUNDAMENTALS '
           iuse(irp0+12) = 0
           iuse(irp0+13) = 0
           iuse(irp0+14) = 0
           iuse(irp0+15) = 0
           do k=1,ntempl(1)
             do j=1,ncos(1)
               iit = itt0(1)+(k-1)*ntem(1)+2*j
               iir = irt0(1)+(k-1)*ntem(1)+2*j
               iuse(iit-1) = 0
               iuse(iit  ) = 0
               iuse(iir-1) = 0
               iuse(iir  ) = 0
               enddo
             enddo
           endif
        

c exclude data from processing
        print*,'excluding bad and unwanted data points '
        call exclude
         
c for the Galaxy (= galaxy ID #1) turn on individual distances, turn off mean distance
c for SMC (=#2) and LMC (=#3) turn on both individual distances and mean distance
c for all other galaxies, turn on mean distance turn off individual distances
c shut of galaxy mean distance
c        iuse(idgal0+1) = 0
c turn off individual distances
c        do k=1,nfit
c          if (igal(k).gt.3) then
c            iuse(id0+k) = 0
c            endif
c          enddo
c don't vary the LMC distance '
c        iuse(idgal0+3) = 0
c set the LMC distance to 18.5
c        distbar1(3)     = 18.5
c        psave(idgal0+3) = distbar1(3)

c turn on distance variables as appropriate to their idcode
c idcode = 0 == vary only mean distance
c idcode = 1 == vary only individual distance
c idcode = 2 == vary both (linked by prior)
c turn on tilts if galaxy has tile
c itcode = 0 == no tilt
c itcode = 1 == tilt possible
c don't vary width at all '
        do k=1,ngal0
          iuse(iwgal0+k)  = 0
          iuse(idgal0+k)  = 0
          iuse(itxgal0+k) = 0
          iuse(itygal0+k) = 0
          enddo
c loop over all cepheids and turn on appropriate distance variables
        do k=1,nfit
          igt         = igal(k)
          iuse(id0+k) = 0
          if (idcode(igt).ne.0) then
            iuse(id0+k) = 1
            endif
          if (idcode(igt).ne.1) then
            iuse(idgal0+igt) = 1
            if (itcode(igt).eq.1) then
              iuse(itxgal0+igt) = 1
              iuse(itygal0+igt) = 1
              endif
            endif  
          enddo
c shut off mean distance to particular galaxies if idfix(k) = 0
        do k=1,ngal0
          if (idfix(k).eq.0) iuse(idgal0+k) = 0
          enddo


        do k=1,nfit
          igt = igal(k)
c vary either mean extinction or individual extinctions: iextgal 0=fix, 1=average, 2=individual
c at present, the mean extinction of a galaxy is just a fixed number -- it is not varied
          iuse(iegal0+igt) = 0
          iuse(iem0+k)     = 0
          if (iextgal(igt).eq.1) then
            iuse(iem0+k)     = 1
            endif
          if (iextgal(igt).eq.2) then
            iuse(iem0+k)    = 1
            endif
          enddo

  

        alam  = 2.0
        print*,'enter alam, ncuse '
        read*,alam,ncuse

        iverb = 0
        print*,'verbose (1) or not (0) '
        read*,iverb

       
        print*,'vary zero point (1/0)'
        read*,izerov
        print*,'vary temperature vector (1/0) '
        read*,itemp0
        print*,'vary extinction vector (1/0) '
        read*,iextc0


        do i=1,ncol
          if (izerov.eq.0) iuse(iv0 +i) = 0
          if (itemp0.eq.0) iuse(ib0 +i) = 0
          if (iextc0.eq.0) iuse(ie0 +i) = 0
          enddo

c count the data
        ntotv = 0
        ntotm = 0
        do j=1,ncol
          nmcol(j) = 0
          enddo
        do k=1,nfit
          ii0 = i0(iceph(k))
          ii1 = i1(iceph(k))
          tmin =  1.0E32
          tmax = -1.0E32
          nmval(k) = 0
          nvval(k) = 0
          do i=ii0,ii1
            if (vr(i).gt.-500) then
              nvval(k) = nvval(k) + 1
              ntotv    = ntotv    + 1
              tmin     = min(tmin,time(i))
              tmax     = max(tmax,time(i))
              endif
            do j=1,ndat(i)
              icol = icdat(i,j)
              if (mag(i,j).gt.-5) then
                nmval(k)    = nmval(k) + 1
                nmcol(icol) = nmcol(icol) + 1
                ntotm       = ntotm    + 1
                tmin        = min(tmin,time(i))
                tmax        = max(tmax,time(i))
                endif
              enddo
            enddo
c if there is no velocity data shut it off
          if (nvval(k).eq.0) then
            iuse(ivr0+k) = 0
            vbar(k)      = -1000.0
            psave(ivr0+k)= -1000.0
            endif
c if there is suddenly new velocity data
          if ((nvval(k).gt.0).and.(vbar(k).lt.-500)) then
            vbar(k)       = 0.0
            psave(ivr0+k) = 0.0
            endif
          nyear    = (tmax-tmin)/360.0
          sense(k) = (per0(k)/(tmax-tmin))**2/18.0
c          if ((nmval(k)+nvval(k).ge.100).and.(nyear.ge.10)) then
c            print*,'activated period derivative for ',k,iceph(k),nyear
c            idpdt1(k)       = 1
c            iuse (idpdt10+k)= 1
c            endif
          enddo

        print*,'vary the phase reference (1/0) '
          read*,ivvp0  
          if (ivvp0.eq.0) then
            print*,'  shutting off the Cepheid phase reference time '
            do i=1,nfit
              iuse(izp0+i) = 0
              enddo
            endif
        print*,'vary the individual distance (1/0) '
          read*,ivvd0
          if (ivvd0.eq.0) then
            print*,'  shutting off the individual Cepheid distances '
            do i=1,nfit
              iuse(id0+i) = 0
              enddo
            endif
        print*,'vary the mean distances (1/0) '
          read*,ivvdbar
          if (ivvdbar.eq.0) then
            print*,'  shutting off the mean galaxy distance variables '
            do k=1,NGMAX
             iuse(idgal0 +k) = 0
             iuse(itxgal0+k) = 0
             iuse(itygal0+k) = 0
             iuse(iwgal0 +k) = 0
             enddo
            endif
        print*,'vary the mean extinction (1/0) '
          read*,ivvebar
          if (ivvebar.eq.0) then
            print*,'  shutting off the mean galaxy extinctions'
            do k=1,NGMAX
              iuse(iegal0+k) = 0
              enddo
            endif
        print*,'vary the temperature (1/0) '
          read*,ivvt0
          if (ivvt0.eq.0) then
            print*,'  shutting off the Cepheid temperatures '
            do i=1,nfit
               iuse(itm0+i) = 0
              enddo 
            endif
        print*,'vary the extinction (1/0) '
          read*,ivve0
          if (ivve0.eq.0) then
            print*,'  shutting off the Cepheid extinction '
            do i=1,nfit
              iuse(iem0+i) = 0
              enddo 
            endif
        print*,'vary the radius (1/0) '
          read*,ivvr0
          if (ivvr0.eq.0) then
             print*,'  shutting off the mean radius '
            do i=1,nfit
              iuse(irr0+i) = 0
              enddo
            endif
        print*,'vary the period (1/0) '
          read*,ivp0
          if (ivp0.eq.0) then
             print*,'  shutting off the period'
            do i=1,nfit
              iuse(ip0+i) = 0
              enddo
            endif
        print*,'vary the radius amplitude (1/0) '
          read*,ivvra0
          if (ivvra0.eq.0) then
            print*,'shutting off the radius amplitude '
            do i=1,nfit
              iuse(ira0+i) = 0
              enddo
            endif
        print*,'vary the mean radial velocity (1/0) '
          read*,ivvv0
          if (ivvv0.eq.0) then
            print*,'  shutting off the mean radial velocity '
            do i=1,nfit
              iuse(ivr0+i) = 0
              enddo
            endif

        print*,'vary the template (1=yes/0=no) '
        read*,itempl
        if (itempl.eq.0) then
          print*,'  shutting off the templates '
          do km=1,2
            do k=1,ntempl(km)
              do j=1,ncos(km)
                iit = itt0(km)+(k-1)*ntem(km)+2*j
                iir = irt0(km)+(k-1)*ntem(km)+2*j
                iuse(iit-1) = 0
                iuse(iit  ) = 0
                iuse(iir-1) = 0
                iuse(iir  ) = 0
                enddo
              enddo
            enddo
          endif

        print*,'vary only the template variables (1/0) '
        read*,itempo
        if (itempo.eq.1) then
          print*,' TURNING EVERYTHING OFF BUT TEMPLATE VARIABLES '
          print*,'setting 1 to ',itt0(1),' off '
          do i=1,itt0(1)
            iuse(i) = 0
            enddo
          print*,'setting ',itt0(1)+1,' to ',irt1(2),' on '
          do i=itt0(1)+1,irt1(2)
            iuse(i) = 1
            enddo
          print*,'setting ',irt1(2)+1,' to ',nvar,' off '
          do i=irt1(2)+1,nvar
            iuse(i) = 0
            enddo
          endif

        print*,'set amplitude to zero, fit only mean properties (1/0) '
        read*,iphase
        if (iphase.eq.1) then
          print*,' setting amplitude to zero, turning amplitude phase off '
          do i=1,nfit
            iuse(ira0+i)  = 0
            psave(ira0+i) = 0.0
            ramp(i)       = 0.0
            iuse(izp0+i)  = 0
            enddo
          endif

        print*,'reinitialize temperature and radius (1/0): threshr,thresht '
        read*,init,threshr,thresht
        if (init.eq.1) then
          print*,' resetting temperatures and radii if off mean relation '
          do k=1,nfit
            ict     = iceph(k)
            tlogp   = lper0(ict)
            if (imode(k).eq.2) then
              print*,'no overtones in this version!!! '
              stop
              endif
c reset to the mean temperature
            tbar(k)       = rvec(7) + rvec(8)*tlogp + rvec(9)*tlogp*tlogp + rvec(10)*zbar1(k)
            delt          = psave(itm0+k)-tbar(k)
            if (abs(delt).gt.thresht) then
              print*,'reset temperature of ',name(ict),tbar(k),psave(itm0+k)
              psave(itm0+k) = tbar(k)
              endif
c reset to the mean radius
            rbar(k)       = rvec(1) + rvec(2)*tlogp + rvec(3)*tlogp*tlogp + rvec(4)*tbar(k) + rvec(5)*zbar1(k)
            delr          = psave(irr0+k)-rbar(k)
            if (abs(delr).gt.threshr) then
              print*,'reset radius of ',name(ict),rbar(k),psave(irr0+k)
              psave(irr0+k) = rbar(k)
              endif
            enddo
          endif
        print*,'reinitialize amplitudes (1/0) '
        read*,init
        if (init.eq.1) then
          print*,'reinitializing amplitudes '
          do k=1,nfit
            ict           = iceph(k)
            ramp0         = 0.5
            psave(ira0+k) = sqrt(ramp0)
            print*,'reinitialized amplitude of ',name(ict),ramp(k)
            enddo
          endif


        print*,'total velocity points   = ',ntotv
        print*,'total photometry points = ',ntotm

c if there is no data on some color, shut it off
        do j=1,ncol
          if (nmcol(j).eq.0) then
            iuse(iv0+j)  = 0
            iuse(ib0+j)  = 0
            iuse(ie0+j)  = 0
            iuse(ibb0+j) = 0 
            iuse(iee0+j) = 0
            iuse(izz0+j) = 0 
            iuse(iza0+j) = 0
            endif
          enddo

c fit covariance matrix  -- keep prior widths fixed
        print*,'enter itmin, icovar '
        read*,itmin,icovar
        icovaroff = 0
        if (icovar.gt.0) then
          iuse(irp0+ 6) = 0
          iuse(irp0+11) = 0
          iuse(irp0+15) = 0
          iuse(irp0+20) = 0
          endif

        print*,'shut off prior width optimization (1/0) '
        read*,ivvpw
        if (ivvpw.eq.0) then
          iuse(irp0+ 6) = 0
          iuse(irp0+11) = 0
          iuse(irp0+15) = 0
          iuse(irp0+20) = 0
          endif
        print*,'shut off priors if = 1 '
        read*,ivvpa
        if (ivvpa.eq.1) then
          do i=1,npriors
            iuse(irp0+i) = 0
            enddo
          endif

c vary calibration zero points
        print*,'vary calibration zero points (1/0) '
        read*,ical
        if (ical.eq.0) then
          do i=izpnt0+1,izpnt1
            iuse(i) = 0
            enddo
          endif

c save the iuse values
        do i=1,nvar
          iuse0(i) = iuse(i) 
          enddo

        print*,'optimize cepheids one by one (1/0) '
        read*,ione
        if (ione.eq.1) then
          do i=1,nvar
            iuse(i) = 0
            enddo
          endif

c check the derivative
         ideriv = 0
         if (ideriv.eq.1) then
           nuse = 0
           do i=1,nvar
             if (iuse(i).eq.1) nuse = nuse + 1
             write(84,'(3(1x,i5),1x,g13.6,1x,a10)')i,iuse(i),nuse,psave(i),vname(i)
             enddo
           do i=1,nvar
             iuse(i) = 1
             enddo
           do i=1,nvar
             psave(i) = 1.1*psave(i)
             enddo
           nuse = 0
           do i=1,nvar
             if (iuse(i).eq.1) then
               nuse = nuse + 1
               p(nuse) = psave(i)
               endif
             enddo
           kstart = 1
           kend   = nfit
           val0 = func(p,dp0)
           delv = 1.0D-6
c plot out chi^2 for one variable
c           do i=114,114
c             print*,'doing variable ',i,' of ',nvar
c             if (iuse0(i).eq.1) then
c               do k=1,720
c                 delvt = 0.0001/72.0
c                 delvt = 1.0/72.0
c                 p(i)  = p(i) + delvt 
c                 valp = func(p,dp1)
c                 if (k.ne.1) then
c                   grad = (valp-pold)/delvt
c                   write(70,'(i5,4(1x,g20.10),1x,i2,1x,a10)')i,p(i),valp,dp1(i),grad,iuse(i),vname(i)
c                   endif
c                 pold = valp
c                 enddo
c               endif
c             enddo
c           stop
c end plotting out chi^2 for one variable
           do i=1,nvar
             print*,'doing variable ',i,' of ',nvar
             if (iuse0(i).eq.1) then
               delvt = delv*max(1.0D0,abs(p(i)))
c for the phase time, use a fraction of a day
               if ((i.ge.izp0).and.(i.lt.izp1)) then
                 delvt = 0.1
                 endif
               p(i) = p(i) + delvt
               valp = func(p,dp1)
               p(i) = p(i) - delvt
               grad = (valp-val0)/delvt
               write(70,'(i5,5(1x,g13.6),1x,i2,1x,a10)')i,dp0(i),dp1(i),grad,p(i)/1.1,val0,iuse(i),vname(i)
               endif
             enddo
           stop
           endif

        
          ftol = 0.001


c do the minimization
          print*,'starting minimization '
          if (ione.eq.0) then
            kstart = 1
            kend   = nfit 
            nuse   = 0
            do i=1,nvar
              write(84,'(3(1x,i5),1x,g13.6,1x,a10)')i,iuse(i),nuse+1,psave(i),vname(i)
              if (iuse(i).eq.1) then
                nuse = nuse + 1
                p(nuse) = psave(i)
                endif
              enddo
            print*,'will fit ',nuse,' variables out of ',nvar
            if ((alam.gt.0.0).and.(nuse.gt.0)) then
              call frprmn(p,nuse,ftol,iter,fret,itmin)
            else
              fret = func(p,dp0)
              print*,'FUNCTION VALUE = ',fret
              endif
          else
            do k=1,nfit
              print*,'Doing cepheid ',k,iceph(k),name(k)
              kstart = k
              kend   = k
              call turnon(k)
              nuse = 0
              do i=1,nvar
                write(84,'(3(1x,i5),1x,g13.6,1x,a10)')i,iuse(i),nuse+1,psave(i),vname(i)
                if (iuse(i).eq.1) then
                  nuse = nuse + 1
                  p(nuse) = psave(i)
                  endif
                enddo
              print*,'will fit ',nuse,' variables out of ',nvar
              if ((alam.gt.0.0).and.(nuse.gt.0)) then
                call frprmn(p,nuse,ftol,iter,fret,itmin)
              else
                fret = func(p,dp0)
                endif
              call shutoff(k)
              enddo
            endif

          print*,'DONE MINIMIZATION '
c convert variable vector back to physical variables
          call ptoy


c          do i=1,1
c            fret   = func(p,dp0)
c            call ptoy
c            print*,'NEW FUNCTION VALUE = ',fret
c            enddo
c            
c          icovar = 0
c          fret   = func(p,dp0)
c          print*,'NEW FUNCTION VALUE = ',fret

c convert angular frequencies back into periods and period derivatives
          do k=1,nfit
            per1(k) = 2.0D0*pi/w0(k)
            enddo

c set up new input files by galaxy, write out summary of residuals, models
          do k=1,20
            if (k.le.9) then
              write(ofile1,'(a10,i1)')'newinput.0',k
            else
              write(ofile1,'(a9,i2)')'newinput.',k
              endif
            open(unit=110+k,file=ofile1,form='formatted',status='unknown')
            write(110+k,*)ncingal(k),' nceph ',k
            write(110+k,'(i3,3(1x,f13.6),1x,i3,2(1x,f13.6),a30)')k,distbar1(k),extbar1(k),extsig(k),
     1          iextgal(k),pserr(idgal0+k),pserr(iegal0+k)
            enddo
          print*,'AT FORT.50 WRITE NFIT = ',NFIT
          do k=1,nfit
            kk  = iceph(k)
            ii0 = i0(kk)
            ii1 = i1(kk)
c kk0 = Cepheid # for galaxy igal(k)
            kk0 = kk - igal0(igal(k))
            write(49,1110)k,kk0,igal(k),per0(k),per1(k),(per1(k)-per0(k))/per0(k)
            write(50,1115)k,kk0,igal(k),per1(k),dbar(k),ebar(k),tbar(k),
     1         vbar(k),rbar(k),lper(k),ramp(k)**2,chisq1(k),chisq2(k),
     1         sprior(k),imode(k),nmval(k),nvval(k),squote,name(k),quote
            write(51,1110)k,kk0,igal(k),pserr(ip0+k),pserr(id0+k),pserr(iem0+k),pserr(itm0+k),
     1         pserr(ivr0+k),pserr(irr0+k),pserr(ita0+k),pserr(ira0+k),pserr(idpdt10+k),pserr(izp0+k)
            write(81,1110)k,kk0,igal(k),per1(k),dpdt1(k),pserr(ip0+k),pserr(idpdt10+k),sense(k)
            write(110+igal(k),1120)kk0,per1(k),mass(k),dbar(k),ebar(k),tbar(k),vbar(k),
     1                       rbar(k),lper(k),ramp(k)**2,dpdt1(k),
     1                       zbar0(k),zbar1(k),time0(k),imode(k),igal(k),name(k), xx(kk), yy(kk), tiltx(3)*xx(kk)+tilty(3)*yy(kk)
            enddo
          do k=1,20
            close(unit=110+k)
            enddo

1110        format(i5,i5,i5,13(1x,g13.6),3(1x,i4))
1115        format(i4,1x,i4,1x,i3,11(1x,g13.6),1x,i2,2(1x,i4),a2,a10,a2)
1120        format(i5,1x,g15.8,11(1x,g13.6),1x,g15.8,1x,2(i3,1x),a10,3(1x,g13.6))
1130        format(5x,i5,4(1x,g13.6))


c read in the galaxy distances etc
        open(unit=13,file='distance.new',form='formatted',status='unknown')
        open(unit=14,file='distance.err',form='formatted',status='unknown')
        write(13,*)ngal0
        do i=1,ngal0
          write(13,'(4(1x,g13.6),4(1x,i2),1x,a1,a10,a1)')distbar1(i),tiltx(i),tilty(i),distsig(i),
     1       idfix(i),idcode(i),itcode(i),i,quote,galname(i),quote
          write(14,'(4(1x,g13.6),4(1x,i2),1x,a1,a10,a1)')pserr(idgal0+i),pserr(itxgal0+i),pserr(itygal0+i),
     1       pserr(iwgal0+i),idfix(i),idcode(i),itcode(i),i,quote,galname(i),quote
          enddo
        close(unit=13)
        close(unit=14)

c write out the priors/PLC relation etc
c write out the new Cepheid dependence vector
	open(unit=13,file='vector.new',form='formatted',status='unknown')
	  write(13,'(1x,g13.6,1x,g13.6)')rhoa,rhob
	  write(13,'(1x,g13.6,1x,g13.6)')taua,taub
          write(13,'(6(1x,g13.6),1x,a1,a30,a1)')(rvec(i),i=1,6),quote,text(1),quote
          write(13,'(5(1x,g13.6),1x,a1,a30,a1)')(rvec(i),i=7,11),quote,text(2),quote
          write(13,'(4(1x,g13.6),1x,a1,a30,a1)')(rvec(i),i=12,15),quote,text(3),quote
          write(13,'(1(1x,g13.6),1x,a1,a30,a1)')rvec(16),quote,text(4),quote
          write(13,'(4(1x,g13.6),1x,a1,a30,a1)')(rvec(i),i=17,20),quote,text(5),quote
c          write(13,'(4(1x,g13.6),1x,a1,a30,a1)')(rvec(i),i=21,24),quote,text(6),quote
          do i=1,ncol
            write(13,'(5(1x,f13.6),1x,a1,a10,a1)')v0(i),b0(i),r0(i),flam(i),fdlam(i),quote,fname(i),quote
            enddo
        close(unit=13)
        open(unit=13,file='vector.err',form='formatted',status='unknown')
          write(13,'(6(1x,g13.6))')(pserr(irp0+i),i=1,6)
          write(13,'(6(1x,g13.6))')(pserr(irp0+i),i=7,11)
          write(13,'(6(1x,g13.6))')(pserr(irp0+i),i=12,15)
          write(13,'(6(1x,g13.6))')pserr(irp0+16)
          write(13,'(6(1x,g13.6))')(pserr(irp0+i),i=17,20)
          write(13,'(6(1x,g13.6))')(pserr(irp0+i),i=21,24)
          do i=1,ncol
            write(13,'(3(1x,f13.6),1x,a10)')pserr(iv0+i),pserr(ib0+i),pserr(ie0+i),fname(i)
            enddo
        close(unit=13)

        open(unit=13,file='vector.covar',form='formatted',status='unknown')
          write(13,'(a40)')'M_V beta_V a_rho b_rho a_tau b_tau'
	  write(13,'(6(1x,g13.6))')pmaterr(iv0+3,iv0+3),pmaterr(iv0+3,ib0+3),pmaterr(iv0+3,irp0+1),
     1                             pmaterr(iv0+3,irp0+2),pmaterr(iv0+3,irp0+7),pmaterr(iv0+3,irp0+8)
	  write(13,'(6(1x,g13.6))')pmaterr(ib0+3,iv0+3),pmaterr(ib0+3,ib0+3),pmaterr(ib0+3,irp0+1),
     1                             pmaterr(ib0+3,irp0+2),pmaterr(ib0+3,irp0+7),pmaterr(ib0+3,irp0+8)
	  write(13,'(6(1x,g13.6))')pmaterr(irp0+1,iv0+3),pmaterr(irp0+1,ib0+3),pmaterr(irp0+1,irp0+1),
     1                             pmaterr(irp0+1,irp0+2),pmaterr(irp0+1,irp0+7),pmaterr(irp0+1,irp0+8)
	  write(13,'(6(1x,g13.6))')pmaterr(irp0+2,iv0+3),pmaterr(irp0+2,ib0+3),pmaterr(irp0+2,irp0+1),
     1                             pmaterr(irp0+2,irp0+2),pmaterr(irp0+2,irp0+7),pmaterr(irp0+2,irp0+8)
	  write(13,'(6(1x,g13.6))')pmaterr(irp0+7,iv0+3),pmaterr(irp0+7,ib0+3),pmaterr(irp0+7,irp0+1),
     1                             pmaterr(irp0+7,irp0+2),pmaterr(irp0+7,irp0+7),pmaterr(irp0+7,irp0+8)
	  write(13,'(6(1x,g13.6))')pmaterr(irp0+8,iv0+3),pmaterr(irp0+8,ib0+3),pmaterr(irp0+8,irp0+1),
     1                             pmaterr(irp0+8,irp0+2),pmaterr(irp0+8,irp0+7),pmaterr(irp0+8,irp0+8)
          write(13,'(6(1x,g13.6))')v0(3),b0(3),rvec(1),rvec(2),rvec(7),rvec(8)
          write(13,'(a40)')'Covariances in reddening vector'
          do i=1,ncol
            write(13,'(36(1x,g13.6))')(pmaterr(ie0+i,ie0+j),j=1,ncol)
            enddo
	  write(13,'(a40)')'Covariances between temperature and reddening, eigenvalue 1, nx1, ny1, eigenvalue 2, nx2, ny2'
           do k=1,nfit
              discriminant = sqrt( (pmaterr(itm0+k,itm0+k)-pmaterr(iem0+k,iem0+k))**2 + 4.0*pmaterr(iem0+k,itm0+k)**2 )
              eqnv1 = (pmaterr(itm0+k,itm0+k)+pmaterr(iem0+k,iem0+k) + discriminant )*0.5
              eqnv2 = (pmaterr(itm0+k,itm0+k)+pmaterr(iem0+k,iem0+k) - discriminant )*0.5
              denomdir = sqrt( pmaterr(itm0+k,iem0+k)**2 + (egnv1 - pmaterr(itm0+k,itm0+k))**2 )
              nx1 = pmaterr(iem0+k,itm0+k)/denomdir
              ny1 = (eqnv1 - pmaterr(itm0+k,itm0+k))/denomdir
              denomdir = sqrt( pmaterr(itm0+k,iem0+k)**2 + (egnv2 - pmaterr(itm0+k,itm0+k))**2 )
              nx2 = pmaterr(iem0+k,itm0+k)/denomdir
              ny2 = (eqnv2 - pmaterr(itm0+k,itm0+k))/denomdir
              write(13,'(i5,1x,i5,1x,g13.6,4(1x,g13.6))')igal(k),iceph(k)-igal0(igal(k)),pmaterr(iem0+k,itm0+k),eqnv1,ny1/nx1, eqnv2,ny2/nx2
c	      write(13,'(2(1x,g13.6))')pmaterr(itm0+k,itm0+k), pmaterr(iem0+k,itm0+k)
c              write(13,'(2(1x,g13.6))')pmaterr(itm0+k,iem0+k), pmaterr(iem0+k,iem0+k)
           enddo
        close(unit=13)


	open(unit=13,file='template.new',form='formatted',status='unknown')
        write(13,*)ncos(1),ntempl(1)
        write(13,'(38(1x,g13.6))')(tlp(1,i),i=1,ntempl(1))
c write fundamental mode templates
        do j=1,ncos(1)
          write(13,'(i5,38(1x,g13.6))')j,(ctt(1,i,j),stt(1,i,j),i=1,ntempl(1))
          enddo
        do j=1,ncos(1)
          write(13,'(i5,38(1x,g13.6))')j,(crt(1,i,j),srt(1,i,j),i=1,ntempl(1))
          enddo
c write overtone mode templates
        write(13,*)ncos(2),ntempl(2)
        write(13,'(38(1x,g13.6))')(tlp(2,i),i=1,ntempl(2))
        do j=1,ncos(2)
          write(13,'(i5,38(1x,g13.6))')j,(ctt(2,i,j),stt(2,i,j),i=1,ntempl(2))
          enddo
        do j=1,ncos(2)
          write(13,'(i5,38(1x,g13.6))')j,(crt(2,i,j),srt(2,i,j),i=1,ntempl(2))
          enddo
        close(unit=13)
        open(unit=13,file='template.err',form='formatted',status='unknown')
c write fundamental mode templates
        km = 1
        write(13,*)ncos(km),ntempl(km)
        do j=1,ncos(1)
          write(13,'(i5,38(1x,g13.6))')j,(pserr(itt0(km)+(i-1)*ntem(km)+2*j-1),
     1                                    pserr(itt0(km)+(i-1)*ntem(km)+2*j  ),i=1,ntempl(km))
          enddo
        do j=1,ncos(1)
          write(13,'(i5,38(1x,g13.6))')j,(pserr(irt0(km)+(i-1)*ntem(km)+2*j-1),
     1                                    pserr(irt0(km)+(i-1)*ntem(km)+2*j  ),i=1,ntempl(km))
          enddo
c write overtone mode templates
        km = 2
        write(13,*)ncos(km),ntempl(km)
        do j=1,ncos(2)
          write(13,'(i5,38(1x,g13.6))')j,(pserr(itt0(km)+(i-1)*ntem(km)+2*j-1),
     1                                    pserr(itt0(km)+(i-1)*ntem(km)+2*j  ),i=1,ntempl(km))
          enddo
        do j=1,ncos(2)
          write(13,'(i5,38(1x,g13.6))')j,(pserr(irt0(km)+(i-1)*ntem(km)+2*j-1),
     1                                    pserr(irt0(km)+(i-1)*ntem(km)+2*j  ),i=1,ntempl(km))
          enddo
        close(unit=13)



        print*,'write out residual etc files (1) or not (0) '
        read*,iwrite
c        if (iwrite.eq.0) then
c          stop
c          endif

c set up new input files by galaxy
          do k=1,nfit
            kk  = iceph(k)
            ii0 = i0(kk)
            ii1 = i1(kk)
c build output files
            kk0 = kk - igal0(igal(k))

            if  (kk0.le.  9)                    write(cnum,'(a4,i1)')'.000',kk0
            if ((kk0.gt.  9).and.(kk0.le.  99)) write(cnum,'(a3,i2)')'.00',kk0
            if ((kk0.gt. 99).and.(kk0.le. 999)) write(cnum,'(a2,i3)')'.0',kk0
            if ((kk0.gt.999).and.(kk0.le.9999)) write(cnum,'(a1,i4)')'.',kk0
            if  (igal(k).le.9)                      write(cgal,'(a1,i1)')'0',igal(k)
            if ((igal(k).gt.9).and.(igal(k).le.99)) write(cgal,'(i2)')igal(k)
            write(ofile1,'(a5,a2,a6,a2,a5)')'curve',cgal,'/data.',cgal,cnum
            write(ofile2,'(a5,a2,a6,a2,a5)')'curve',cgal,'/errs.',cgal,cnum
            write(ofile3,'(a5,a2,a6,a2,a5)')'curve',cgal,'/modl.',cgal,cnum
            write(ofile4,'(a5,a2,a6,a2,a5)')'curve',cgal,'/resd.',cgal,cnum
            write(ofile5,'(a5,a2,a5,a2,a5)')'curve',cgal,'/per.',cgal,cnum
            if (iwrite.eq.1) open(unit=100,file=ofile1,form='formatted',status='unknown')
            if (iwrite.eq.1) open(unit=200,file=ofile2,form='formatted',status='unknown')
            if (iwrite.eq.1) open(unit=300,file=ofile3,form='formatted',status='unknown')
            if (iwrite.eq.1) open(unit=400,file=ofile4,form='formatted',status='unknown')
	    if (iwrite.eq.1) open(unit=550,file=ofile5,form='formatted',status='unknown')
c plot the template model curves
            km      = imode(k)
            aramp   = ramp(k)
            call settemplate(km,lper(k))
c now produce the output files
            do i=1,nplt
              ang = pi*tphi(i)/180.0
              call gettemplate(0,km,k,-1,ang)
              call getvelocity(0,km,k,i,aramp,vrad(i),tvr)
              do j=1,ncol
                call getmagnitude(idoderiv,km,k,i,j,dbar(k),ebar(k),aramp,tv(j),ttot,rtot)
                enddo
              if (vbar(k).lt.-500.0) vrad(i)=vrad(i)-vbar(k)
              if (iwrite.eq.1) write(300,'(f7.3,2(1x,g13.6),1x,f9.4,37(1x,f7.3))')tphi(i),
     1            temp(i),rad(i),vrad(i),(tv(j),j=1,ncol),glog
              enddo
1210        format(i5,18(1x,g13.6))
            do i=ii0,ii1
c write out the data points shifted to the template phase - phase in units of 0 < aphi < 1 here
              delt    = time(i)-time0(k)
              aphi(i) = w0(k)*delt/(2.0*pi)
              aphi0   = aphi(i)
              iphi    = int(aphi(i))
              aphi(i) = aphi(i)-float(iphi)
              if (aphi(i).lt.0.0) aphi(i) = aphi(i) + 1.0             
              if (aphi(i).gt.1.0) aphi(i) = aphi(i) - 1.0             
              aphi(i) = 360.0*aphi(i)
              do jj=1,NCMAX
                mout(jj) = -10.0
                eout(jj) = -10.0
                enddo
              do jj=1,ndat(i)
                icol = icdat(i,jj)
                mout(icol) = mag(i,jj)
                eout(icol) = emag(i,jj)
                enddo
              if ((iwrite.eq.1).and.(mout(3).ne.-10)) write(550,*)time(i),mout(3),'V'
c              if ((iwrite.eq.1).and.(mout(20).ne.-10)) write(550,*)time(i),mout(20),'Ve'
              if ((iwrite.eq.1).and.(mout(5).ne.-10)) write(550,*)time(i),mout(5),'Ic'
c              if ((iwrite.eq.1).and.(mout(10).ne.-10)) write(550,*)time(i),mout(10),'Ij'
c              if ((iwrite.eq.1).and.(mout(22).ne.-10)) write(550,*)time(i),mout(22),'Ie'
              if ((iwrite.eq.1).and.(mout(4).ne.-10)) write(550,*)time(i),mout(4),'Rc'
c              if ((iwrite.eq.1).and.(mout(9).ne.-10)) write(550,*)time(i),mout(9),'Rj'
c              if ((iwrite.eq.1).and.(mout(21).ne.-10)) write(550,*)time(i),mout(21),'Re'

              if (iwrite.eq.1) write(100,'(f7.3,1x,f9.3,36(1x,f7.3),1x,i3)')aphi(i),
     1               vr(i),(mout(j),j=1,ncol),iref(i)
              if (iwrite.eq.1) write(200,'(f7.3,1x,f9.3,36(1x,f7.3),1x,i3)')aphi(i),
     1               evr(i),(eout(j),j=1,ncol),iref(i)
              if (vr(i).gt.-500.0) then
                if (iwrite.eq.1) write(400,'(i2,1x,i3,1x,f7.3,1x,f10.3,2(1x,f7.3),2(1x,f8.3),2(1x,f13.4))')0,iref(i),per1(k),
     1            time(i),aphi(i),vres(i),vr(i),evr(i),aphi0,delt
                endif
              do j=1,ndat(i)
                icol = icdat(i,j)
                if (mag(i,j).gt.-5.0) then
                  if (iwrite.eq.1) write(400,'(i2,1x,i3,1x,f7.3,1x,f10.3,2(1x,f7.3),2(1x,f8.3),2(1x,f13.4))')icol,iref(i),per1(k),
     1               time(i),aphi(i),res(i,j),mag(i,j),emag(i,j),aphi0,delt
                  endif
                enddo
c write out points with high residuals
              irf = iref(i)
              if (vr(i).gt.-500.0) then
                refvbar(irf) = refvbar(irf) + vres(i)
                refvres(irf) = refvres(irf) + vres(i)**2
                refvnum(irf) = refvnum(irf) + 1
c DEBUG
c                if (irf.eq.1) print*,'velocity in reference 1 ',name(kk)
c                if (irf.eq.2) print*,'velocity in reference 2 ',name(kk)
c                if (irf.eq.4) print*,'velocity in reference 4 ',name(kk)
c                if (irf.eq.9) print*,'velocity in reference 9 ',name(kk)
c                if (irf.eq.43) print*,'velocity in reference 43 ',name(kk)
c                if (irf.eq.150) print*,'velocity in reference 150 ',name(kk)
                if (abs(vres(i)).gt.2.0) then
c                  if (iwrite.eq.1) write(71,'(i3,1x,f9.3,3(1x,f7.3),1x,i4,1x,a9)')
c     1               k,time(i),aphi(i),vr(i),vres(i),iref(i),name(kk)
                   if (iwrite.eq.1) write(71,'(4(i4,1x),f9.3,1x,a9,2(1x,f7.3))')igal(k),
     1                iceph(k),0,iref(i),time(i),
     1                name(kk),vr(i),vres(i)
                  endif
                endif
              do jj=1,ndat(i)
                icol = icdat(i,jj)
                if (mag(i,jj).gt.-5.0) then
                  refpbar(irf,icol) = refpbar(irf,icol) + res(i,jj)
                  refpres(irf,icol) = refpres(irf,icol) + res(i,jj)**2
                  refpnum(irf,icol) = refpnum(irf,icol) + 1
                  if ((emag(i,jj).lt.0.01).or.(emag(i,jj).gt.1.0)) then
                    if (iwrite.eq.1) write(80,'(4(i4,1x),f9.3,1x,a9,2(1x,f7.3))')igal(k),
     1                iceph(k),icol,iref(i),time(i),
     1                name(kk),aphi(i),emag(i,jj)
                    endif
                  if (abs(res(i,jj)).gt.0.20) then
c  if (iwrite.eq.1) write(70,'(i3,1x,i2,1x,f9.3,3(1x,f7.3),1x,i4,1x,a9)')k,icol,
c     1               time(i),aphi(i),mag(i,jj),res(i,jj),iref(i),name(kk)
                      if (iwrite.eq.1) write(70,'(4(i4,1x),f9.3,1x,f9.3,1x,a9,2(1x,f7.3))')igal(k),
     1                   iceph(k),icol,iref(i),time(i),aphi(i)/360.,
     1                   name(kk),mag(i,jj),res(i,jj)
                    endif
                  endif
                enddo
              enddo
            close(unit=100)
            close(unit=200)
            close(unit=300)
            close(unit=400)
            enddo

c this records the errors for every reference by band/velocity -- if copied
c to errors.dat, these are then used at the errors for the chisquare
        open(unit=13,file='errors.new',form='formatted',status='unknown')
        do i=1,NREFMAX
          if (refvnum(i).gt.0) then
            evbar(i) = refvbar(i)/float(refvnum(i))
            evrms(i) = refvres(i)/float(refvnum(i))
            evrms(i) = sqrt(abs(evrms(i)-evbar(i)**2))
            write(73,'(i5,1x,i5,2(1x,f13.6),1x,i6,1x,f6.3)')i,0,evrms(i),evbar(i),refvnum(i),-1.0
            write(13,'(i5,1x,i5,2(1x,f13.6),1x,i6,1x,f6.3)')i,0,evrms(i),evbar(i),refvnum(i),-1.0
            endif
          enddo
        do i=1,NREFMAX
          do j=1,NCMAX
            if (refpnum(i,j).gt.0) then
              embar(i,j) = refpbar(i,j)/float(refpnum(i,j))
              emrms(i,j) = refpres(i,j)/float(refpnum(i,j))
              emrms(i,j) = sqrt(abs(emrms(i,j)-embar(i,j)**2))
              write(74,'(i5,1x,i5,2(1x,f13.6),1x,i6,1x,f6.3)')i,j,emrms(i,j),embar(i,j),refpnum(i,j),flam(j)
              write(13,'(i5,1x,i5,2(1x,f13.6),1x,i6,1x,f6.3)')i,j,emrms(i,j),embar(i,j),refpnum(i,j),flam(j)
              endif
            enddo
          enddo
        close(unit=13)

c write out the updated estimates of the zero points for each reference
        open(unit=13,file='zeros.new',form='formatted',status='unknown')
        do i=1,NREFMAX
          if (refvnum(i).gt.0) then
            write(13,'(2(1x,i5),2(1x,f13.6))')i,0,mzero(i,NCMAX+1),mzeroe(i,NCMAX+1)
            endif
          do j=1,NCMAX
            if (refpnum(i,j).gt.0) then
              write(13,'(2(1x,i5),2(1x,f13.6))')i,j,mzero(i,j),mzeroe(i,j)
              endif
            enddo
          enddo
        close(unit=13)

        end



