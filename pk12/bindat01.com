        

        parameter (NCMAX=36,NPMAX=8900000,NFMAX=12000)
        parameter (NVMAX=20000,NBMAX=50)
        parameter (NGMAX=20)
        parameter (NREFMAX=600)
        parameter (NCOVAR=4000)
c        parameter (NCOVAR=1800)
c        parameter (NCOVAR=3000)





c INPUT DATA
c input cepheid period
        real period(NFMAX)
c log of input period
        real lper0(NFMAX)
c per0(i) = period(iceph(i))
        real per0(NFMAX)
c coordinates
        real xx(NFMAX),yy(NFMAX)
c reference epoch
        real time0(NFMAX)
c phase of observation
        real aphi(NPMAX)
c time of observation
        real time(NPMAX)
c number of colors available at time i
        integer*2 ndat(NPMAX)
c correspondance between time i, entry j and color = icdat(i,j)
        integer*2 icdat(NPMAX,5)
c A -10 ENTRY IS USED FOR EMPTY OR EXCLUDED MAGNITUDES
c mag(i,j)  = mag at epoch time(i), in color icdat(i,j)
c emag(i,j) = errors in mag at epoch time(i), in color icdat(i,j)
c res(i,j)  = model residual at time(i), in color icdat(i,j)
        real mag(NPMAX,5),res(NPMAX,5),emag(NPMAX,5)
c A -1000 ENTRY IS USED FOR EMPTY OR EXCLUDED VELOCITIES
c vr(i)  = radial velocity at time(i)
c evr(i) = error in velocity at time i
c vres(i)= model residual in velocity at time(i)
        real vr(NPMAX),evr(NPMAX),vres(NPMAX)
c the reference from which data point i comes is iref(i)
        integer*2 iref(NPMAX)
c Cepheid ID for data point i is ic(i)
        integer*2 ic(NPMAX)
c total number of input Cepheids, total number of data points
        integer nc,nd
c data on Cepheid i starts at i0(i) and ends at i1(i)
        integer i0(NFMAX),i1(NFMAX)
c character name of Cepheid
        character*9  name(NFMAX)
c idnum0(i) = Cepheid id # in parent galaxy for i-th Cepheid being fit
c idsys0(i) = galaxy id # for i-th cepheid being fit
        integer idnum0(NFMAX),idsys0(NFMAX)

        common /data1a/period,lper0,per0,xx,yy,aphi,time,time0
        common /data1b/mag,emag,res
        common /data1c/vr,evr,vres
        common /data1d/icdat,ndat
        common /data1e/iref
        common /data1f/i0,i1,nc,nd
        common /data1g/name
        common /data1h/idnum0,idsys0


c galaxy ID number: cepheid i is in galaxy igal(i)
        integer igal(NFMAX)
c Companion flag: cepheid i has a companion if icom(i)=1
        integer icom(NFMAX)
c Template ID#:  cepheid is between template idtm(i) and idtm(i)+1
        integer idtm(NFMAX)
c Binary flag: cepheid i is in a binary if ibin(i)=1
        integer ibin(NFMAX)
c Period Derivative flag: cepheid i has a period derivative if idpdt1(i)=1
        integer idpdt1(NFMAX)
c Period Second Derivative flag: cepheid i has period 2nd deriv if idpdt2(i)=1
        integer idpdt2(NFMAX)
c mode flag: 1=fundamental, 2=overtone
        integer imode(NFMAX)
c the number of data points on Cepheid i is npt(i)
        integer npt(NFMAX)
        common /data3/igal,icom,ibin,idpdt1,idpdt2,imode,idtm,npt


c Indices for transforming from minimization vector p back to physical variables
        common /eyes/iv0,iv1,ib0,ib1,ie0,ie1,id0,id1,
     1       itm0,itm1,iem0,iem1,ivr0,ivr1,irr0,irr1,
     1       it0,it1,ir0,ir1,ibb0,ibb1,iee0,iee1,ig00,ig01,ig10,ig11,
     1       ita0,ita1,ira0,ira1,izp0,izp1,itt0(2),itt1(2),
     1       irt0(2),irt1(2),ip0,ip1,ipt0,ipt1,irp0,irp1,izz0,izz1,iza0,iza1,
     1       izc0,izc1,idgal0,idgal1,icomt0,icomt1,icomr0,icomr1,im0,im1,
     1       ibint0,ibint1,ibine0,ibine1,ibinv0,ibinv1,ibinp0,ibinp1,
     1       ibinw0,ibinw1,idpdt10,idpdt11,idpdt20,idpdt21,iegal0,iegal1,
     1       izpnt0,izpnt1,itxgal0,itxgal1,itygal0,itygal1,iwgal0,iwgal1
c Operate on a restrictred range of Cepheids from kstart to kend
        common /data2/kstart,kend


c dependence vectors
c filter names
        character*10 fname(NCMAX)
c filter wavelengths and widths
        real*8 flam(NCMAX),fdlam(NCMAX)
        common /filters1/fname
        common /filters2/flam,fdlam
c Zero point vector 
        real*8 v0(NCMAX)
c prior estimated value v0, and its uncertainty 
        real*8 v00(NCMAX),ev00(NCMAX)
c Extinction: constant part r0 temperature dependent part r1
        real*8 r0(NCMAX),r1(NCMAX)
c prior estimated value r0, and its uncertainty 
        real*8 r00(NCMAX),er00(NCMAX)
c prior estimated value r1, and its uncertainty 
        real*8 r10(NCMAX),er10(NCMAX)
c Temperature: constant part b0, temperature dependent part b1
        real*8 b0(NCMAX),b1(NCMAX)
c prior estimated value b0, and its uncertainty 
        real*8 b00(NCMAX),eb00(NCMAX)
c prior estimated value b1, and its uncertainty 
        real*8 b10(NCMAX),eb10(NCMAX)
c Log(g): constant part g0, temperature dependent part g1
        real*8 g0(NCMAX),g1(NCMAX)
c prior estimated value of g0 and its uncertainty
        real*8 g00(NCMAX),eg00(NCMAX)
c prior estimated value of g1 and its uncertainty
        real*8 g10(NCMAX),eg10(NCMAX)
c Composition: zero point correction z0, temperature correction z1
        real*8 z0(NCMAX),z1(NCMAX)
c prior estimated value z0, and its uncertainty 
        real*8 z00(NCMAX),ez00(NCMAX)
c prior estimated value z1, and its uncertainty 
        real*8 z10(NCMAX),ez10(NCMAX)
        common /model1/v0,r0,r1,b0,b1,z0,z1,g0,g1
        common /model2a/v00,r00,r10,b00,b10,z00,z10,g00,g10
        common /model2b/ev00,er00,er10,eb00,eb10,ez00,ez10,eg00,eg10
c Prior limits on radius and temperature == PLC/PL relations
        integer npriors
        real taua, taub, rhoa, rhob
        real*8 rvec(50)
        common /model5/rvec,npriors,taua,taub,rhoa,rhob
c Kurucz stellar spectra models -- third order polynomial fit as a function
c of T expanded around 5500 K -- kurz0 is difference of zeropoint from model
c with only linear terms 
        real*8 kurz0(NCMAX),kurz1(NCMAX),kurz2(NCMAX),kurz3(NCMAX)
        common /kurucz/kurz0,kurz1,kurz2,kurz3

         
c BINARY VARIABLES
c   number of binaries
        integer nbinary
c   reference time, period, ellipticity, K, line of nodes
      real*8 bint(NBMAX),binp(NBMAX),bine(NBMAX),binv(NBMAX),binw(NBMAX)
c   prior constraints on period and ellipticity from other sources
        real*8 binp0(NBMAX),bine0(NBMAX)
c   real binary or just an acceleration
        integer btype(NBMAX)
        common /binary1/bint,binp,bine,binv,binw,binp0,bine0,btype,nbinary


c PARALLAX MEASUREMENTS
c   number of parallaxes = npara
c   ipara(i)=1 means there is a parallax for Cepheid i
        integer ipara(NFMAX),npara
c   parallax value and its error for Cepheid i
        real*8 vpara(NFMAX),epara(NFMAX)
        common /parallax/vpara,epara,ipara,npara


c Statistical Weight factors
c   ntotv = total number of velocities, ntotm = total number of photometric points
c   used to give equal weight to velocity and photometry data in chi^2
        integer ntotv,ntotm
        common /weights/ntotv,ntotm

c individual mean cepheid properties: 
c    distance, extinction, temperature
        real*8 dbar(NFMAX),ebar(NFMAX),tbar(NFMAX)
c    radial velocity, radius, period
        real*8 vbar(NFMAX),rbar(NFMAX),per1(NFMAX)
c mean temperature amplitude, radius amplitude phase zero point
        real*8 tamp(NFMAX),ramp(NFMAX)
c composition prior and estimated value
        real*8 zbar0(NFMAX),zbar1(NFMAX)
c radius and temperature of companion
        real*8 comr(NFMAX),comt(NFMAX)
c period derivative
        real*8 dpdt1(NFMAX)
c period second derivative
        real*8 dpdt2(NFMAX)
c angular frequencdy and its first and second derivatives
        real*8 w0(NFMAX),w1(NFMAX),w2(NFMAX)
c log period
        real*8 lper(NFMAX)
c log mass
        real*8 mass(NFMAX)
c icarry=1 --> use Cepheid only in individual fits
        integer icarry(NFMAX)
        common /model3/dbar,ebar,tbar,vbar,rbar,tamp,ramp,per1,
     1     zbar0,zbar1,comr,comt,dpdt1,dpdt2,w0,w1,w2,lper,mass,icarry


c individual cepheids being fit
c   ncol = total number of colors being fit
c   ntem = total number of variables associted with radius or temperature template
c   nfit = total number of Cepheids to be fit
c   iceph(i) = cepheid ID # for the i-th cepheid being fit (e.g. its name is name(iceph(i)) )
        integer iceph(NFMAX)
c   nuse = total number of variables to be fit
c   nvar = total number of variables in model
c   iverb = verbose mode if =1
c   ncos  = Order of Fourier series in phase for template
c   ntempl= Order of polynomial series in period for template
        integer ncos(2),ntempl(2),ntem(2)

        common /model4/iceph,ncol,ntem,nfit,nuse,nvar,
     1                    iverb,ncos,ntempl
        
c if ideriv=1 then the program is systematically checking the derivative expressions
        common /debug1/ideriv

c temperature and radius templates
c   ctt=cosine part of temperature, stt=sine part of temperature
c   crt=cosine part of radius     , srt=sine part of radius
c 1=fundamental, 2=overtone
        real*8 ctt(2,20,20),stt(2,20,20)
        real*8 crt(2,20,20),srt(2,20,20)
c log(P/10) for template i of mode k = tlp(k,i) 
        real*8 tlp(2,20)
        common /template1/ctt,stt,crt,srt,tlp

c Properties of individual galaxies
c   number of galaxies currently being fit
        integer ngal
c   number of galaxies in distance table
        integer ngal0
c   number of first Cepheid in galaxy, number of cepheids in that galaxy
        integer igal0(NGMAX),ncingal(NGMAX)
c   ID # of galaxy k = ignum(k)
c   input number of galaxy with ID #k = ignum0(k)
        integer ignum(NGMAX),ignum0(NGMAX)
c   mean distance to that galaxy, mean extinction
        real*8  distbar1(NGMAX),extbar1(NGMAX) 
c   tilt of galaxy
        real*8 tiltx(NGMAX),tilty(NGMAX)
c   width of galaxy
        real*8 distsig(NGMAX)
c   extinction spread of a galaxy
        real*8 extsig(NGMAX)
c   name of galaxy
        character*20 galname(NGMAX)
c   code for how to treat distances to galaxy i: 
c     idcode = 0 at fixed distance
c     idcode = 1 at individual distance
c     idcode = 2 at individual distance constrained to a range around fixed distance
        integer idcode(NGMAX)
c   keep mean distance fixed if idfix(i)=1
        integer idfix(NGMAX)
c   code for whether galaxy is tilted
c      itcode = 0 no tilt
c      itcode = 1 tilted
         integer itcode(NGMAX)
c   vary mean distance/extinction or individual distance/extinction
        integer idistgal(NGMAX),iextgal(NGMAX)
        common /galaxies/distbar1,extbar1,idistgal,iextgal,igal0,ignum,ignum0,ngal
        common /galaxies2/tiltx,tilty,distsig,extsig,galname,idcode,itcode,idfix,ngal0


c Minimization variables: handles transfer of physical variables in and out of
c the minimization vector
c   the 1D vector corresponding to the current model
        real*8 psave(NVMAX)
c   variable names
        character*10 vname(NVMAX)
c   the 1D vector of the errors in the parameters
        real*8 pserr(NVMAX)
c   matrix of covariances
        real*8 pmaterr(NVMAX,NVMAX)
c   the estimated first derivative of the chi^2 with respect to variable i
c   the Marquandt method ``second derivative'' of chi^2 with respect to i
        real*8 dps(NVMAX)
        real*8 dqs(NVMAX)
c   if iuse(i)=1 then entry i in psave(i) will be optimized, otherwise it is held fixed
c   iuse0 = original model of what is on and off
        integer iuse(NVMAX),iuse0(NVMAX)
        real*8 alam
c   flag to build covariance matrix
        integer icovar
        common /fit1/psave,pserr,pmaterr,dps,dqs,iuse,iuse0,alam,icovar
        common /fit2/vname
c   the number of magnitude and velocity measurements for Cepheid i
        integer nmval(NFMAX),nvval(NFMAX)
c   prior probability values
        real*8 sprior(NFMAX),rprior(NFMAX)
c   the velocity and magnitude chi^2 values
        real*8 chisq1(NFMAX),chisq2(NFMAX)
        common /fitqual/sprior,rprior,chisq1,chisq2,nmval,nvval

c useful constants: ln(10) and pi
        real*8 ln10,pi
c expansion velocity factor, conversion of template derivative to velocity
        real*8 pexp,vcon
c gravity: r0/t0^2=acon, mass r0^3/G/t0^2=mcon
        real*8 acon,mcon
        common /constants/ln10,pi,pexp,vcon,acon,mcon


c errors in velocity: 
c   mean error for reference i = evbar(i)
c   mean error for reference i in color j = embar(i,j)
c   rms error for reference i = evrms(i)
c   rms error for reference i in color j = emrms(i,j)
c   zero point shift for reference i = mzero(i,j)
c   uncertainty in zero point for reference i = mzeroe(i,j)
c      velocity in NCMAX+1 entry
c   variable code for zero point = icode(i,j)
        real*8  mzero(NREFMAX,NCMAX+1),mzeroe(NREFMAX,NCMAX+1)
        integer izero(NREFMAX,NCMAX+1)
        real*8  evrms(NREFMAX),emrms(NREFMAX,NCMAX)
        real*8  evbar(NREFMAX),embar(NREFMAX,NCMAX)
        common  /errors/evrms,emrms,evbar,embar
        common  /zeros/mzero,mzeroe,izero


c index of minimized variable to full variable vector
        integer ptops(NVMAX)
c index of variable vector to minimized vector
        integer pstop(NVMAX)
        common /copies/ptops,pstop

c covariance matrix
        real*8 covar(NCOVAR,NCOVAR),covari(NCOVAR,NCOVAR),diag(NCOVAR)
        common /covarmat/covar,covari,diag
c shut off covariance evaluation if it won't be used
        common /nocovar/icovaroff


c template related variables
c   rms amplitudes of radius and temperature templates
        real*8 rmsr,rmst
c   derivative of rms amplitude with respect to period
        real*8 drmsr,drmst
c   derivatives of rmsr and rmst with respect to template coefficients
        real*8 rdelc(20),rdels(20),tdelc(20),tdels(20)
c   storage of cosines and sines
        real*8 cost(20),sint(20)
c   template weights with period, range of period index used for current cepheid
c   derivative of template weights with period
        real*8 tw(20),dtw(20)
        integer itemp0,itemp1
c   radius/temperature template value and its first three phase derivatives
        real*8 rt(4),tt(4)
c   period derivative radius/temperature template value and its first three phase derivatives
        real*8 drt(4),dtt(4)
c   time derivative of phase
        real*8 dphidt1,dphidt2
c   derivatives of temperature, radius and velocity with period
        real*8 dtp(5),drp(5),dvp(5),dap(5)
c   derivatives of temperature, radius and velocity with phase
        real*8 dta(5),dra(5),dva(5),daa(5)
c   derivatives of temperature, radius and velocity with period derivative
c        real*8 dtb(2),drb(2),dvb(2),dab(2)
c   derivatives of temperature, radius and velocity with period second derivative
c        real*8 dtc(2),drc(2),dvc(2),dac(2)
c   radius = 10^rbar(k)
        real*8 rad0
        common /templatevars1/rdelc,rdels,cost,sint,tdelc,tdels,tw,dtw,rt,tt,drt,dtt
        common /templatevars2/itemp0,itemp1
        common /templatevars3/dtp,drp,dvp,dap,dta,dra,dva,daa,dtb,drb,dvb,dab,dtc,drc,dvc,dac,dphidt1,dphidt2
        common /templatevars4/rmsr,rmst,drmsr,drmst,rad0
           
