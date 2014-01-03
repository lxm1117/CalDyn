      implicit double precision(a-h,o-z)
      parameter (nux=9000,nc=2,nn=2)
      dimension  a(nux),b(nux),cc(nux),r(nux),gam(nux),u(nux),
     c   ub(nux),
     c   vol(nux),tnum(nux),area(nux),open1(nux),open2(nux),
     c   boun1(nux),boun2(nux),
     c   desn11(nux),desn21(nux),desns1(nux),
     c   desn12(nux),desn22(nux),desns2(nux),
     c   boun21(nux),
     c   boun22(nux),
     c   rel(90),
     c   b1lc(nux),blc(nux),brc(nux),b1rc(nux),up(nux),
     c   ar1(nux),ar2(nux),a2r1(nux),a2r2(nux),a2d1(nux),a2d2(nux),
     c   a2r1s(nux),a2r2s(nux),
     c   ad1(nux),a2d1s(nux),
     c   ad2(nux),a2d2s(nux),
     c   r1(nux),r2(nux),r1t(nux),r2t(nux),
     c   alf(2),be(2),akoff(2),akon(2),akoff2(2),akon2(2),
     c   akd1(2),akd2(2),akr1(2),akr2(2),
     c   be1(2),be2(2),be3(2),be4(2),
     c   alf1(2),alf2(2),alf3(2),alf4(2),
     c   akon3(2),akoff3(2),
     c   pR(nux)
c        dad(nux),ddad(nux),da2ds(nux),dda2ds(nux)
c
c
      character*40 filout,fillst,ufile
c
      open(file='inputfilesh20',unit=39,status='old')
 1    read(39,59,end=9999) fillst
 59   format(a40)
c
      ir=32
      iw=31
      iu=34
      inul=33
      open(unit=inul,status='unknown',file='syn2.null')
c
c      nul file contains the write statements useful for debugginhg
c
      open(unit=ir,status='old',file=fillst)
c
c
       write(inul,*) ' enter non-NMDA',
     c   'alpha,beta,koff,kon,koff2,kon2,koff3,kon3,alpha1,beta1,alpha2,
     c   beta2,alpha3,beta3,alpha4,beta4'
       read(ir,*) alf(1),be(1),akoff(1),akon(1),akoff2(1),akon2(1),
     c    akoff3(1),akon3(1),alf1(1),be1(1),alf2(1),
     c    be2(1),alf3(1),be3(1),alf4(1),be4(1)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
       write(inul,*) ' enter NMDA',
     c   'alpha,beta,koff,kon,koff2,kon2,koff3,kon3,alpha1,beta1,alpha2,
     c   beta2,alpha3,beta3,alpha4,beta4'
       read(ir,*) alf(2),be(2),akoff(2),akon(2),akoff2(2),akon2(2),
     c    akoff3(2),akon3(2),alf1(2),be1(2),alf2(2),
     c    be2(2),alf3(2),be3(2),alf4(2),be4(2)
c
      write(inul,*) ' enter output filename '
      read(ir,22) filout
 22   format(a40)
      write(inul,*) ' output file name is ',filout
      open(unit=iw,status='unknown',file=filout)
c
      write(inul,*) ' enter ufile  filename '
      read(ir,22) ufile

c        ufile contians glu concentration at various locations
c
      write(inul,*) ' output file name is ',ufile
c     open(unit=iu,status='unknown',file=ufile)
c
      write(inul,*) ' enter max time point '
      read(ir,*) tmax
      write(inul,*) ' maxt= ',tmax
c
      isw=1
      pi=4.d0*datan(1.d0)
      write(inul,*) ' enter dt .001 is default '
      read(ir,*) dt
      write(inul,*) ' dt = ',dt
      dr=.001d0
      write(inul,*) ' enter dr .001 is default '
      read(ir,*) dr
      write(inul,*) ' dr = ',dr
      nt=tmax/dt
      write(inul,*) ' enter print interval for first 1 ms, rest '
      read(ir,*) tprt1,tprt2
      write(inul,*) ' tprt1,tprt2= ',tprt1,tprt2
      prt1=tprt1/dt+1.d-4
      prt2=tprt2/dt+1.d-4
      iprt=prt1
      write(inul,*) ' nt= ',nt
c23456789012345678901234567890123456789012345678901234567890123456789012
      n=4.0d0/dr +.01
      write(inul,*) ' n= ',n
      if(n.GT.9000) then
        write(6,*) ' overflow of n = ',n
        goto 9998
      endif
c23456789012345678901234567890123456789012345678901234567890123456789012
      d=6.0d-1
      write(inul,*) ' enter diffusion constant, i.e. .6 '
      read(ir,*) d
      write(inul,*) ' d= ',d
      clefth=0.02d0
      write(inul,*) ' enter cleft height '
      read(ir,*) clefth
      write(inul,*) ' enter uptake end '
      read(ir,*) uend
      write(inul,*) ' uend = ',uend
c     u0=0.d0
      v=.025d0
c
c       v is vesicle radius
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c     nvt2=0.5d0*vt/dt
c     nvt=2*nvt2
      nvr=v/dr+.01
      nu=0.5/dr+.01
      nuend=uend/dr+.01
      write(inul,*) ' nu = ',nu
      write(inul,*) ' nuend = ',nuend
      nvr=int(v/dr+0.1d0) +1
c
c    why repeat nvr with +1? 
c
c     nvr1 = nvr +1
      write(inul,*) ' enter psd radius, typically psdr=0.1d0 '
      read(ir,*) psdr
      npsd=psdr/dr +.01
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      gluuu=0.0d0
c     gluuu=2.9d0
      write(inul,*) ' enter gluu '
      read(ir,*) gluuu
c
      write(inul,*) ' enter number of NT molecules '
      read(ir,*) nmolnt
      write(inul,*) ' # NT molecules is ',nmolnt
      write(inul,*) ' enter receptor densities--non-NMDA, NMDA '
      write(inul,*) ' density*pi*psdr*psdr =  #   receptors '
      read(ir,*) r1den,r2den
      write(inul,*) ' non-NMDA density= ',r1den,' NMDA= ',r2den
      write(inul,*) ' enter uptake Vmax and Km '
      read(ir,*) vmax,aku
      read(ir,*) vmax2,aku2
      read(ir,*) vmax3,aku3
      alfa=dfloat(nmolnt)*1.d-3*1.5273d+6 
c        1.5273d+6 scales to 1000 molecules
c        1.d+6 gives 655 molecules of NT with r=.001
c        and the number scales linearly  
c        alfa = 3*nmolnt/(pi*vesicle_radius**2) from Wathey et al
c        initial glu is considered to be contained within the 
c        vesicle radius as a density. alfa = #/um2
c
      write(inul,*) ' enter number of time points to release nt '
      read(ir,*) nrel
      write(inul,*) nrel
      write(inul,*) ' enter ',nrel, ' time points '
      read(ir,*) (rel(i),i=1,nrel)
      write(inul,98) (rel(i),i=1,nrel)
 98   format(8f7.0)
      do 5 i=1,n
c            # per um2 / cleft height = # per volume in um3
c            #/um3  x 1 mole/6.02e23 x 10e12 um3/cm3 = moles/ml
c            moles/ml x 1000 ml/l = M
c            M x 1e6 = uM
c          So den/um2  *  (1/ch  /6.02e23   *1e12   *1e3  *1e6 )
c          = den/cleft height /602
         area(i)=2.d0*dfloat(i-1)*dr*dr*pi
         vol(i)=clefth*area(i)                 
c
c          values for i=1 are fixed later
C
         ub(i)=0.d0
c
c         ub not used currently, so set to 0
c
         rov=dfloat(i-1)*dr/v
         if(rov.le.1.d0) then
            prop=rov**4-2.d0*rov**2+1.d0
c                ^^^^^^^^^^^^^^^^^^^^^^^^
c                  is 1 at 0 and 0 at rov=1
            u(i)=alfa*prop/clefth /602.d0
            tnum(i)=alfa*prop*2.d0*dfloat(i-1)*pi*dr*dr
            up(i)=0.d0
         else
            u(i)=0.d-0
            up(i)=(uden/clefth) /602.d0
c
c             up not used currently, uden not defined, so 0
c
         endif
         if(i.gt.npsd) u(i)=gluuu
 5    continue
      do 6 i=1,npsd
         r1(i)= (r1den/clefth) /602.d0
         r1t(i)=r1(i)
         r2(i)=(r2den/clefth) /602.d0
         r2t(i)=r2(i)
         ar1(i)=0.d0
         ar2(i)=0.d0
         a2r1(i)=0.d0
         a2r2(i)=0.d0
         ad1(i)=0.d0
         ad2(i)=0.d0
         a2d1(i)=0.d0
         a2d1s(i)=0.d0
         a2d2s(i)=0.d0
         a2d2(i)=0.d0
         a2r1s(i)=0.d0
         a2r2s(i)=0.d0
 6    continue
      anumto=0.d0
      tr1=0.d0
      tr2=0.d0
c
c        fixing section 1 area and volume
c
      area(1)=pi*dr*dr/4.d0
      vol(1)=clefth*area(1)
      tnum(1)=alfa*pi*dr*dr/4.d0
      do  62 ii=1,npsd
            tr1=tr1+r1(ii)*vol(ii)*602.d0
            tr2=tr2+r2(ii)*vol(ii)*602.d0
            anumto=anumto+tnum(ii)
            write(inul,*) tnum(ii)
            if (ii.lt.30) write(6,61) ii,u(ii),r1(ii)
 61         format(i4,2f12.5)
 62      continue
         write(6,*) ' tot nt= ',anumto,nmolnt, ' should be equal '
c     read(6,*) fff
c
      k=0
      time=dfloat(k)
           sumnt=anumto
              sum1r=r1den*pi*psdr*psdr
              sum2r=r2den*pi*psdr*psdr
           sumntt=nmolnt
           sumbnt=0.d0
           sum1b=0.d0
           sum2b=0.d0
           sum1b2=0.d0
           sum2b2=0.d0
           sum1o=0.d0
           sum2o=0.d0
           sum1d1=0.d0
           sum2d1=0.d0
           sum1d2=0.d0
           sum2d2=0.d0
           sum1ds=0.d0
           sum2ds=0.d0
      write(6,*) tr1,tr2,sum1r,sum2r
           write(iw,131) time,sumnt,sumntt,sumbnt,
     c          sum1r,sum1b,sum1b2,sum1d1,sum1d2,sum1ds,sum1o,
     c          sum2r,sum2b,sum2b2,sum2d1,sum2d2,sum2ds,sum2o
      write(iu,66) time,u(1),u(30),u(250),u(400),u(1000),u(1999)
c
      icount=0
      ak=.1d0
      akuu=.1d0
      write(inul,*) ' enter uptake rate constants '
      read(ir,*) ak,akuu
c     close(unit=inul)
      dtinv=1.d0/dt
      dr2=d/(dr*dr)
c
C
      do 900 i=1,n
         u02k=ak*up(i)*0.5d0

         b1lc(i)=dtinv+2.d0*dr2+u02k
         blc(i)=dtinv+1.d0*dr2+u02k
         brc(i)=dtinv-1.d0*dr2-u02k
         b1rc(i)=dtinv-2.d0*dr2-u02k
 900  continue
      close(unit=inul)
c
c              main loop
c
      do 1000 k=1,nt
c
         time=dfloat(k)*dt
c        dtt=dt
c        write(6,*) ' time= ',time
            if(isw.eq.1) then
               dtt=0.5d0*dt
               isw=0
            else
               dtt=dt
            endif
c
         do 40 i=1,nrel
            ubrel=rel(i)+dtt*1.1d0
            ulbrel=rel(i)+dtt*0.9d0
            if(time.gt.ulbrel.and.time.lt.ubrel) irsw=1
 40      continue
         if(irsw.eq.1) then
            do 50 i=1,n
               rov=dfloat(i-1)*dr/v
               if(rov.le.1.d0) then
                  prop=rov**4-2.d0*rov**2+1.d0
                  u(i)=u(i)+alfa*prop/clefth/602.d0
               endif
 50         continue
            irsw=0
         endif
c
c   u(i) gives glu concentration at locations i
c
c                 advance receptor conc by dt/2 
c
c
         do 1100 i=1,npsd
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
          R(i) = r1t(i)-a2r1(i)-a2r1s(i)-a2d1(i)-ar1(i)-ad1(i)-a2d1s(i)
c
c              R gives receptor concentration here
c
            dar = dtt*(ad1(i)*alf1(1)
     c               + a2r1(i)*akoff2(1)
     c               + u(i)*R(i)*akon(1)
     c               - ar1(i)*(u(i)*akon2(1)+akoff(1)+be1(1)))
c
            da2r = dtt*(a2d1(i)*alf2(1)
     c                + a2r1s(i)*alf(1)
     c                + ar1(i)*u(i)*akon2(1)
     c                - a2r1(i)*(akoff2(1)+be(1)+be2(1)))
c
            da2rs = dtt*(a2d1s(i)*alf3(1)
     c                 + a2r1(i)*be(1)
     c                 -a2r1s(i)*(alf(1)+be3(1)))
c
            dad = dtt*(a2d1(i)*akoff3(1)
     c               + ar1(i)*be1(1)
     c               - ad1(i)*(u(i)*akon3(1)+alf1(1)))
c
            da2d = dtt*(a2d1s(i)*alf4(1)
     c                + a2r1(i)*be2(1)
     c                + ad1(i)*u(i)*akon3(1)
     c                - a2d1(i)*(akoff3(1)+alf2(1)+be4(1)))
c
            da2ds = dtt*(a2r1s(i)*be3(1)
     c                 + a2d1(i)*be4(1)
     c                 - a2d1s(i)*(alf3(1)+alf4(1)))
c
            par = ar1(i) + dar
            pa2r = a2r1(i) + da2r
            pa2rs = a2r1s(i) + da2rs
            pad = ad1(i) + dad
            pa2d = a2d1(i) + da2d
            pa2ds = a2d1s(i) + da2ds
c
c
            pR(i) = r1t(i)-par-pa2r-pa2rs-pad-pa2d-pa2ds
c
            ddar = dtt*(pad*alf1(1)
     c               + pa2r*akoff2(1)
     c               + u(i)*pR(i)*akon(1)
     c               - par*(u(i)*akon2(1)+akoff(1)+be1(1)))
c
            dda2r = dtt*(pa2d*alf2(1)
     c                + pa2rs*alf(1)
     c                + par*u(i)*akon2(1)
     c                - pa2r*(akoff2(1)+be(1)+be2(1)))
c
            dda2rs = dtt*(pa2ds*alf3(1)
     c                 + pa2r*be(1)
     c                 - pa2rs*(alf(1)+be3(1)))
c
            ddad = dtt*(pa2d*akoff3(1)
     c               + par*be1(1)
     c               - pad*(u(i)*akon3(1)+alf1(1)))
c
            dda2d = dtt*(pa2ds*alf4(1)
     c                + pa2r*be2(1)
     c                + pad*u(i)*akon3(1)
     c                - pa2d*(akoff3(1)+alf2(1)+be4(1)))
c
            dda2ds = dtt*(pa2rs*be3(1)
     c                 + pa2d*be4(1)
     c                 - pa2ds*(alf3(1)+alf4(1)))
c
            ar1(i)=ar1(i)+0.5d0*(dar+ddar)
            a2r1(i)=a2r1(i)+0.5d0*(da2r+dda2r)
            a2r1s(i)=a2r1s(i)+0.5d0*(da2rs+dda2rs)
            ad1(i)=ad1(i)+0.5d0*(dad+ddad)
            a2d1(i)=a2d1(i)+0.5d0*(da2d+dda2d)
            a2d1s(i)=a2d1s(i)+0.5d0*(da2ds+dda2ds)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
          R(i) = r2t(i)-a2r2(i)-a2r2s(i)-a2d2(i)-ar2(i)-ad2(i)-a2d2s(i)
c
c              R gives receptor concentration here
c
            dar = dtt*(ad2(i)*alf1(2)
     c               + a2r2(i)*akoff2(2)
     c               + u(i)*R(i)*akon(2)
     c               - ar2(i)*(u(i)*akon2(2)+akoff(2)+be1(2)))
c
            da2r = dtt*(a2d2(i)*alf2(2)
     c                + a2r2s(i)*alf(2)
     c                + ar2(i)*u(i)*akon2(2)
     c                - a2r2(i)*(akoff2(2)+be(2)+be2(2)))
c
            da2rs = dtt*(a2d2s(i)*alf3(2)
     c                 + a2r2(i)*be(2)
     c                 -a2r2s(i)*(alf(2)+be3(2)))
c
            dad = dtt*(a2d2(i)*akoff3(2)
     c               + ar2(i)*be1(2)
     c               - ad2(i)*(u(i)*akon3(2)+alf1(2)))
c
            da2d = dtt*(a2d2s(i)*alf4(2)
     c                + a2r2(i)*be2(2)
     c                + ad2(i)*u(i)*akon3(2)
     c                - a2d2(i)*(akoff3(2)+alf2(2)+be4(2)))
c
            da2ds = dtt*(a2r2s(i)*be3(2)
     c                 + a2d2(i)*be4(2)
     c                 - a2d2s(i)*(alf3(2)+alf4(2)))
c
            par = ar2(i) + dar
            pa2r = a2r2(i) + da2r
            pa2rs = a2r2s(i) + da2rs
            pad = ad2(i) + dad
            pa2d = a2d2(i) + da2d
            pa2ds = a2d2s(i) + da2ds
c
c
            pR(i) = r2t(i)-par-pa2r-pa2rs-pad-pa2d-pa2ds
c
            ddar = dtt*(pad*alf1(2)
     c               + pa2r*akoff2(2)
     c               + u(i)*pR(i)*akon(2)
     c               - par*(u(i)*akon2(2)+akoff(2)+be1(2)))
c
            dda2r = dtt*(pa2d*alf2(2)
     c                + pa2rs*alf(2)
     c                + par*u(i)*akon2(2)
     c                - pa2r*(akoff2(2)+be(2)+be2(2)))
c
            dda2rs = dtt*(pa2ds*alf3(2)
     c                 + pa2r*be(2)
     c                 - pa2rs*(alf(2)+be3(2)))
c
            ddad = dtt*(pa2d*akoff3(2)
     c               + par*be1(2)
     c               - pad*(u(i)*akon3(2)+alf1(2)))
c
            dda2d = dtt*(pa2ds*alf4(2)
     c                + pa2r*be2(2)
     c                + pad*u(i)*akon3(2)
     c                - pa2d*(akoff3(2)+alf2(2)+be4(2)))
c
            dda2ds = dtt*(pa2rs*be3(2)
     c                 + pa2d*be4(2)
     c                 - pa2ds*(alf3(2)+alf4(2)))
c
            ar2(i)=ar2(i)+0.5d0*(dar+ddar)
            a2r2(i)=a2r2(i)+0.5d0*(da2r+dda2r)
            a2r2s(i)=a2r2s(i)+0.5d0*(da2rs+dda2rs)
            ad2(i)=ad2(i)+0.5d0*(dad+ddad)
            a2d2(i)=a2d2(i)+0.5d0*(da2d+dda2d)
            a2d2s(i)=a2d2s(i)+0.5d0*(da2ds+dda2ds)
cccccccc
c          new receptor concentrations
c
       r1(i)=r1t(i)-ar1(i)-a2r1(i)-a2r1s(i)-ad1(i)-a2d1(i)-a2d1s(i)
c
       r2(i)=r2t(i)-ar2(i)-a2r2(i)-a2r2s(i)-ad2(i)-a2d2(i)-a2d2s(i)
c           if(r1(i).ne.0.d0.and.i.gt.npsd) 
c    c         write(6,*) ar1(i),a2r1s(i),ar2(i),a2r2s(i)
c           write(6,*) r1(i),r2(i)
c
 1100    continue
c
c       do 1200 i=1,n
c          dub=dtt*(ak*u(i)*(up(i)-ub(i))-akuu*ub(i))
c          pub=ub(i)+dub
c          ddub=dtt*(ak*u(i)*(up(i)-pub)-akuu*pub)
c          ub(i)=ub(i)+0.5d0*(dub+ddub)
c1200   continue
c
c                 load tridiag matrix
c
         i=1
         a(1)=0.d0
         b(1)=b1lc(1)+0.5d0*(-ak*ub(1)+
     c      akon(1)*r1(1)+akon(2)*r2(1)+
     c      akon2(1)*ar1(1)+akon2(2)*ar2(1)+
     c      akon3(1)*ad1(1)+akon3(2)*ad2(1))
         cc(1)=-2.d0*dr2
         r(1)=(b1rc(1)+0.5d0*(ak*ub(1)-
     c          akon(1)*r1(1)-akon(2)*r2(1)-
     c          akon2(1)*ar1(1)-akon2(2)*ar2(1)-
     c          akon3(1)*ad1(1)-akon3(2)*ad2(1)))*u(1) -
     c        cc(1)*u(2)+
     c        akoff(1)*ar1(1)+akoff(2)*ar2(1)+
     c        akoff2(1)*a2r1(1)+akoff2(2)*a2r2(1)+
     c        akoff3(1)*a2d1(1)+akoff3(2)*a2d2(1)
c           write(6,86) i,a(1),b(1),cc(1),r(1)
 86         format(i4,4e12.5)
c
         do  97 i=2,nvr
            fi=dfloat(i-1)
            a(i)  = dr2*0.5d0*(0.5d0/fi-1.d0)
            b(i)  = blc(i)+0.5d0*(-ak*ub(i)+akon(1)*r1(i)+
     c         akon(2)*r2(i)+akon2(1)*ar1(i)+akon2(2)*ar2(i)+
     c         akon3(1)*ad1(i)+akon3(2)*ad2(i))
            cc(i)  =-dr2*0.5d0*(0.5d0/fi+1.d0)
            r(i)  =-a(i)*u(i-1)+
     c              (brc(i)+0.5d0*(ak*ub(i)-
     c          akon(1)*r1(i)-akon(2)*r2(i)-
     c          akon2(1)*ar1(i)-akon2(2)*ar2(i)-
     c          akon3(1)*ad1(i)-akon3(2)*ad2(i)))*u(i)-
     c              cc(i)*u(i+1)+
     c        akoff(1)*ar1(i)+akoff(2)*ar2(i)+
     c        akoff2(1)*a2r1(i)+akoff2(2)*a2r2(i)+
     c        akoff3(1)*a2d1(i)+akoff3(2)*a2d2(i)
c           write(6,86) i, a(i),b(i),cc(i),r(i)
  97      continue
c
         do 100 i=nvr+1,npsd
            fi=dfloat(i-1)
            a(i)  = dr2*0.5d0*(0.5d0/fi-1.d0)
            b(i)  = blc(i)+0.5d0*(-ak*ub(i)+akon(1)*r1(i)+
     c         akon(2)*r2(i)+akon2(1)*ar1(i)+akon2(2)*ar2(i)+
     c         akon3(1)*ad1(i)+akon3(2)*ad2(i))
            cc(i)  =-dr2*0.5d0*(0.5d0/fi+1.d0)
            r(i)  =-a(i)*u(i-1)+
     c              (brc(i)+0.5d0*(ak*ub(i)-
     c          akon(1)*r1(i)-akon(2)*r2(i)-
     c          akon2(1)*ar1(i)-akon2(2)*ar2(i)-
     c          akon3(1)*ad1(i)-akon3(2)*ad2(i)))*u(i)-
     c              cc(i)*u(i+1)+
     c        akoff(1)*ar1(i)+akoff(2)*ar2(i)+
     c        akoff2(1)*a2r1(i)+akoff2(2)*a2r2(i)+
     c        akoff3(1)*a2d1(i)+akoff3(2)*a2d2(i)
     c         -vmax*u(i)/(aku+u(i))
c           write(6,86) i, a(i),b(i),cc(i),r(i)
 100      continue
c
         do 120 i=npsd+1,nu
            fi=dfloat(i-1)
            a(i)  = dr2*0.5d0*(0.5d0/fi-1.d0)
            b(i)  = blc(i)+0.5d0*(-ak*ub(i))
            cc(i)  =-dr2*0.5d0*(0.5d0/fi+1.d0)
            r(i)  =-a(i)*u(i-1)+
     c              (brc(i)+0.5d0*(ak*ub(i)))*u(i)-
     c              cc(i)*u(i+1)
     c         -vmax*u(i)/(aku+u(i))
c           write(6,86) i, a(i),b(i),cc(i),r(i)
 120      continue
c
         do 130 i=nu+1,nuend
            fi=dfloat(i-1)
            a(i)  = dr2*0.5d0*(0.5d0/fi-1.d0)
            b(i)  = blc(i)+0.5d0*(-ak*ub(i))
            cc(i)  =-dr2*0.5d0*(0.5d0/fi+1.d0)
            r(i)  =-a(i)*u(i-1)+
     c              (brc(i)+0.5d0*(ak*ub(i)))*u(i)-
     c              cc(i)*u(i+1)
     c         -vmax2*u(i)/(aku2+u(i))
c           write(6,86) i, a(i),b(i),cc(i),r(i)
 130      continue
c
c
         do 150 i=nuend+1,n-1
            fi=dfloat(i-1)
            a(i)  = dr2*0.5d0*(0.5d0/fi-1.d0)
            b(i)  = blc(i)+0.5d0*(-ak*ub(i))
            cc(i)  =-dr2*0.5d0*(0.5d0/fi+1.d0)
            r(i)  =-a(i)*u(i-1)+
     c              (brc(i)+0.5d0*(ak*ub(i)))*u(i)-
     c              cc(i)*u(i+1)
     c         -vmax3*u(i)/(aku3+u(i))
c           write(6,86) i, a(i),b(i),cc(i),r(i)
 150      continue
c
          a(n)  = dr2*0.5d0*(0.5d0/dfloat(n-1)-1.d0)
          b(n)  = blc(n)+0.5d0*(-ak*ub(n))
          cc(n)  = 0.d0
          r(n)  =-a(n)*u(n-1)+
     c            (brc(n)+0.5d0*(ak*ub(n)))*u(n)
     c         -vmax3*u(n)/(aku3 +u(n))
     c         +gluuu*(+dr2*0.5d0*(0.5d0/dfloat(n)+1.d0))*2.d0
c                write(6,86) i,a(n),b(n),cc(n),r(n)
c
c                solve tridiag system
c
         IF(B(1).EQ.0.)PAUSE
         BET=B(1)
         U(1)=R(1)/BET
         DO 11 J=2,N
           GAM(J)=CC(J-1)/BET
           BET=B(J)-A(J)*GAM(J)
           IF(BET.EQ.0.)PAUSE
           U(J)=(R(J)-A(J)*U(J-1))/BET
11       CONTINUE
         DO 12 J=N-1,1,-1
           U(J)=U(J)-GAM(J+1)*U(J+1)
12       CONTINUE
c
c
        time=dfloat(k)*dt
        icount=icount+1
        if(icount.ge.iprt) then
c          write(iu,66) time,u(1),u(30),u(250),u(400),  
c    c        u(1000),u(1999)    
 66        format(f8.4,6e11.4)
           icount=0
           sumnt=0.d0
           sum1r=0.d0
           sum2r=0.d0
           sum1o=0.d0
           sum2o=0.d0
           sum1b=0.d0
           sum2b=0.d0
           sum1b2=0.d0
           sum2b2=0.d0
           sum1d1=0.d0
           sum2d1=0.d0
           sum1d2=0.d0
           sum2d2=0.d0
           sum1ds=0.d0
           sum2ds=0.d0
           do 2000 i=1,npsd
              tadd=u(i)*clefth*602.d0*area(i)
c             tdif=(tnum(i)-tadd)
c             write(inul,*) ' tdif= ',tdif
              sumnt=sumnt+u(i)*clefth*602.d0*area(i)
              sum1r=sum1r+r1(i)*clefth*602.d0*area(i)
              sum2r=sum2r+r2(i)*clefth*602.d0*area(i)
              open1(i)=a2r1s(i)*clefth*602.d0*area(i)
              open2(i)=a2r2s(i)*clefth*602.d0*area(i)
              desn11(i)=ad1(i)*clefth*602.d0*area(i)
              desn12(i)=ad2(i)*clefth*602.d0*area(i)
              desn21(i)=a2d1(i)*clefth*602.d0*area(i)
              desn22(i)=a2d2(i)*clefth*602.d0*area(i)
              desns1(i)=a2d1s(i)*clefth*602.d0*area(i)
              desns2(i)=a2d2s(i)*clefth*602.d0*area(i)
              boun1(i)=ar1(i)*clefth*602.d0*area(i)
              boun2(i)=ar2(i)*clefth*602.d0*area(i)
              boun21(i)=a2r1(i)*clefth*602.d0*area(i)
              boun22(i)=a2r2(i)*clefth*602.d0*area(i)
              sum1o=sum1o+open1(i)
              sum2o=sum2o+open2(i)
              sum1b=sum1b+boun1(i)
              sum2b=sum2b+boun2(i)
              sum1b2=sum1b2+boun21(i)
              sum2b2=sum2b2+boun22(i)
              sum1d1=sum1d1+desn11(i)
              sum2d1=sum2d1+desn12(i)
              sum1d2=sum1d2+desn21(i)
              sum2d2=sum2d2+desn22(i)
              sum1ds=sum1ds+desns1(i)
              sum2ds=sum2ds+desns2(i)
 2000      continue
c                write(6,86) i,a(n),b(n),cc(n),r(n)
           zumnt=0.d0
c          do 3000 ik=npsd+1,n
c             zumnt=zumnt+u(ik)*clefth*602.d0*area(ik)
c3000      continue
c          write(inul,*) zumnt
c23456789012345678901234567890123456789012345678901234567890123456789012
           sumbnt=sum1o+sum2o+sum1b+sum2b+
     c           sum1d1+sum2d1+sum2d2+sum1d2+sum1ds+sum2ds+
     c           sum1b2+sum2b2
           sumntt=sumnt+sumbnt
c          write(6,88) time,sumnt,sumntt,sumbnt
 88        format(f8.4,3f9.3)
c          write(6,66) time,sum1r,sum2r,      
c    c         sum1b,sum2b,sum1o,sum2o
           write(iw,131) time,sumnt,sumntt,sumbnt,
     c          sum1r,sum1b,sum1b2,sum1d1,sum1d2,sum1ds,sum1o,
     c          sum2r,sum2b,sum2b2,sum2d1,sum2d2,sum2ds,sum2o
 131       format(f9.4,3f9.2,2(3f9.3,4f8.3))
c23456789012345678901234567890123456789012345678901234567890123456789012
           if(time.ge.1.d0) iprt=prt2
        endif
 1000 continue
 9998 close(unit=ir)
      close(unit=iw)
c     close(unit=inul)
c     close(unit=iu)
      goto 1
 9999 continue
      stop
      end
