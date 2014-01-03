        PROGRAM DIFFUS
	implicit double precision(a-h,o-z)
c
	COMMON /PATH/KMAX,KOUNT,DXSAV,XP(3000),
     c	 pi,ct1,ct2,bt1,bt2,cf(4),cb(4),bb,bf,y0,
     c   ca(2),t(2),col(14),fff,zfvol,tend,
     c   ncol,icc,ica,iu,irepet
c
	dimension a(60,60),al(60),d(60),ar(60),v(60),tn(60),r(60),
     c	 p(60),b(60,60)
c
	ir=34
	id=39
	open(file='dif0504br.dat',unit=ir,status='old')
c
c	write(6,*) ' enter number of compartments '
	read(ir,*) n
	nvar=7*n+14
        write(6,*) ' nvar= ',nvar
c
c       write(6,*) ' enter itholo, isu ')
 1      read(ir,*) itholo,isu
c
	call maine(n,nvar,a,b,al,d,ar,v,tn,r,p,itholo,isu,
     c    id,ir)
c
c	write(6,*) ' any more files 1=yes '
	read(ir,*) nnn
	if(nnn.eq.1) goto 1
	stop
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine maine(n,nvar,a,b,al,d,ar,v,tn,r,p,itholo,isu,
     c    id,ir)
	implicit double precision(a-h,o-z)
c
 	EXTERNAL DERIVS,rnd
c       parameter (itholo=166,isu=6)
c
	COMMON /PATH/KMAX,KOUNT,DXSAV,XP(3000),
     c	 pi,ct1,ct2,bt1,bt2,cf(4),cb(4),bb,bf,y0,
     c   ca(2),t(2),col(14),fff,zfvol,tend,
     c   ncol,icc,ica,iu,irepet
c
	DIMENSION a(n,n),al(n),d(n),ar(n),v(n),tn(n),r(n),p(n),
     c	   y(nvar),dydx(nvar),ystart(nvar),
     c	   yscal(nvar),YOUT(nvar),yt(nvar),dyt(nvar),
     c	   dym(nvar),ytemp(nvar),ysav(nvar),dysav(nvar),
     c     b(n,n),
     c     ap(4),am(4),bp(4),bm(4),cp(4),cm(4),
     c     yp(5),ym(5),zp(5),zm(5),
     c	   im(16), zk(12), pr(12),icamk(itholo,isu+1,2),
     c	   icamb(itholo), icamtp(itholo), icama(itholo),
     c     icamcp(itholo),iphos(itholo),icamf(itholo),
     c	   ipo4(3*isu+1), itp(isu+1),ibd(isu+1),iau(isu+1),
     c     icp(isu+1),ifr(isu+1),iThis(5000)
c
	character*64 filen(100)
        character*40 file2,file3,fileic
	pi=4.d0*datan(1.d0)
	iu=31
c
 1	do 100 i=1,n
	   do 100 j=1,n
	      a(i,j)=0.d0
 100	continue
c
	call geta(n,al,d,a,tn,r,ar,v,ir,id,lsh,lsn,isd,
     c		dif,cnvrt)
        amol=cnvrt
c
c            calmodulin diffusion is 10% that of calcium
c       write(6,*) ' enter calmodulin diffusion const '
        read(ir,*) caldif
  	do 105 i=1,n
	   do 105 j=1,n
	      b(i,j)=caldif/dif*a(i,j)
 105	continue
c
  	X=0.d0
c	      write(6,*) ' enter end time point '
	read(ir,*) x2
        tend=x2
c	      write(6,*) ' enter 5 compartments to monitor '
	read(ir,*) (im(k),k=1,5)
	   im(6)=im(1)+n
	   im(7)=im(1)+2*n
	   im(8)=im(1)+3*n
	   im(9)=im(1)+4*n
	   im(10)=im(2)+n
	   im(11)=im(2)+2*n
	   im(12)=im(2)+3*n
	   im(13)=im(2)+4*n
	   im(14)=im(1)+5*n
	   im(15)=im(2)+5*n
	write(6,88) im
 88	format(12I6)
c	     write(6,*) ' enter which compartment gets the current '
	read(ir,*) ica
	kmax=3000
	     write(6,*) ' enter # filenames '
        read(ir,*) ifilen
        do 89 i=1,ifilen
	   read(ir,101) filen(i)
           write(6,101) filen(i)
   89   continue
  101	format(a64)
	    write(6,*) ' enter # of columns after time '
	read(ir,*) ncol
c	    write(6,*) ' which col has calcium in it '
	read(ir,*) icc
c	    write(6,*) ' enter ca current multiplier '
	read(ir,*) fff
c
c	    write(6,*) ' enter calmodulin buffer ct1 and ct2 '
	read(ir,*) ct1,ct2
c	    write(6,*) ' enter fast buffer bt1 and bt2 '
	read(ir,*) bt1,bt2
c	    write(6,*) ' enter cf,cb   rate const for calmodulin - Ca'
	read(ir,*) cf(1),cf(2),cf(3),cf(4),cb(1),cb(2),cb(3),cb(4)
c	    write(6,*) ' enter bf,bb   rate const for  fast buffer - Ca'
	read(ir,*) bf,bb
 	    write(6,*) ' enter pump velocity 1.4 x 10-3 um/ms '
	read(ir,*) pv
	do 400 i=1,n
	   p(i)=pv*2.0d-3/r(i)
 400	continue
	p(1)=p(1)+pv*1.0d-3/al(1)
	p(lsh)=p(lsh)+pv*(1.0d-3/al(lsh))*(1.d0-(r(lsh+1)/r(lsh))**2)
	p(lsn+1)=p(lsn+1)+pv*1.0d-3/al(lsn+1)
	p(n)=p(n)+pv*1.0d-3/al(n)
c	write(6,*) p
	y0=.02d0
c
c           write(6,*) ' enter pka, I1-total, PP1-total '
        read(ir,*) pka, ai1t, pp1t
        write(6,*) ' pka, ai1t, pp1t= ',pka, ai1t, pp1t
c           write(6,*) ' enter k2i1, kmi1, k2i1p, kmi1p, kp1pp, km1pp '
        read(ir,*) ak2i1, akmi1, ak2i1p, akmi1p, akp1pp, akm1pp 
        write(6,*) ' k2i1,kmi1,k2i1p,kmi1p,kp1-pp,km1-pp '
        write(6,*)  ak2i1, akmi1, ak2i1p, akmi1p, akp1pp, akm1pp 
c           write(6,*) ' enter k2pp, ppkm '
        read(ir,*) ak2pp, ppkm
        write(6,*) ak2pp, ppkm

c
	DO 3 I=1,n
	  YSTART(I)=0.02d0
 3	CONTINUE
 	    write(6,*) ' enter buffer factor '
	read(ir,*) bffac
c
        x31=ct1/(1.d0+cf(4)/cb(4)*2.0d-2+cb(3)/(cf(3)*2.0d-2)+
     c     cb(2)/cf(2)*cb(3)/(cf(3)*4.0d-4)+
     c     cb(1)/cf(1)*cb(2)/cf(2)*cb(3)/cf(3)/8.0d-6)
        x21=x31*cb(3)/(cf(3)*2.0d-2)
        x11=x21*cb(2)/(cf(2)*2.0d-2)
        x01=x11*cb(1)/(cf(1)*2.0d-2)
        x41=ct1-x31-x21-x11-x01
c
        x32=ct2/(1.d0+cf(4)/cb(4)*2.0d-2+cb(3)/(cf(3)*2.0d-2)+
     c     cb(2)/cf(2)*cb(3)/(cf(3)*4.0d-4)+
     c     cb(1)/cf(1)*cb(2)/cf(2)*cb(3)/cf(3)/8.0d-6)
        x22=x32*cb(3)/(cf(3)*2.0d-2)
        x12=x22*cb(2)/(cf(2)*2.0d-2)
        x02=x12*cb(1)/(cf(1)*2.0d-2)
        x42=ct2-x32-x22-x12-x02
        write(6,*) x01,x11,x21,x31,x41
        write(6,*) x02,x12,x22,x32,x42
	YSTART(n+1)=x01
	YSTART(n+2)=x02
	DO 4 I=n+3,n+n
	  YSTART(I)=x02
 4	CONTINUE
        do 7 i=2*n+1,3*n
           ystart(i)=ystart(i-n)*cf(1)/cb(1)*2.0d-2
           ystart(i+n)=ystart(i-n)*cf(1)*cf(2)*4.0d-4/
     c        (cb(1)*cb(2))
           ystart(i+2*n)=ystart(i-n)*8.0d-6*cf(1)*cf(2)*cf(3)/
     c        (cb(1)*cb(2)*cb(3))
           ystart(i+3*n)=ystart(i-n)*16.0d-8*cf(1)*cf(2)*cf(3)*cf(4)/
     c        (cb(1)*cb(2)*cb(3)*cb(4))
           ystart(i+4*n)=bt2/(1.d0+bf/bb*2.0d-2)
 7      continue 
c
        ystart(1+6*n)=bt1/(1.d0+bf/bb*2.0d-2)
        do 8 i=7*n+1,nvar
           ystart(i)=0.d0
 8      continue
        ystart(7*n+12)=0
        qb=pp1t+akm1pp/akp1pp-ai1t
        qc=-akm1pp/akp1pp*ai1t
        qdisc=qb*qb-4.d0*qc
        if(qdisc.lt.0) write(6,*) 'bad concs ic' 
	qu=0.5*(-qb+dsqrt(qdisc))
        write(6,*) qb,qc,qu
        if(qu.gt.ai1t) qu=0.5*(-qb-dsqrt(qb*qb-4.d0*qc))
        if(qu.lt.0) write(6,*) ' neg concs ic'
        write(6,*) qb,qc,qu
        ystart(7*n+13)=qu  
        ystart(7*n+14)=akm1pp*pp1t/(akp1pp*qu+akm1pp)
c
        write(6,*) ystart(7*n+13),ystart(7*n+14)
       write(6,*) ystart(n+1),ystart(2*n+1),ystart(3*n+1),
     c    ystart(4*n+1),ystart(5*n+1),ystart(6*n+1),ystart(7*n+1)
c
c***  time loop replacing former  ODEint, RKqc, RK4  portion of program
c***
	Do 11 ii= 1,nvar
	  y(ii)=ystart(ii)
 11	continue
c
	bb41 =bt1-y(6*n+1)	
	bb42 =bt2-y(6*n+2)	
        pp1ina=pp1t-y(7*n+14)
c		write(6,*) ' enter output filename '	
	     read(ir,102) file2
 102         format(a40)
	     open(file=file2,unit=36,status='unknown')
c
	iend=0
	Kount=0
       do 198 iTx=1,5000
          iThis(iTx)=0
 198   continue
c
        iTcan=10
        iTcanb=0
       write(6,*) ' enter calcineurin conc '
       read(ir,*) can
       cnt=can
       iTcan=can/cnvrt
c
	zk(1)=0.1  
	zk(2)=.0025     
	zk(3)=.005
	zk(4)=1d-5
	zk(5)=.00003
	zk(6)=.005
	zk(7)=.00001
	zk(8)=.00003
	zk(9)=.033
	zk(10)=.000015
        zk(11)=.1
        zk(12)=.0006
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       1=  0->b       5= a->0        9= a->t
c       2=  b->0       6= a->c       10= c->a
c       3=  b->t       7= c->0       11= bind to CaN
c       4=  t->a       8= t->b       12= unbind from CaN
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(6,*) ' enter zk values '
       read(ir,*) (zk(ll),ll=1,12)
       write(6,*) zk
       zk4=zk(4)
       zk5=zk(5)
       zk7=zk(7)
       zk8=zk(8)
       zk10=zk(10)
c
c      write(6,*) ' enter ap '
       read(ir,*) ap
c      write(6,*) ' enter am '
       read(ir,*) am
c      write(6,*) ' enter bp '
       read(ir,*) bp
c      write(6,*) ' enter bm '
       read(ir,*) bm
c      write(6,*) ' enter cp '
       read(ir,*) cp
c      write(6,*) ' enter cm '
       read(ir,*) cm
c      write(6,*) ' enter yp '
       read(ir,*) yp
c      write(6,*) ' enter ym '
       read(ir,*) ym
c      write(6,*) ' enter zp '
       read(ir,*) zp
c      write(6,*) ' enter zm '
       read(ir,*) zm
c
	autom=1.0
c
c
	dcal4=0.0
c
        iTcamf=itholo*isu
        ckt=dfloat(iTcamf)*cnvrt
 	iTcamb=0
 	iTcamtp=0
        iTcama=0
        iTcamcp=0
	Do 110 ih=1, iTholo
            icamf(ih)=isu
	    icamb(ih)=0
	    icamtp(ih)=0
            icama(ih)=0
            icamcp(ih)=0
	    Do 110 i=1, isu+1
	       icamk(ih,i,2)=0
	       icamk(ih,i,1)=0
 110	continue
c
c           A+C = 8%  A=8 and C=80
c
c       iTcama=8
c       iTcamcp=80
c       write(6,*) ' enter iTcama, iTcamcp
        read(ir,*) iTcama, iTcamcp
c       iTcamf=iTcamf-iTcama-iTcamcp
        iTTT=0
        if(iTcama.gt.0) iTTT=1
        if(iTTT.eq.1) then
           iTcamf=0
           iTcamb=0
           iTcama=0
           iTcamtp=0
           iTcamcp=0
c          write(6,*) ' enter IC file name '
           read(ir,102) fileic
           open(unit=71,file=fileic,status='old')
           do 82 i=1,itholo
              icamf(i)=0
              icamb(i)=0
              icamtp(i)=0
              icama(i)=0
              icamcp(i)=0
              read (71,72) (icamk(i,iss,2),iss=2,isu+1)
              do 83 j=2,isu+1
                 icamk(i,j,1)=icamk(i,j,2)
                 if(icamk(i,j,2).eq.0) then
                    iTcamf=iTcamf+1
                    icamf(i)=icamf(i)+1
                 endif
                 if(icamk(i,j,2).eq.1) then
                    iTcamb=iTcamb+1
                    icamb(i)=icamb(i)+1
                 endif
                 if(icamk(i,j,2).eq.2) then
                    iTcamtp=iTcamtp+1
                    icamtp(i)=icamtp(i)+1
                 endif
                 if(icamk(i,j,2).eq.3) then
                    iTcama=iTcama+1
                    icama(i)=icama(i)+1
                 endif
                 if(icamk(i,j,2).eq.4) then
                    iTcamcp=iTcamcp+1
                    icamcp(i)=icamcp(i)+1
                 endif

 83           continue
              icamk(i,1,1)=icamk(i,isu+1,1)
              icamk(i,1,2)=icamk(i,isu+1,2)
 82        continue
        endif
c
c
	x=0.d0
	H= 0.0005d0
        write(6,*) ' enter dt .0005 '
        read(ir,*) H
        write(6,*) H
	maxstp=(x2/H + 1.d-4)
        ca(2)=0.d0
        t(2)=0.d0
        nrepet=1
        write(6,*) ' enter nrepet, final mult value '
        read(ir,*) nrepet,nextra
c
c       pintvl=1.d0
        pintvl=1.d0/H
	camkintvl=.5d0/H
        write(6,*) ' enter camk interval, binding interval '
        read(ir,*) civl,bivl
        camkintvl=civl/H
        bintvl=bivl/H
        bintvl=1.0
        write(6,*) ' camkintvl, bintvl ',camkintvl,bintvl
        write(6,*) ' enter ampf, iago '
        read(ir,*) ampf,ago
        write(6,*)  ampf, ago 
        iago=int(ago)
        if(iago.gt.5000) write(6,*) ' iago too big '
        
c************************ repeat loop   *************
c
       do 17 irepet=1,nrepet
	open(file=filen(irepet),unit=iu,status='old')
	read(iu,*) ty,(col(i),i=1,ncol)
        tt=ty+(irepet-1)*tend
	cc=col(icc)*fff
c
c  1e-12C/sec*1e6(umol/mol)/(zF*vol(pi*.275**2*.05e-15liters)*1e3ms/sec)
c
	zfvol=2.d0*9.648456d4*ar(ica)*al(ica)
c       cc=cc/(2.d0*9.648456d4*pi*r(1)*r(1)*al(1))*1d6
	cc=cc/zfvol*1d6
        ca(1)=ca(2)
	ca(2)=cc
        if(ca(2).lt.ca(1).and.irepet.gt.1.and.ty.lt.1) ca(2)=ca(1)
        if(tt.ne.t(2)) t(1)=t(2)           
	t(2)=tt
c
        itbcnt=0
c
        if(irepet.eq.nrepet) then
           write(6,*) maxstp, nextra
           newstp=maxstp*nextra
           maxstp=newstp
           write(6,*) maxstp
           write(6,*) ' tbcnt= ',itbcnt
           itbcnt=0
        endif
c
c
c************************     time loop     ****************
c
	Do 16  nstp=1, maxstp
c
	   prt= kount/pintvl
	   prt=prt-int(prt)
	   If (prt.eq.0.d0) then
c	      write(6,44) x,
 	      write(36,44) x,
     c           y(1),y(n+1),y(2*n+1),y(3*n+1),y(4*n+1),y(5*n+1),dcal4,
     c           y(7*n+3),y(7*n+4),y(7*n+5),y(7*n+6),y(7*n+1),
     c           y(7*n+7),y(7*n+8),y(7*n+9),y(7*n+10),y(7*n+2),
     c           y(7),y(12),
     c	         iTcamf,iTcamb,iTcamtp,iTcama,iTcamcp,iTcan,iTcanb,
     c           y(7*n+12),y(7*n+13),y(7*n+14),zk(5)
  44	      format(f9.2, 19f8.3, 5i5,2i5,4e12.4)
c
c    c           y(4),y(n+4),y(2*n+4),y(3*n+4),y(4*n+4),y(5*n+4),
c    c	         y(6*n+1), bb41,
c    c	         y(6*n+2), bb42,y(7*n+1),y(7*n+2),
c
c              ixxxx=ixxxx+1
c              if(ixxxx.eq.10) then
c                  ixxxx=0
c             write(55,*) t
c             do 56 i=1,itholo
c             write(55,55) i,(icamk(i,j,1),j=1,isu)
c55           format(i4,10i3)
c56           continue
c             endif
c
c	      write(36,44) x,
c    c           y(1),y(n+1)+y(2*n+1)+y(3*n+1)+y(4*n+1)+y(5*n+1),
c    c               y(5*n+1),dcal4,
c    c           y(n+4)+y(2*n+4)+y(3*n+4)+y(4*n+4)+y(5*n+4),y(5*n+4),
c    c		     iTcamf,iTcamb,iTcamtp,iTcama,iTcamcp,iTcan,iTcanb
c 44	      format(f8.2, 6f8.3, 5i5,2i5)
                 do 48 ihis=1,4999     
                    iThis(ihis)=iThis(ihis+1)
 48              continue
                 iThis(5000)=iTcanb
                 iThis1=iThis(5000-iago)
	   endif
c
	   x=x+H
c
           if(nstp.gt.50000.and.
     c        y(1).lt.0.021.and.y(5*n+1).lt.(0.9*cnvrt)) then
                 continue
           else
	      call derivs(nvar,x,y,dydx,zk,amol,iTcan,iTcamf,iTcama,
     c  	 ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c	         dcal4,
     c           ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c		 iend)
c
     	      if (iend.eq.1) goto  321
c 
c***           need to add iend in  all (call derivs)  to continue camk etc rxn
c***** 
	      call RK4(Y,dydx,zk,Nvar,x,H,derivs,
     c              ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c		    dcal4,amol,iTcan,iTcamf,iTcama,
     c           ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c		    iend)
           endif
c
	   bb41=bt1-y(6*n+1)
	   bb42=bt2-y(6*n+2)
c
          do 500 iff=1,iTholo
             do 501 ifff=1,isu+1
                icamk(iff,ifff,1)=icamk(iff,ifff,2)
 501         continue
 500      continue
c
         do 9876 iazx=1,itholo
             ifreec=0
             ibounc=0
             itrapc=0
             iautoc=0
             icappc=0
            do 9875 iazc=2,isu+1
               if(icamk(iazx,iazc,2).eq.0) ifreec=ifreec+1
               if(icamk(iazx,iazc,2).eq.1) ibounc=ibounc+1
               if(icamk(iazx,iazc,2).eq.2) itrapc=itrapc+1
               if(icamk(iazx,iazc,2).eq.3) iautoc=iautoc+1
               if(icamk(iazx,iazc,2).eq.4) icappc=icappc+1
 9875      continue
            if(ifreec.ne.icamf(iazx)) write(6,*) 
     c          x,' free ',iazx,ifreec,icamf(iazx)
            if(ibounc.ne.icamb(iazx)) write(6,*) 
     c          x,' boun ',iazx,ibounc,icamb(iazx)
            if(itrapc.ne.icamtp(iazx)) write(6,*) 
     c          x,' trap ',iazx,itrapc,icamtp(iazx)
            if(iautoc.ne.icama(iazx)) write(6,*) 
     c          x,' auto ',iazx,iautoc,icama(iazx)
            if(icappc.ne.icamcp(iazx)) write(6,*) 
     c          x,' capp ',iazx,icappc,icamcp(iazx)
 9876     continue
              
c
          newba=0
          nb=iTcamb
c     do 8763 jklm=1,itholo
c         if(icamb(jklm).lt.0) then
c            write(6,*) ' 8763   pre   ',x,jklm,icamb(jklm),iTcamb
c         endif
c8763 continue
          if(y(7*n+1).gt.amol.or.(y(7*n+1).lt.amol.and.nb.gt.0)) then
             newb=y(7*n+1)/amol
             if(newb.gt.nb) then
                newba=newb-nb
                do 600 kb=1,newba
                   call RND(rr)
                   nsite=int(rr*iTcamf)+ 1
                   ksite=0
                   ihh=1
                   iss=1
                   do 111 while (ksite.lt.nsite)
                      iss=iss+1
                      if(iss.gt.isu+1) then
                         iss=2
                         ihh=ihh+1
                      endif
                      if(ihh.gt.itholo) write(6,*) ' too big ih',
     c                    ihh,iss,x,pr(1),pr(11),iTcamf,nsite,ksite,
     c                    y(7*n+1),amol,newb,nb,newba
                      if(icamk(ihh,iss,2).eq.0) ksite=ksite+1
 111               end do
c                  write(6,*) ' add ',ksite,nsite,iTcamf,rr
                   icamk(ihh,iss,2)=1
                   if(iss.eq.isu+1) icamk(ihh,1,2)=1
c                       write(6,*) ' 0->b ',x,pr(1),rr
                   iTcamf=iTcamf-1
                   iTcamb=iTcamb+1
                   icamb(ihh)=icamb(ihh)+1
                   icamf(ihh)=icamf(ihh)-1
 600            continue
             endif
      do 8762 jklm=1,itholo
          if(icamb(jklm).lt.0) then
             write(6,*) ' 8762   pre   ',x,jklm,icamb(jklm),iTcamb
          endif
 8762 continue
             if(newb.lt.nb) then
                newd=nb-newb
                newba=-newd
                do 700 ku=1,newd
                   call RND(rr)
                   nsite=int(rr*iTcamb)+ 1
                   ksite=0
                   ihh=1
                   iss=1
                   do 112 while (ksite.lt.nsite)
                      iss=iss+1
                      if(iss.gt.isu+1) then
                         iss=2
                         ihh=ihh+1
                      endif
                      if(ihh.gt.itholo) write(6,*) ' too big ih',
     c                    ihh,iss,x,pr(1),pr(11),iTcamf,nsite,ksite
                      if(icamk(ihh,iss,2).eq.1) ksite=ksite+1
c                     if(x.gt.243.422.and.ihh.eq.29) then
c                        write(6,*) nsite,ksite,ihh,iss
c                        write(6,*) newd,nb,newb,icamb(ihh) 
c                        write(6,*) icamk(ihh,iss,1),icamk(ihh,iss,2)
c                     endif
 112               end do
c                  write(6,*) ' sub ',ksite,nsite,iTcamb,rr
                   icamk(ihh,iss,2)=0
                   if(iss.eq.isu+1) icamk(ihh,1,2)=0
c        if(x.gt.243.422) then 
c             write(6,*) ' b->0 ',x,pr(1),rr,
c    c                    icamk(ihh,iss,1),
c    c                    icamk(ihh,iss,2),ihh,iss
c           write(6,*) nsite,ksite,iTcamb
c           endif
                   iTcamf=iTcamf+1
                   iTcamb=iTcamb-1
                   icamb(ihh)=icamb(ihh)-1
                   icamf(ihh)=icamf(ihh)+1
                   if(icamb(ih).lt.0) then
                      write(6,*) x,ihh,icamb(ihh),iss,iTcamb,
     c                    icamk(ihh,iss,1),
     c                    icamk(ihh,iss,2),ihh,iss
                   endif
 700            continue
             endif
          endif
c                                      auton to trapped
          tnew=-y(7*n+11)
          if(tnew.ge.amol.and.iTcama.gt.0) then
             newt=tnew/amol
             if(newt.gt.iTcama) newt=iTcama
             do 800 kb=1,newt
                call RND(rr)
                nsite=int(rr*iTcama)+ 1
                ksite=0
                ihh=1
                iss=1
                do 811 while (ksite.lt.nsite)
                   iss=iss+1
                   if(iss.gt.isu+1) then
                      iss=2
                      ihh=ihh+1
                   endif
                   if(ihh.gt.itholo) write(6,*) ' too big ih',
     c                 ihh,iss,x,pr(1),pr(9),iTcama,nsite,ksite
                   if(icamk(ihh,iss,2).eq.3) ksite=ksite+1
 811            end do
c               write(6,*) ' add ',ksite,nsite,iTcamf,rr
                icamk(ihh,iss,2)=2
                if(iss.eq.isu+1) icamk(ihh,1,2)=2
c                    write(6,*) ' a->t ',x,pr(9),rr
                iTcama=iTcama-1
                iTcamtp=iTcamtp+1
                icamtp(ihh)=icamtp(ihh)+1
                icama(ihh)=icama(ihh)-1
                y(7*n+11)=y(7*n+11)+amol
                y(5*n+1)=y(5*n+1)-amol
 800         continue
          endif
c
          iTcanb=y(7*n+2)/amol
          iTcan=(cnt/amol)-iTcanb

      do 8764 jklm=1,itholo
          if(icamb(jklm).lt.0) then
             write(6,*) ' camkin pre   ',x,jklm,icamb(jklm),iTcamb
          endif
 8764 continue
c
	   Kount=kount+1
c
           iccall=0
	   bincall= dfloat(kount)/bintvl
	   bincall=bincall-int(bincall+.01*h)
	   camkcall= dfloat(kount)/camkintvl
	   camkcall=camkcall-int(camkcall+.01*h)
c	   If (bincall.lt.1.d-9) then
	      If (camkcall.lt.1.d-9) then
                 iccall=1
                 camkph=dfloat(iTcama+iTcamtp+iTcamc)*amol    
c                write(6,*) ' time= ',x
                 zk(4)=zk4*1d-4
c                zk(4)=zk4*1d-3/(.00228*((y(1)*1d3)**1.6919)+9.88)
c                zk(4)=zk4*dexp(-(y(1)-.01d0)/.43d0)
		 zk(5)=ak2pp*y(7*n+14)*camkph/(ppkm+camkph)
		 zk(7)=ak2pp*y(7*n+14)*camkph/(ppkm+camkph)
		 zk(8)=ak2pp*y(7*n+14)*camkph/(ppkm+camkph)
		 zk(10)=ak2pp*y(7*n+14)*camkph/(ppkm+camkph)
cc               zk(5)=zk5*(1.d0+dfloat(iThis1)/(cnt/amol))**ampf
cc   c              -zk5
cc               zk(7)=zk7*(1.d0+dfloat(iThis1)/(cnt/amol))**ampf
cc   c              -zk7
cc               zk(8)=zk8*(1.d0+dfloat(iThis1)/(cnt/amol))**ampf
cc   c              -zk8
cc               zk(10)=zk10*(1.d0+dfloat(iThis1)/(cnt/amol))**ampf
cc   c              -zk10
c                zk(5)=zk5*dfloat(iThis1)/(cnt/amol)*ampf+zk5
c                zk(7)=zk7*dfloat(iThis1)/(cnt/amol)*ampf+zk7
c                zk(8)=zk8*dfloat(iThis1)/(cnt/amol)*ampf+zk8
c                zk(10)=zk10*dfloat(iThis1)/(cnt/amol)*ampf
c    c              +zk10
              if(iThis1.gt.int(cnt/amol)) then
                 write(6,*) x,iThis1, iTcanb,iTcan,ampf
                 write(6,*) ' compiler problem check results '
              endif
	      hh=h*camkintvl
              hhh=h*bintvl
c234567890123456789012345678901234567890123456789012345678901234567890
              call camkin(rnd,y(5*n+1),y(7*n+1),hh,hhh,zk,pr,itholo,isu,
     c           autom,
     c   	 icamf,icamb,icamtp,icama,icamcp,iphos,icamk,
     c           iTcan,iTcanb,iTcamf,iTcamb,iTcamtp,iTcama,iTcamcp,
     c		 avePh,aveTp,iactiho,ipo4,itp,
     c           ibd,iau,icp,ifr,iccall,itbcnt,
     c           cnvrt,x)
              endif
	      dcal4=cnvrt*dfloat(iTcamb+iTcamtp+iTcanb)
c	  endif 


 16	continue
c
c      *******        end of time loop   ****************************
 321	close(unit=iu)  
c
           write(6,*) ' tbcnt= ',itbcnt
 17     continue
c     *******       end of repeat loop   ***********************
c 
	close(unit=36)
c       write(6,*) ' enter output file name '
        read(ir,102) file3
        open(unit=62, file=file3,status='unknown')
        do 62 i=1,itholo
           write(62,72) (icamk(i,iss,2),iss=2,isu+1)
 72        format(6i4)
 62     continue
        close(unit=62)
           write(6,*) ' tbcnt= ',itbcnt
        write(41,52) file2, (j,j=0,isu)
 52     format(a12,7i4)
        write(41,53) ifr
        write(41,53) ibd
        write(41,53) itp
        write(41,53) iau
        write(41,53) icp
 53     format('            ',7i4)
	return
	END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE RK4(Y,DYDX,zk,Nvar,X,H,DERIVS,
     c	 ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c		    dcal4,amol,iTcan,iTcamf,iTcama,
     c       ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c		    iend)
	implicit double precision(a-h,o-z)
c     PARAMETER (NMAX=120)
      DIMENSION Y(nvar),DYDX(nvar),
     c	 YT(nvar),DYT(nvar),DYM(nvar),
     c     ap(4),am(4),bp(4),bm(4),cp(4),cm(4),
     c     yp(5),ym(5),zp(5),zm(5),
     c	 p(n),a(n,n),b(n,n),zk(12)
      HH=H*0.5d0
      H6=H/6.d0
      XH=X+HH
      DO 11 I=1,Nvar
	YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(nvar,XH,YT,DYT,zk,amol,iTcan,iTcamf,iTcama,
     c	 ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c		    dcal4,
     c      ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c		    iend)
      DO 12 I=1,Nvar
	YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(nvar,XH,YT,DYM,zk,amol,iTcan,iTcamf,iTcama,
     c	 ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c		    dcal4,
     c      ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c		    iend)
      DO 13 I=1,Nvar
	YT(I)=Y(I)+H*DYM(I)
	DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(nvar,X+H,YT,DYT,zk,amol,iTcan,iTcamf,iTcama,
     c	 ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c		    dcal4,
     c      ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c		    iend)
      DO 14 I=1,Nvar
	Y(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.d0*DYM(I))
14    CONTINUE
c     if(x.gt.117.7) write(6,*) ' in RKQ$ '
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE DERIVS(nvar,X,Y,DYDX,zk,amol,iTcan,iTcamf,iTcama,
     c	 ap,am,bp,bm,cp,cm,yp,ym,zp,zm,n,p,a,b,lsh,lsn,isd,
     c   dcal4,
     c       ak2i1,ai1t,pka,akmi1,ak2i1p,akmi1p,akp1pp,akm1pp,pp1t,
     c	 iend)
	implicit double precision(a-h,o-z)
c
	COMMON /PATH/KMAX,KOUNT,DXSAV,XP(3000),
     c	 pi,ct1,ct2,bt1,bt2,cf(4),cb(4),bb,bf,y0,
     c   ca(2),t(2),col(14),fff,zfvol,tend,
     c   ncol,icc,ica,iu,irepet
c
	DIMENSION Y(nvar),DYDX(nvar),p(n),a(n,n),b(n,n),
     c     ap(4),am(4),bp(4),bm(4),cp(4),cm(4),
     c     yp(5),ym(5),zp(5),zm(5),
     c    zk(12)
c
	pi=4.d0*datan(1.d0)
 1	if(x.ge.t(2)) then
	   read(iu,*,end=99) ty,(col(i),i=1,ncol)
            tt=ty+dfloat(irepet-1)*tend
	   cc=col(icc)*fff
	   t(1)=t(2)
	   t(2)=tt
c
c    1e-12C/sec*1e6(umol/mol)/(zF*vol(pi*.275**2*.05e-15liters)*1e3ms/sec)
c
c          cc=cc/(2.0*9.648456e4*pi*.275*.275*.05)*1e6
	   cc=cc/zfvol*1d6
	   ca(1)=ca(2)
	   ca(2)=cc
        if(ca(2).lt.ca(1).and.irepet.gt.1.and.ty.lt.1) ca(2)=ca(1)
	endif
	if(x.gt.t(2)) goto 1
c
	tx=x-t(1)
	cca=ca(1)+tx/(t(2)-t(1))*(ca(2)-ca(1))
c       if(x.lt.0.02) write(53,*) ' cca= ',cca,x,ca,col(icc)
c       if(x.gt.99.98.and.x.lt.100.02) write(53,*) x,cca,ca,col(icc)
c
c	write(6,*) ' y= '
c	write(6,64) a
c64	format(11f7.2)
c	write(6,65) y
c65	format(11f7.3)
c	write(6,66) p
c66	format(11f7.5)
c
        freeK=dfloat(iTcamf)*amol-y(7*n+3)-y(7*n+4)-y(7*n+5)-y(7*n+6)
        freeN=dfloat(iTcan)*amol-y(7*n+7)-y(7*n+8)-y(7*n+9)-y(7*n+10)
        autonm=dfloat(iTcama)*amol
c
	DYDX(1)=-A(1,2)*y(1)+A(1,2)*y(2)-
     c	   y(1)*(cf(1)*y(n+1)+cf(2)*y(2*n+1)+
     c       cf(3)*y(3*n+1)+cf(4)*y(4*n+1))+
     c     cb(1)*y(2*n+1)+cb(2)*y(3*n+1)+cb(3)*y(4*n+1)+
     c       cb(4)*y(5*n+1)-
     c       bf*y(1)*y(6*n+1)+bb*(bt1-y(6*n+1))-
     c	   p(1)*(y(1)-y0)
     c     +am(1)*y(7*n+4)+am(2)*y(7*n+5)+am(3)*y(7*n+6)
     c     +am(4)*y(7*n+1)
     c     -y(1)*(ap(1)*y(7*n+3)+ap(2)*y(7*n+4)+ap(3)*y(7*n+5)
     c       +ap(4)*y(7*n+6))
     c     +cm(1)*y(7*n+8)+cm(2)*y(7*n+9)+cm(3)*y(7*n+10)
     c     +cm(4)*y(7*n+2)
     c     -y(1)*(cp(1)*y(7*n+7)+cp(2)*y(7*n+8)+cp(3)*y(7*n+9)
     c       +cp(4)*y(7*n+10))

c***************   in the 1st cmpt above  need to add  cam kinase -cal4 rxn 


	DYDX(2)=-A(2,1)*y(2)+A(2,1)*y(1)-A(2,3)*y(2)+A(2,3)*y(3)-
     c	   y(2)*(cf(1)*y(n+2)+cf(2)*y(2*n+2)+
     c       cf(3)*y(3*n+2)+cf(4)*y(4*n+2))+
     c     cb(1)*y(2*n+2)+cb(2)*y(3*n+2)+cb(3)*y(4*n+2)+
     c       cb(4)*y(5*n+2)-
     c       bf*y(2)*y(6*n+2)+bb*(bt2-y(6*n+2))-
     c	   p(2)*(y(2)-y0)
c
	do 100 i=3,lsn-1
	   DYDX(i)=-a(i,i-1)*y(i)+a(i,i-1)*y(i-1)-
     c	      a(i,i+1)*y(i)+a(i,i+1)*y(i+1)-
     c          y(i)*(cf(1)*y(n+i)+cf(2)*y(2*n+i)+
     c          cf(3)*y(3*n+i)+cf(4)*y(4*n+i))+
     c        cb(1)*y(2*n+i)+cb(2)*y(3*n+i)+cb(3)*y(4*n+i)+
     c          cb(4)*y(5*n+i)-
     c       bf*y(i)*y(6*n+i)+bb*(bt2-y(6*n+i))-
     c	         p(i)*(y(i)-y0)
 100	continue
c
	dydx(lsn)=-a(lsn,lsn-1)*y(lsn)+a(lsn,lsn-1)*y(lsn-1)-
     c	   a(lsn,isd)*y(lsn)+a(lsn,isd)*y(isd)-
     c	   a(lsn,isd+1)*y(lsn)+a(lsn,isd+1)*y(isd+1)-
     c        y(lsn)*(cf(1)*y(n+lsn)+cf(2)*y(2*n+lsn)+
     c          cf(3)*y(3*n+lsn)+cf(4)*y(4*n+lsn))+
     c        cb(1)*y(2*n+lsn)+cb(2)*y(3*n+lsn)+cb(3)*y(4*n+lsn)+
     c          cb(4)*y(5*n+lsn)-
     c       bf*y(lsn)*y(6*n+lsn)+bb*(bt2-y(6*n+lsn))-
     c	   p(lsn)*(y(lsn)-y0)
	DYDX(lsn+1)=-A(lsn+1,lsn+2)*y(lsn+1)+
     c	   A(lsn+1,lsn+2)*y(lsn+2)-
     c        y(lsn+1)*(cf(1)*y(n+lsn+1)+cf(2)*y(2*n+lsn+1)+
     c          cf(3)*y(3*n+lsn+1)+cf(4)*y(4*n+lsn+1))+
     c        cb(1)*y(2*n+lsn+1)+cb(2)*y(3*n+lsn+1)+cb(3)*y(4*n+lsn+1)+
     c          cb(4)*y(5*n+lsn+1)-
     c       bf*y(lsn+1)*y(6*n+lsn+1)+bb*(bt2-y(6*n+lsn+1))-
     c	   p(lsn+1)*(y(lsn+1)-y0)
c
	Do 200 i=lsn+2,n-1
	   DYDX(i)=-a(i,i-1)*y(i)+a(i,i-1)*y(i-1)-
     c	      a(i,i+1)*y(i)+a(i,i+1)*y(i+1)-
     c          y(i)*(cf(1)*y(n+i)+cf(2)*y(2*n+i)+
     c          cf(3)*y(3*n+i)+cf(4)*y(4*n+i))+
     c        cb(1)*y(2*n+i)+cb(2)*y(3*n+i)+cb(3)*y(4*n+i)+
     c          cb(4)*y(5*n+i)-
     c       bf*y(i)*y(6*n+i)+bb*(bt2-y(6*n+i))-
     c	      p(i)*(y(i)-y0)
 200	continue
c
	DYDX(n)=-A(n,n-1)*y(n)+A(n,n-1)*y(n-1)-
     c          y(n)*(cf(1)*y(n+n)+cf(2)*y(2*n+n)+
     c          cf(3)*y(3*n+n)+cf(4)*y(4*n+n))+
     c        cb(1)*y(2*n+n)+cb(2)*y(3*n+n)+cb(3)*y(4*n+n)+
     c          cb(4)*y(5*n+n)-
     c       bf*y(n)*y(6*n+n)+bb*(bt2-y(6*n+n))-
     c	   p(n)*(y(n)-y0)
c
	dydx(isd)=dydx(isd)-
     c	   a(isd,lsn)*y(isd)+a(isd,lsn)*y(lsn)
	dydx(isd+1)=dydx(isd+1)-
     c	   a(isd+1,lsn)*y(isd+1)+a(isd+1,lsn)*y(lsn)
c
        do 240 i=1,n
           if(i.eq.1) btt=bt1
           if(i.ne.1) btt=bt2
           dydx(n+i)=-cf(1)*y(i)*y(n+i)+cb(1)*y(2*n+i)
           dydx(2*n+i)=cf(1)*y(i)*y(n+i)-y(2*n+i)*(cb(1)+cf(2)*y(i))
     c        +cb(2)*y(3*n+i)
           dydx(3*n+i)=cf(2)*y(i)*y(2*n+i)-y(3*n+i)*(cb(2)+cf(3)*y(i))
     c        +cb(3)*y(4*n+i)
           dydx(4*n+i)=cf(3)*y(i)*y(3*n+i)-y(4*n+i)*(cb(3)+cf(4)*y(i))
     c        +cb(4)*y(5*n+i) 
           dydx(5*n+i)=cf(4)*y(i)*y(4*n+i)-y(5*n+i)*cb(4)
           dydx(6*n+i)=-bf*y(i)*y(6*n+i)+bb*(btt-y(6*n+i))
 240    continue
c
        dydx(7*n+1)=zk(1)*y(5*n+1)*freeK -
     c       zk(2)*y(7*n+1)
     c       +ap(4)*y(7*n+6)*y(1)
     c       -y(7*n+1)*am(4)
        dydx(7*n+2)=zk(11)*y(5*n+1)*freeN -
     c       zk(12)*y(7*n+2)
     c       +cp(4)*y(7*n+10)*y(1)
     c       -y(7*n+2)*cm(4)
c
        dydx(7*n+3)=yp(1)*freeK*y(n+1) + am(1)*y(7*n+4)
     c       -y(7*n+3)*(ym(1)+ap(1)*y(1))
        dydx(7*n+4)=yp(2)*freeK*y(2*n+1) + am(2)*y(7*n+5)
     c        + ap(1)*y(7*n+3)*y(1)
     c       -y(7*n+4)*(ym(2)+ap(2)*y(1)+am(1))
        dydx(7*n+5)=yp(3)*freeK*y(3*n+1) + am(3)*y(7*n+6)
     c        + ap(2)*y(7*n+4)*y(1)
     c       -y(7*n+5)*(ym(3)+ap(3)*y(1)+am(2))
        dydx(7*n+6)=yp(4)*freeK*y(4*n+1) + am(4)*y(7*n+1)
     c        + ap(3)*y(7*n+5)*y(1)
     c       -y(7*n+6)*(ym(4)+ap(4)*y(1)+am(3))
c
        dydx(7*n+7)=zp(1)*freeN*y(n+1) + cm(1)*y(7*n+8)
     c       -y(7*n+7)*(zm(1)+cp(1)*y(1))
        dydx(7*n+8)=zp(2)*freeN*y(2*n+1) + cm(2)*y(7*n+9)
     c        + cp(1)*y(7*n+7)*y(1)
     c       -y(7*n+8)*(zm(2)+cp(2)*y(1)+cm(1))
        dydx(7*n+9)=zp(3)*freeN*y(3*n+1) + cm(3)*y(7*n+10)
     c        + cp(2)*y(7*n+8)*y(1)
     c       -y(7*n+9)*(zm(3)+cp(3)*y(1)+cm(2))
        dydx(7*n+10)=zp(4)*freeN*y(4*n+1) + cm(4)*y(7*n+2)
     c        + cp(3)*y(7*n+9)*y(1)
     c       -y(7*n+10)*(zm(4)+cp(4)*y(1)+cm(3))
c
        dydx(7*n+11)=-zk(9)*autonm*y(5*n+1)
c
        dydx(5*n+1)=dydx(5*n+1)+zk(2)*y(7*n+1)+zk(12)*y(7*n+2) -
     c       y(5*n+1)*(zk(1)*freeK + zk(11)*freeN +zk(9)*autonm)
        dydx(4*n+1)=dydx(4*n+1)
     c       +ym(4)*y(7*n+6)+zm(4)*y(7*n+10)
     c       -y(4*n+1)*(yp(4)*freeK+zp(4)*freeN)
        dydx(3*n+1)=dydx(3*n+1)
     c       +ym(3)*y(7*n+5)+zm(3)*y(7*n+9)
     c       -y(3*n+1)*(yp(3)*freeK+zp(3)*freeN)
        dydx(2*n+1)=dydx(2*n+1)
     c       +ym(2)*y(7*n+4)+zm(2)*y(7*n+8)
     c       -y(2*n+1)*(yp(2)*freeK+zp(2)*freeN)
        dydx(n+1)=dydx(n+1)
     c       +ym(1)*y(7*n+3)+zm(1)*y(7*n+7)
     c       -y(n+1)*(yp(1)*freeK+zp(1)*freeN)
c
        do 250 j=1,5
           dydx(j*n+1)=dydx(j*n+1)-b(1,2)*y(j*n+1)+b(1,2)*y(j*n+2)
           dydx(j*n+2)=dydx(j*n+2)-b(2,1)*y(j*n+2)+b(2,1)*y(j*n+1)-
     c       b(2,3)*y(j*n+2)+b(2,3)*y(j*n+3)
           do 255 k=3,lsn-1
              dydx(j*n+k)=dydx(j*n+k)
     c           -b(k,k-1)*y(j*n+k)+b(k,k-1)*y(j*n+k-1)
     c           -b(k,k+1)*y(j*n+k)+b(k,k+1)*y(j*n+k+1)
 255       continue
           dydx(j*n+lsn)=dydx(j*n+lsn)
     c        -b(lsn,lsn-1)*y(j*n+lsn)+b(lsn,lsn-1)*y(j*n+lsn-1)
     c        -b(lsn,isd)*y(j*n+lsn)+b(lsn,isd)*y(j*n+isd)
     c        -b(lsn,isd+1)*y(j*n+lsn)+b(lsn,isd+1)*y(j*n+isd+1)
           dydx(j*n+lsn+1)=dydx(j*n+lsn+1)
     c        -b(lsn+1,lsn+2)*y(j*n+lsn+1)+b(lsn+1,lsn+2)*y(j*n+lsn+2)
           do 258 k=lsn+2,n-1
              dydx(j*n+k)=dydx(j*n+k)
     c           -b(k,k-1)*y(j*n+k)+b(k,k-1)*y(j*n+k-1)
     c           -b(k,k+1)*y(j*n+k)+b(k,k+1)*y(j*n+k+1)
 258       continue
           dydx(j*n+n)=dydx(j*n+n)
     c        -b(n,n-1)*y(j*n+n)+b(n,n-1)*y(j*n+n-1)
c
           dydx(j*n+isd)=dydx(j*n+isd)-
     c        b(isd,lsn)*y(j*n+isd)+b(isd,lsn)*y(j*n+lsn)
           dydx(j*n+isd+1)=dydx(j*n+isd+1)-
     c        b(isd+1,lsn)*y(j*n+isd+1)+b(isd+1,lsn)*y(j*n+lsn)
c
 250    continue
c
	dydx(ica)=dydx(ica)+cca
c
c  12=I1
c  13=I1P
c  14=PP1
c
c        if(x.lt.0.001) write(6,*) pka,ai1t,pp1t,akp1pp,akm1pp
        dydx(7*n+12)=-ak2i1*pka*y(7*n+12)/(akmi1+y(7*n+12)) +
     c      ak2i1p*y(7*n+2)*y(7*n+13)/(akmi1p+y(7*n+13))
        dydx(7*n+13)=-ak2i1p*y(7*n+2)*y(7*n+13)/(akmi1p+y(7*n+13)) +
     c      ak2i1*pka*y(7*n+12)/(akmi1+y(7*n+12)) -
     c      akp1pp*y(7*n+14)*y(7*n+13)+
     c      akm1pp*(ai1t-y(7*n+12)-y(7*n+13))
        dydx(7*n+14)=-akp1pp*y(7*n+14)*y(7*n+13) +
     c      akm1pp*(pp1t-y(7*n+14))
c
c
c	write(6,44) dydx(7),dydx(12),dydx(7*n+1)
 	if(x.lt.0.001) 
     c     write(6,44) dydx(7*n+12),dydx(7*n+13),dydx(7*n+14)
 44	format(3e12.4)
 	RETURN
c
 99	iend=1
	RETURN
	END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine geta(n,al,d,a,tn,r,ar,v,ir,id,lsh,lsn,isd,
     c		dif,cnvrt)
	implicit double precision(a-h,o-z)
       dimension al(n),d(n),a(n,n),tn(n),r(n),ar(n),v(n)
       character*40 dimfil
       pi=4.d0*datan(1.d0)
c
	An=6.023e23
c
c	write(6,*) ' enter dimension file '
       read(ir,10) dimfil
 10    format(a40)
       open(file=dimfil,unit=id,status='old')
       do 200 i=1,n
	  read(id,*) ic,al(i),d(i)
 200   continue
c      write(6,*) ' enter last spine head compartment '
       read(id,*) lsh
c      write(6,*) ' enter last spine neck compartment '
       read(id,*) lsn
c      write(6,*) ' enter first compartment neck is connected to '
       read(id,*) isd
c
       close(unit=id)
c
       do 300 i=1,n
	  r(i)=d(i)*.5d0
	  tn(i)=al(i)*.5d0
	  ar(i)=pi*r(i)*r(i)
	  v(i)=ar(i)*al(i)
 300   continue
c
	cnvrt=1/(An*v(1))*(1e21)
       write(6,*) ' cnvrt (amol) = ',cnvrt
c
       dif=.6d0
c      write(6,*) ' enter calcium diffusion const '
       read(ir,*) dif
	  do 310 i=1,lsn-1
	     c=dif/v(i)
	     th=tn(i)+tn(i+1)
	     th2=th*th
	     a(i,i+1)=c*(ar(i)*tn(i)+ar(i+1)*tn(i+1))/th2
 310	  continue
	  do 312 i=lsn+1,n-1
	     c=dif/v(i)
	     th=tn(i)+tn(i+1)
	     th2=th*th
	     a(i,i+1)=c*(ar(i)*tn(i)+ar(i+1)*tn(i+1))/th2
 312	  continue
	  do 320 i=2,lsn
	     c=dif/v(i)
	     th=tn(i)+tn(i-1)
	     th2=th*th
	     a(i,i-1)=c*(ar(i)*tn(i)+ar(i-1)*tn(i-1))/th2
 320	  continue
	  do 322 i=lsn+2,n
	     c=dif/v(i)
	     th=tn(i)+tn(i-1)
	     th2=th*th
	     a(i,i-1)=c*(ar(i)*tn(i)+ar(i-1)*tn(i-1))/th2
 322	  continue
c
	  cn=dif/v(lsn)
	  cd1=dif/v(isd)
	  cd2=dif/v(isd+1)
	  thn1=tn(lsn)+tn(isd)
	  thn2=tn(lsn)+tn(isd+1)
	  th12=thn1*thn1
	  th22=thn2*thn2
	  a(lsn,isd)=cn*(ar(lsn)*tn(lsn)+ar(isd)*tn(isd))/th12
	  a(lsn,isd+1)=cn*(ar(lsn)*tn(lsn)+ar(isd+1)*tn(isd+1))/th22
	  a(isd,lsn)=cd1*(ar(lsn)*tn(lsn)+ar(isd)*tn(isd))/th12
	  a(isd+1,lsn)=cd2*(ar(lsn)*tn(lsn)+ar(isd+1)*tn(isd+1))/th22
c

 
ccccccccccccccccccccccccccccccccccccccccccccccc
           c=dif/v(lsh)
           th=tn(lsh)+tn(lsh+1)
           a(lsh,lsh+1)=c*ar(lsh+1)/th
           c=dif/v(lsh+1)
           th=tn(lsh+1)+tn(lsh)
           a(lsh+1,lsh)=c*ar(lsh+1)/th
           cn=dif/v(lsn)
           cd1=dif/v(isd)
           cd2=dif/v(isd+1)
           thn1=tn(lsn)+ tn(isd)
           thn2=tn(lsn)+ tn(isd+1)
           th12=thn1*thn1
           th22=thn2*thn2
           a(lsn,isd)=cn*ar(lsn)*.5/thn1
           a(lsn,isd+1)=cn*ar(lsn)*.5/thn2
           a(isd,lsn)=cd1*ar(lsn)*.5/thn1
           a(isd+1,lsn)=cd2*ar(lsn)*.5/thn2
ccccccccccccccccccccccccccccccccccccccccccccccc
c
c	 do 1000 i=1,n
c	    write(6,21) (a(i,j),j=1,n)
c 21	    format(12f8.2)
c1000	 continue
c
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine camkin(rnd,cal4,ckb,hh,hhh,zk,pr,itholo,isu,autom,
     c		icamf,icamb,icamtp,icama,icamcp,iphos,icamk,
     c          iTcan,iTcanb,iTcamf,iTcamb,iTcamtp,iTcama,iTcamcp,
     c		avePh, aveTp, iactiho,ipo4,itp,
     c          ibd,iau,icp,ifr,iccall,itbcnt,
     c		cnvrt,x)
	implicit double precision (a-h,o-z)
c
	DIMENSION zk(12), pr(12), h1(itholo), h2(itholo), h3(itholo),
     c		icamk(itholo,isu+1,2), icamb(itholo),icamf(itholo),
     c		ipo4(3*isu+1), itp(isu+1), iphos(itholo),
     c		icamtp(itholo),icama(itholo),icamcp(itholo),
     c          iau(isu+1),icp(isu+1),ibd(isu+1),ifr(isu+1)
c
c     do 8765 jklm=1,itholo
c         if(icamb(jklm).lt.0) then
c            write(6,*) ' camkin start ',x,jklm,icamb(jklm),iTcamb
c         endif
c8765 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  rnd   	external subroutine for random number generator
c  cal4		y(5*n+1)  [Cal-Ca4]
c  hh		time interval for cal-ca4 reactions
c  zk(12)	rate constants
c  pr(12)	rate constants * time interval  (probabilities)
c  itholo	# of holoenzymes. Set to 100
c  isu		# subunits per holoenzyme.  Set to 10
c  autom	Set to 0.4  reduction in probability
c  icamf(itholo)	#subunits free within a holoenzyme
c  icamb(itholo)	# subunits bound with cal-ca4 in a holoenzyme
c  icamtp(itholo)	# subunits bound and trapped in a holoenzyme
c  icama(itholo)	# subunits in autonomous state in a holoenzyme
c  icamcp(itholo)	# subunits in capped state in a holoenzyme
c  iphos(itholo)	Not used
c  icamk(itholo,isu+1,2)  Status of subunits in holoenzymes
c                         1=present status, 2=future status
c  iTcan	# free calcineurin molecules
c  iTcanb	# calcineurin molecules bound with cal-ca4
c  iTcamf	# caM subunits free total
c  iTcamb	# caM subunits bound with Cal-Ca4 total
c  iTcamtp	# caM subunits in trapped state total
c  iTcama	# caM subunits in autonomous state total
c  iTcamcp	# caM subunits in capped state total
c  avePh	not used      
c  aveTp	iTcamtp/itholo
c  iactiho	Not used
c  ipo4(3*isu+1)	Not used
c  itp(isu+1)	# holoenzymes with 0-10 subunits trapped
c  ibd(isu+1)	# holoenzymes with 0-10 subunits bound
c  iau(isu+1)	# holoenzymes with 0-10 subunits autonomous
c  icp(isu+1)	# holoenzymes with 0-10 subunits capped
c  ifr(isu+1)	# holoenzymes with 0-10 subunits free
c  cnvrt	conversion factor conc to number of molecules
c  x		time
c  h1(itholo)	Not used
c  h2(itholo)	Not used
c  h3(itholo)	Not used
ccccccccccccccccccccccccccccccccccccccc
	pi=4.0d0*DATAN(1.d0)
        amol=cnvrt
c
c            probabilities = rate constants * time interval
c
 	Do 20 irxn=1, 12
 	   pr(irxn)=1.d0-dexp(-zk(irxn)*hh)
c          if(pr(irxn).lt.1d-8.and.irxn.ne.5.and.irxn.ne.7.
c    c       and.irxn.ne.8.and.irxn.ne.10) 
c    c       write(6,*) ' hi ',x,irxn,pr(irxn)
 20 	continue
c       if(x.gt.47.7.and.x.lt.48) write(96,*) ' pr= ',pr
c
c            set present status from previous future status
c
	Do 69 ih=1, iTholo
	   Do 70 i=1,isu+1 
	      icamk(ih,i,1)=icamk(ih,i,2)
 70	   continue
 69	continue
c
c
c
c*************************  does a cal-ca4 bind?
c              cal-ca4 can bind to  autonom subunit
c
cc      idum=int(cal4/cnvrt)
c
cc      do 100 while (idum.gt.0.and.iTcama.gt.0)
c               if(x.gt.117) write(6,*) ' 100 loop ',x,idum
cc        pr(9)=1.d0-dexp(-zk(9)*hhh*cnvrt*iTcama)
cc        if(pr(9).gt.0.5) write(6,*) ' pr(9) too big ',pr(9),x
cc        call RND(rr)
c
c                                      bind to autom camK  
c                                      i=#free sites, n=bind site
cc        if(rr.lt.pr(9)) then
c            if(x.gt.117.7) write(6,*) ' bind a->t ',pr(9)
cc           nsite=int(rr/pr(9)*iTcama) + 1
cc           ksite=0
cc           ihh=1
cc           iss=1
cc           do 120 while (ksite.lt.nsite)
c               if(x.gt.117) write(6,*) ' 120 loop ',x,nsite,ksite
cc              iss=iss+1
cc              if(iss.gt.isu+1) then
cc                 iss=2
cc                 ihh=ihh+1
cc              endif
cc              if(ihh.gt.itholo) write(6,*) ' too biga ih',ihh,iss,x,
cc   c              pr(9),iTcama,nsite,ksite
cc              if(icamk(ihh,iss,1).eq.3) ksite=ksite+1
cc120         end do
c            write(6,*) ksite,nsite
cc           icamk(ihh,iss,2)=2
cc           if(iss.eq.isu+1) icamk(ihh,1,2)=2
cc           iTcama=iTcama-1
cc           iTcamtp=iTcamtp+1
cc           icama(ihh)=icama(ihh)-1
cc           icamtp(ihh)=icamtp(ihh)+1
cc           cal4=cal4-cnvrt
cc        endif
cc        idum=idum-1
cc100   continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          if interval is not for unbinding skip the next code
c
       if(iccall.ne.1) goto 5000
c       write(96,*) ' in camkin ',x
c
c********************************************
	Do 11 ih=1, iTholo
	   idbnd=0
           idtrap=0
           idau=0
           idcap=0
           idfree=0
	   icamk(ih,1,1)=icamk(ih,isu+1,1)
c          if(x.gt.50.and.x.lt.50.2) then
c             write(66,*) x
c             write(66,66) ((icamk(ix,iy,1),iy=1,isu+1),ix=1,itholo)
c66           format(11i5)
c          endif
c
           Do 21  i=2, isu+1
c
              if(icamb(ih).lt.0) then
                 write(6,*) x,i,ih,icamb(ih),idbnd
              endif
	      call RND(rr)
	      ici=icamk(ih,i,1)
	      icim1=icamk(ih,i-1,1)
              if(i.eq.isu+1) icip1=icamk(ih,2,1)
              if(i.lt.isu+1) icip1=icamk(ih,i+1,1)
              if(ici.ne.icamk(ih,i,2)) goto 21
c
	      if (ici.eq.1) then
                 if ((icim1.eq.1.or.icim1.eq.2.or.icip1.eq.1.
     c                or.icip1.eq.2).and.rr.lt.pr(3)) then
c                       if(x.gt.47.7.and.x.lt.48) write(6,*) ' b->t ',
c    c                       x,pr(3),rr
c                      write(6,*) ' testing a b ',x,rr,pr(3)
		    icamk(ih,i,2) = 2
		    idbnd = idbnd - 1
		    idtrap = idtrap + 1
                    ckb=ckb-amol
		 endif	
	         if ((icim1.eq.3.or.icim1.eq.4.or.icip1.eq.3.
     c                 or.icip1.eq.4).and.icamk(ih,i,2).eq.1.and.
     c                   (rr.lt.pr(3)*autom)) then
c                        write(6,*) ' testing a b2 ',x,rr,pr(3)
c                       if(x.gt.117.7) write(6,*) 
c    c                             ' b->t2 ',x,pr(2),pr(2)+pr(3),rr
	            icamk(ih,i,2) = 2
		    idbnd = idbnd - 1
		    idtrap = idtrap + 1
                    ckb=ckb-amol
	         endif
                 if (icim1.eq.0.and.icip1.eq.0.and.rr.lt.pr(3)*0.01d0
     c              .and.icamk(ih,i,2).eq.1) then
c                       if(x.gt.117.7) write(6,*) 
c                            write(6,*) ' testing a b3 ',x,rr,pr(3)
c    c                             ' b->tr3',x,pr(2),pr(2)+pr(3),rr,ih,i
	            icamk(ih,i,2) = 2
	            idbnd = idbnd - 1
	            idtrap = idtrap + 1
                    ckb=ckb-amol
	         endif
	      endif
c
	      if (ici.eq.2) then
		   if (rr.lt.pr(4)) then
c                       write(6,*) ' t->a ',x,pr(4),rr
		      icamk(ih,i,2) = 3
		      idtrap = idtrap - 1	
		      idau = idau + 1	
		      cal4=cal4+cnvrt
		   endif
                   if(rr.ge.pr(4).and.rr.lt.(pr(4)+pr(8))) then
                        if(x.gt.7.7) then 
                            write(6,*) ' t->b  ',x,pr(4),pr(4)+pr(8),rr
                            itbcnt=itbcnt+1
                        endif
                      icamk(ih,i,2) = 1
                      idtrap = idtrap - 1
		      idbnd = idbnd + 1
                      ckb=ckb+amol
                   endif
	      endif
c
	      if (ici.eq.3.) then
		  if (rr.gt.pr(5)) then
c                    pp=(icamb(ih)+icamtp(ih)+icama(ih)+icamcp(ih))/5
                     if (rr.lt.(pr(5)+pr(6)).and.
     c                      (icim1.ne.0.or.icip1.ne.0)) then
c                       write(6,*) ' a->c ',x,pr(5),pr(5)+pr(6),rr
		         icamk(ih,i,2) = 4
		         idau = idau - 1	
		         idcap = idcap + 1	
		     endif
		  else
                        write(6,*) ' a->0 ',x,pr(5),rr
		     icamk(ih,i,2) = 0
		     idau = idau - 1	
                     idfree=idfree+1
		  endif
	      endif
c
   	      if (ici.eq.4) then
		   if(rr.lt.pr(7)) then 
                        write(6,*) ' c->0 ',x,pr(7),rr
	              icamk(ih,i,2) = 0
		      idcap = idcap - 1	
                      idfree=idfree+1
		   endif
		   if(rr.ge.pr(7).and.rr.lt.pr(7)+pr(10)) then 
c                       write(6,*) ' c->a ',x,pr(7),pr(7)+pr(10),rr
	              icamk(ih,i,2) = 3
		      idcap = idcap - 1	
		      idau  = idau  + 1	
		   endif
       	      endif 
c
c             if (ici.eq.0) then
c                pp=(icamb(ih)+icamtp(ih)+icama(ih)+icamcp(ih))/5
c                if(rr.lt.pr(6)*pp*0.1) then
c                     write(6,*) ' 0->c ',x,pr(6),rr
c                   icamk(ih,i,2)=4
c                   idfree=idfree-1
c                   idcap=idcap + 1
c                endif
c             endif
c
c
 21        continue
c
 	iTcamf = iTcamf + idfree
 	icamf(ih) = icamf(ih) + idfree
 	iTcamb = iTcamb + idbnd
 	icamb(ih) = icamb(ih) + idbnd
 	iTcamtp = iTcamtp + idtrap
 	icamtp(ih) = icamtp(ih) + idtrap 
 	iTcama = iTcama + idau
 	icama(ih) = icama(ih) + idau 
 	iTcamcp = iTcamcp + idcap
 	icamcp(ih) = icamcp(ih) + idcap 
c       if(x.gt.373.499.and.ih.lt.50) then
c           write(6,*) x,ih,icamb(ih),iTcamb,idbnd   
c       endif
c
 11	continue
c
 5000   continue
c
ccccccccccccccccccccccc   iccall goto goes to here ccccccccccccc	
c
 	aveTp=iTcamtp/iTholo
c
        Do 93 i=1,isu+1
 	   ibd(i)=0
 	   itp(i)=0
 	   iau(i)=0
 	   icp(i)=0
           ifr(i)=0
 93     continue
c
c       if(x.gt.117.7) write(6,*) ' start 747   '
        Do 747 ih=1, iTholo 
 	   iqtp=0
           iqcp=0
           iqau=0
           iqbd=0
           iqfr=0
c
           i=isu
c           if(x.gt.117.7) write(6,*) ' start free ',ih,icamf
  	   Do while (iqfr.eq.0)
 	      if (icamf(ih).eq.i) then
  	         ifr(i+1)=ifr(i+1)+1
  	         iqfr=1
 	      endif
 	      i=i-1
              if(i.lt.-1) write(6,*) ' qfr ',i,icamf(ih)
 	   enddo 
c           if(x.gt.117.7) write(6,*) ' end   free ',ih,icamb
c
           i=0
  	   Do while (iqbd.eq.0)
 	      if (icamb(ih).eq.i) then
  	         ibd(i+1)=ibd(i+1)+1
  	         iqbd=1
 	      endif
 	      i=i+1
              if(i.gt.isu+1) then
                   write(6,*) x,' qbd ',i,icamb(ih)
                   pause
              endif
 	   enddo 
c           if(x.gt.117.7) write(6,*) ' end   bound ',ih,icamtp
c
           i=0
 	   Do while (iqtp.eq.0)
 	      if (icamtp(ih).eq.i) then
 	         itp(i+1)=itp(i+1)+1
 	         iqtp=1
 	      endif
 	      i=i+1
              if(i.gt.isu+1) write(6,*) ' qtp ',i,icamtp(ih)
 	   enddo 
c           if(x.gt.117.7) write(6,*) ' end   trap  ',ih,icama
c
           i=0
           Do while (iqau.eq.0)
              if (icama(ih).eq.i) then
 	         iau(i+1)=iau(i+1)+1
 	         iqau=1
 	      endif
 	      i=i+1
              if(i.gt.isu+1) write(6,*) ' qau ',i,icama(ih)
 	   enddo 
c           if(x.gt.117.7) write(6,*) ' end   auto  ',ih,icamcp
c
           i=0
 	   Do while (iqcp.eq.0)
 	       if (icamcp(ih).eq.i) then
 	          icp(i+1)=icp(i+1)+1
 		  iqcp=1
 	       endif
 	       i=i+1
              if(i.gt.isu+1) write(6,*) ' qcp ',i,icamcp(ih)
 	   enddo 
c           if(x.gt.117.7) write(6,*) ' end   capp  '
 747	continue
c
c       if(x.gt.117.7) write(6,*) ' exit camkin '
c     do 8766 jklm=1,itholo
c         if(icamb(jklm).lt.0) then
c            write(6,*) ' camkin end ',x,jklm,icamb(jklm),iTcamb
c         endif
c8766 continue
	return
	end
c
c ********************************
        subroutine RND(rr)
        real*8 rr, norm 
        integer m1,m2,a12,q12,r12,a13,q13,r13,a21,q21,r21,
     c     a23,q23,r23, h, p12, p13, p21, p23, 
     c     x10, x11, x12, x20, x21, x22, irr
        parameter (M1=2147483647, M2=2145483479, 
     c     a12 =   63308,  q12 = 33921,  r12 = 12979,
     c     a13 = -183326,  q13 = 11714,  r13 =  2883,
     c     a21 =   86098,  q21 = 24919,  r21 =  7417,
     c     a23 = -539608,  q23 =  3976,  r23 =  2071,
     c     norm=4.656612873077393d-10) 
        data  x10, x11, x12, x20, x21, x22 / 72345, 
     c     985, 12345,  9128, 85, 34567/
        SAVE
c
        h=x10/q13
        p13=-a13*(x10 - h * q13) - h*r13
        h=x11/q12
        p12=a12*(x11 - h * q12) - h*r12
        if(p13.lt.0) p13=p13+m1
        if(p12.lt.0) p12=p12+m1
        x10=x11
        x11=x12
        x12=p12-p13
        if(x12.lt.0) x12=x12+m1
c
        h=x20/q23
        p23=-a23*(x20 - h * q23) - h*r23
        h=x22/q21
        p21=a21*(x22 - h * q21) - h*r21
        if(p23.lt.0) p23=p23+m2
        if(p21.lt.0) p21=p21+m2
        x20=x21
        x21=x22
        x22=p21-p23
        if(x22.lt.0) x22=x22+m2
c
        if(x12.le.x22) then
           irr=x12-x22+m1
        else
           irr=x12-x22
        endif
        rr=irr*norm
c
        return
        END
c**********************************


