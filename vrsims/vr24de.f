      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION V(16),G(991),AL(990),R(990),alamda(990),
     C  D(990),RM(990),RI(990),CM(990),GG(991),G0(991),ipt(9),
     c   JA(3990),IA(992),IRR(991),ICC(991),IC(991),
     c   inode(991),NDAU(990),IPRT(990),ID(4,990),ISOMAD(20),
     C   SIGNUM(990),RMM(990),GP(990),GPR(990),irepl(990),
     c   spta(990),ar(990),fac(990),block(990),
     C   IGPC(6,990),VRP(6),tp(6,300),ITP(6),ITYP(6),
     c   ALF(6),BE(6),A1(6),A2(6),A1A(6),A2A(6),A1BK(6),A2BK(6),
     C   SA2B2(6),X0(6),X0MAX(6),
     c   SCC(6),icstyp(6),vrs(6),gbar(6,990),
     c   exm(6),exh(6),taum(6,1500),tauh(6,1500),
     c   aminf(6,1500),ahinf(6,1500),igs(6,990),am(6,990),
     c   ah(6,990),cai(990),tabpc(6,1500000),sopen(6,990),sblk(6,990),
     C   ARX(6),ARS(6),AR0(6),AR0S(6),YG(6,9),xg(6,9),inf(3,990),
     C   arxr(6,990),ar0r(6,990),ar0sr(6,990),arsr(6,990),tar(3,990),
     c   irc(6,3),stt(3),ett(3),fr(3),kks(3,990),nsyn(3),nas(3)
c      real random
       external GETDM,CHGSPC,CURIP,CONDIP,prtel,linv,GETPQR,
     c   GETRIN,uprcon
       CHARACTER*10 FILE1,file2,file3,file4,fillst
       real rr
c
       open(file='runls6',unit=39,status='old')
 1     read(39,59,end=9999) fillst
 59    format(a10)
c
       MAXIEL=3990
       NC=6
       NSC=6
       NRI=3
       iw=32
       ir=33
       open(file=fillst,unit=ir,status='old')
       open(file='VLASOT',unit=iw,status='unknown')
c      open(file='RAN',unit=iw,status='unknown')
       write(iw,*) ' max number of segments is 990 '
       write(iw,*) ' max # of segment conductances is ',NSC
       write(iw,*) ' max # of point conductances is ',NC
       is1=31
c      is2=37
c      is3=38
       is4=51
c      is5=52
       is6=53
       is7=54
       is8=55
       is9=56
       told=0.d0
       t=0.d0
       write(iw,*) ' enter number of segments '
       read(ir,*) nseg
       neq=nseg+1
       INMAX=nseg+1
       write(iw,*) ' enter n (even) for inverse Laplace ',
     c   ' Transform routine (14 is a good choice) '
       read(ir,11) n
       write(iw,*) ' n= ',n
       m=n+1
       write(iw,*) ' enter MG conc in mM (if NMDA conductance exists)'
       read(ir,*) amg
       write(iw,*) ' enter initial resting potential '
       read(ir,*) brp
       write(iw,*) ' initial resting potential is ',brp
       PI=4.d0*datan(1.d0)
       NVC=0
       NIP=0
       NCIP=0
       NCURLO=0
       NCOND=0
       ICOND=0
       IICOND=0
       ISAVE=0
       h=0.d0
       DO 881 J=1,NC
          do 1234 xxxy=1,1500000
             tabpc(j,xxxy)=0.d0
 1234    continue

          VRP(J)=0.d0
          VRS(J)=0.d0
          ITYP(J)=0
          ITP(J)=0
          ICSTYP(J)=0
          AR0(J)=0.d0
          A1A(J)=0.d0
          A2A(J)=0.d0
          A1(J)=0.d0
          A2(J)=0.d0
          ALF(J)=0.d0
          BE(J)=0.d0
          AR0S(J)=0.d0
          A1BK(J)=0.d0
          A2BK(J)=0.d0
          SA2B2(J)=0.d0
          X0(J)=0.d0
          X0MAX(J)=0.d0
          SCC(J)=0.d0
          ARX(J)=0.d0
          ARS(J)=0.d0
          DO 818 I=1,990
             ARXR(J,I)=0
             AR0R(J,I)=0
             AR0SR(J,I)=0
             ARSR(J,I)=0
 818      CONTINUE
          do 821 m=1,nri
             irc(j,m)=0
             nsyn(m)=0
             nas(m)=0
 821      continue
 881   CONTINUE
       DO 819 I=1,990
          do 820 m=1,nri
           INF(m,I)=0
           kks(m,i)=0
           TAR(m,I)=0
 820      continue
 819   CONTINUE
       do 882 j=1,nc
          do 882 i=1,9
             xg(j,i)=0.d0
             yg(j,i)=0.d0
 882   continue
 879   DO 898 I=1,nseg
           G(I)=0.d0
           GP(I)=0.d0
           GPR(I)=0.d0
           RMM(I)=0.d0
           SIGNUM(I)=0.d0
           ALAMDA(I)=0.d0
           G0(I)=0.d0
           GG(I)=0.d0
           DO 700 J=1,NC
               IGPC(J,I)=0
               igs(J,I)=0
               am(j,i)=0.d0
               ah(j,i)=0.d0
               gbar(j,i)=0.d0
 700       CONTINUE
 898   CONTINUE
       g(neq)=0.d0
       gg(neq)=0.d0
       g0(neq)=0.d0
c
       WRITE(iw,*) ' SAVE FILE? 1=YES '
       READ(ir,11) NN
       WRITE(iW,11) NN
       IF(NN.EQ.1) THEN
          WRITE (iw,*) ' ENTER SAVE FILE NAMES XXXXXX'
c         READ(ir,799) FILE1,file2,file3
c799      FORMAT(3A10)
          READ(ir,799) FILE1
 799      FORMAT(A10)
          write(iw,*) ' enter 9 locations to save '
          read(ir,*) ipt
          write(iw,*) ipt
c          write(iw,*) ' enter location to monitor conductance '
c          read(ir,*) imo
c          write(iw,*) imo
c         OPEN (UNIT=is1,FILE=FILE1,status='unknown')
c         OPEN (UNIT=is2,FILE=FILE2,status='unknown')
c         OPEN (UNIT=is3,FILE=FILE3,status='unknown')
          WRITE (iw,*) ' ENTER loc1 curr FILE NAMES '
          READ(ir,798) FILE1,file2
 798      format(2a10)
c         OPEN (UNIT=is4,FILE=FILE1,status='unknown')
c         OPEN (UNIT=is5,FILE=FILE2,status='unknown')
          WRITE (iw,*) ' ENTER loc2 curr FILE NAMES '
          READ(ir,797) FILE1,file2,file3,file4
 797      FORMAT(4A10)
c         OPEN (UNIT=is6,FILE=FILE1,status='unknown')
          OPEN (UNIT=is7,FILE=FILE2,status='unknown')
c         WRITE (iw,*) ' ENTER loc3 curr FILE NAMES '
c         READ(ir,798) FILE1,file2
c         OPEN (UNIT=is8,FILE=FILE3,status='unknown')
c         OPEN (UNIT=is9,FILE=FILE4,status='unknown')
          ISAVE=1
       ENDIF
c
       CALL GETDM(D,AL,RI,RM,CM,NSEG,iw,ir,isn,ish,idsp,
     c   JA,IA,IRR,ICC,IC,NEQ,MAXIEL,IEL,fac,ar,spta,irepl,
     c   inode,inmax,ndau,IPRT,ID,ISOMAD,IND)
c
 122   WRITE(iw,822)
 822   FORMAT(' DO YOU WANT TO CHANGE ANY SPECS 1=YES')
       READ(ir,11) NN
       WRITE(iW,11) NN
 11    FORMAT(I2)
       IF (NN.EQ.1) CALL CHGSPC(D,AL,RI,RM,CM,ar,fac,spta,irepl,
     c   id,nseG,IW,IR,PI,isn,ish,idsp)
c
       write(iw,91)
  91   format(' do you want TO CHANGE SEGMENT conductance input 1=yes')
       read(ir,11) nn
       WRITE(iW,11) nn
       if(nn.eq.1) CALL CONDIS(NCOND,NSEG,NSC,igs,vrs,icstyp,NperMv,
     c   gbar,exm,exh,taum,tauh,aminf,ahinf,brp,cai,iw,ir)
       dNper=dfloat(NperMv)
       temp=dNper*(brp+100.d0)+1.d0
       ibrp=temp
c
       write(iw,92)
  92   format(' do you want TO CHANGE POINT conductance input 1=yes')
       read(ir,11) nn
       WRITE(iW,11) nn
       if(nn.eq.1) CALL CONDIP(NCOND,IGPC,VRP,NSEG,NC,tp,brp,
     c   ITP,ITYP,ALF,BE,A1,A2,A1A,A2A,A1BK,A2BK,SA2B2,X0,
     c   X0MAX,SCC,tabpc,iw,ir)
c
       write(iw,93)
 93    format('do you want to have RANDOM input 1=yes')
       read(ir,11) nnr
       write(iw,11) nnr
       if(nnr.eq.0) goto 709
       call getrin(ir,iw,nri,stt,ett,fr,nc,igpc,kks,nsyn,irc,nseg)
c
 709   DO 1910 I=1,NC
          AR0(I)=X0(I)
          AR0S(I)=0.d0
          do 1916 m=1,nri
            DO 1915 J=1,NSYN(m)
              AR0R(I,J)=0.D0
              AR0SR(I,J)=0.D0
 1915       CONTINUE
 1916     continue
 1910  CONTINUE
c
C
 2000  WRITE(iw,1400)
 1400  FORMAT(' DOES CONDUCTANCE INPUT TOTALLY DEFINE RM 1=YES')
       READ(ir,11) NDETRM
       WRITE(iW,11) NDETRM
       DO 1410 I=1,NSEG
          R(I)=RI(I)*4.d0/(PI*D(I)*D(I))
          SIGNUM(I)=0.0d0
          CTOT=0.0d0
          DO 1420 J=1,NSC
             if(igs(j,i).eq.0) goto 1420
             if(icstyp(j).eq.1) cond=gbar(j,i)
             if(icstyp(j).eq.2.or.icstyp(j).eq.3.or.
     c         icstyp(j).eq.5.or.icstyp(j).eq.6) then
                gtab=g0(i)*dNper
                it=gtab+ibrp
c               write(6,*) 'gtab,it= ',gtab,it
c  ibrp is computed from the resting potential assumed uniform
c
                am(j,i)=aminf(j,it)
                if(exh(j).lt.1.d-4) then
                   cond=gbar(j,i)*(aminf(j,it)**exm(j))
                else
                   cond=gbar(j,i)*(aminf(j,it)**exm(j)*
     c                        ahinf(j,it)**exh(j))
                   ah(j,i)=ahinf(j,it)
                endif
c               write(6,*) ' am,cond= ',am(j,i),cond
             endif
             if(icstyp(j).eq.4) then
                dlcar=dlog10(0.01d0)
                gtab=(dlog10(cai(i))-dlcar)/4.d0*1500.d0     
                it=idint(gtab)+1
                if(it.lt.1) it=1
                am(j,i)=aminf(j,it)
                cond=gbar(j,i)*(aminf(j,it)**exm(j))
c               write(6,*) ' cam,cond= ',am(j,i),cond,i
             endif
             if(icstyp(j).eq.6) then 
                cminca=cai(i)/2.5
                if(cminca.gt.1.d0) cminca=1.d0
                cond=cond*cminca
             endif
c
             CTOT=CTOT+cond   
             SIGNUM(I)=SIGNUM(I)+VRS(J)*cond
 1420        CONTINUE
          IF(NDETRM.NE.1) RMM(I)=1.0d0/(1.0d0/RM(I)+CTOT)
          IF(NDETRM.EQ.1) RMM(I)=1.0d0/CTOT
          ALAMDA(I)=dSQRT((RMM(I)*D(I))/(4.d0*RI(I)))
          SIGNUM(I)=SIGNUM(I)*RMM(I)
 1410     CONTINUE
C
       WRITE(iw,876)
 876   FORMAT(' PRINT ELECTROTONIC LENGTHS 1=YES')
       READ(ir,11) NN
       WRITE(iW,11) NN
       IF(NN.EQ.1) CALL PRTEL(AL,ALAMDA,RMM,NSEG,IPRT,INODE,IND,IW,IR)
c
       write(iw,2110)
 2110  format(' do you want CURRENT IP OR VOL CLAMP 1=yes')
       READ(ir,11) NCIP
       WRITE(iw,11) NCIP
       h=0.d0
       IF(NCIP.EQ.1) CALL CURIP(H,h2,h3,CUR,cur2,cur3,
     c   NVC,NIP,TEXP,NCURLO,CSTR,CEND,CEND2,IW,IR)
c
       write(iw,2100)
 2100  format(' enter tmax, time incr, cond time inc,neval  ')
       read(ir,*) tmax,timinc,dtc,neval
       WRITE(iw,*) tmax,timinc,dtc,neval
c      close(unit=iw)
       tg2=10.d0
c
       IF(M.NE.N) call linv(n,v,m)
       hh=h
       ISW=0
c
       ress=r(1)*al(1)
c
 1980  continue
c      write(6,*) ' at 1980 '
       t=told+timinc
       if(t.gt.tmax+.0001d0) then
c         write(iw,*) ' end? 1=yes '
          read(ir,11) nn
          WRITE(iW,11) nn
          if(nn.eq.1) goto 9998
c         write(iw,*) ' enter new tmax, time incr '
          read(ir,*) tmax,timinc
          WRITE(iW,*) tmax,timinc
          t=told+timinc
       endif
       h=0.d0
       if(t.lt.cstr+1.d-4) h=0.d0
       if(t.ge.cstr+1.d-4) h=hh
       if(t.ge.cend+1.d-4) h=h2  
       if(t.ge.cend2+1.d-4) h=h3  
c      write(6,*) ' h= ',h
C
       call UPSCON(nseg,nsc,neq,igs,vrs,icstyp,dNper,al,
     c   gbar,exm,exh,taum,tauh,aminf,ahinf,am,ah,ibrp,cai,
     c   g0,signum,alamda,RMM,RM,D,RI,NDETRM,timinc,iprt,iw,ir)
c
       call UPPCON(nc,nseg,neq,gp,gpr,ityp,t,tp,itp,tg2,
     c   arx,ar0,a1a,ar0s,alf,a1,a2a,a2,sa2b2,ars,a1bk,be,a2bk,
     c   x0,x0max,scc,igpc,vrp,xg,yg,g,block,amg,ipt,brp,
     c   tabpc,sopen,sblk,rr,timinc)
c
c
        if(nnr.eq.0) goto 711
       call UPRCON(kks,fr,nsyn,probf,inf,tar,t,be,x0,x0max,
     c      arxr,ar0r,a1a,ar0sr,alf,a1,a2a,a2,sa2b2,a1bk,a2bk,
     C      gnowt,scc,gp,gpr,vrp,block,g,xg,yg,brp,nseg,neq,amg,nc,
     c      arsr,p,k,nri,timinc,irc,ityp,stt,ett,nas,ipt)
c
 711   continue
       TT=T
       IF (IICOND.EQ.1) T=TT-TOLD
       DO 222 I=1,neq
          G(I)=0.d0
  222     CONTINUE
       a=dlog(2.d0)/t
c      a=0.693147180559945309417d0/t
       vcso=0.d0
 4411  do 6 i=1,n
          n3=i
          ai=a*dfloat(i)
 4412     call getpqr(ai,GG,H,alamda,al,RMM,CM,r,cur,
     C         NVC,NIP,TEXP,NCOND,NSEG,NEQ,
     C         G0,vvso,IICOND,NCURLO,NCIP,iw,
     C         SIGNUM,GP,GPR,irepl,
     C         IA,JA,IRR,IC,ICC,IEL,n3,ISW)
          DO 223 II=1,neq
             G(II)=G(II)+V(I)*GG(II)
  223        CONTINUE
          vcso=vcso+v(i)*vvso
  6    continue
       DO 224 II=1,neq
          G(II)=A*G(II)
  224     CONTINUE
       vcso=a*vcso
c
c1210  write(iw,1205) Tt
c1205  format('  time=  ',f10.4)
c
       vccur=(g(neq)-g(1))/ress
       if(tt.lt.cur2) then
          vcss=vccur
          vccur1=0.d0
       else
          vccur1=(vccur-vcss)*1.0d+6
       endif
c      write(6,*) al(1),alamda(1),ress,r(1)
c
 100   CONTINUE
       IF(ISAVE.EQ.1) THEN
c          WRITE(iw,*) ' SAVE THIS DATA? 1=YES '
c          READ(ir,11) NN
c          IF(NN.EQ.1) THEN
c         WRITE(is1,103) TT,g(neq),g(ipt(1)),g(ipt(2)),g(ipt(3)),
c    c      g(ipt(4)),g(ipt(5)),g(ipt(6)),g(ipt(7)),
c    c      g(ipt(8)),g(ipt(9))
c         WRITE(is2,103) TT,g(ipt(1))
c         WRITE(is3,103) TT,g(ipt(2))
          gca1=gbar(3,ipt(1))*am(3,ipt(1))**exm(3)*
     c      ah(3,ipt(1))**exh(3)*(g(ipt(1))-vrs(3))*
     c      pi*d(ipt(1))*al(ipt(1))*1.d+6
          gca2=gbar(3,ipt(2))*am(3,ipt(2))**exm(3)*
     c      ah(3,ipt(2))**exh(3)*(g(ipt(2))-vrs(3))*
     c      pi*d(ipt(2))*al(ipt(2))*1.d+6
          gca3=gbar(3,ipt(3))*am(3,ipt(3))**exm(3)*
     c      ah(3,ipt(3))**exh(3)*(g(ipt(3))-vrs(3))*
     c      pi*d(ipt(3))*al(ipt(3))*1.d+6
          gca4=gbar(3,ipt(4))*am(3,ipt(4))**exm(3)*
     c      ah(3,ipt(4))**exh(3)*(g(ipt(4))-vrs(3))*
     c      pi*d(ipt(4))*al(ipt(4))*1.d+6
          gca5=gbar(3,ipt(5))*am(3,ipt(5))**exm(3)*
     c      ah(3,ipt(5))**exh(3)*(g(ipt(5))-vrs(3))*
     c      pi*d(ipt(5))*al(ipt(5))*1.d+6
          gca6=gbar(3,ipt(6))*am(3,ipt(6))**exm(3)*
     c      ah(3,ipt(6))**exh(3)*(g(ipt(6))-vrs(3))*
     c      pi*d(ipt(6))*al(ipt(6))*1.d+6
          gnas=gbar(1,1)*am(1,1)**exm(1)*ah(1,1)**exh(1)*
     c      (g(1)-vrs(1))*
     c      pi*d(1)*al(1)
          gks=gbar(2,1)*am(2,1)**exm(2)*(g(1)-vrs(2))*
     c      pi*d(1)*al(1)
          gcas=gbar(3,1)*am(3,1)**exm(3)*ah(3,1)**exh(3)*
     c      (g(1)-vrs(3))*
     c      pi*d(1)*al(1)
          gks2=gbar(4,1)*am(4,1)**exm(4)*(g(1)-vrs(4))*
     c      pi*d(1)*al(1)
c         write(6,*) ' after gks2 write '
c         gkkk=gbar(4,1)*am(4,1)**exm(4)
c         write(6,*) gkkk
c         gnaa=gbar(1,ipt(1))*am(1,ipt(1))**exm(1)*
c    c      ah(1,ipt(1))**exh(1)*(g(ipt(1))-vrs(1))
c         gka=gbar(2,ipt(1))*am(2,ipt(1))**exm(2)*(g(ipt(1))-vrs(2))
c         gnad=gbar(1,ipt(2))*am(1,ipt(2))**exm(1)*
c    c      ah(1,ipt(2))**exh(1)*(g(ipt(2))-vrs(1))
c         gkd=gbar(2,ipt(2))*am(2,ipt(2))**exm(2)*(g(ipt(2))-vrs(2))
c         WRITE(is4,104) TT,gca1,gca2,gca3,gca4,gca5,gca6
c         WRITE(is4,104) TT,gnas,gks,gks2,gcas,cai(1)
c         WRITE(is5,104) TT,gnaa,gka,gnad,gkd,cai(4)
c         WRITE(is6,104) TT,gks2
c         WRITE(is7,104) TT,gcas 
c         WRITE(is8,104) TT,cai(1)
c         WRITE(is9,104) TT,gka 
 104      format(f10.3,6e13.5)
          syncur=yg(1,1)+yg(2,1)
c         WRITE(is6,105) TT,g(neq),g(1),vcso,vccur1,vccur,G(ipt(1)),
c    c         xg(1,1),xg(2,1),yg(1,1),yg(2,1),syncur
c         WRITE(is6,105) TT,g(neq),G(ipt(1)),xg(1,1),xg(2,1),yg(2,1),
c    c         G(ipt(2)),xg(3,2),xg(4,2),yg(4,2)
          WRITE(is7,105) TT,g(neq),G(ipt(3)),xg(1,3),xg(3,3),yg(3,3),
     c         G(ipt(4)),xg(1,4),xg(2,4),yg(2,4)
c         WRITE(is8,105) TT,g(neq),G(ipt(5)),xg(1,5),xg(2,5),yg(2,5),
c    c         G(ipt(6)),xg(1,6),xg(2,6),yg(2,6)
c         WRITE(is9,105) TT,g(neq),G(ipt(7)),xg(1,7),xg(2,7),yg(2,7),
c    c         G(ipt(8)),xg(1,8),xg(2,8),yg(2,8)
c         WRITE(is7,105) TT,g(neq),(G(ipt(i)),xg(1,i),yg(1,i),I=4,6)
c         WRITE(is8,105) TT,g(neq),(G(ipt(i)),xg(4,i),yg(4,i),I=7,9)
 103         FORMAT(F9.3,10F9.3)
 105         FORMAT(F9.3,2F9.4,3e13.4,f9.4,5e13.5)
c          ENDIF
       ENDIF
c
c1220  WRITE(iw,3020)
c3020  FORMAT(' USE THESE VALUES AS IC FOR NEXT RUN 1=YES')
c      READ(ir,11) ICOND
c      IF (ICOND.NE.1) GOTO 990
          icond=1
          TOLD=TT
          IICOND=ICOND
          DO 8000 I=1,NEQ
             G0(I)=G(I)
 8000        CONTINUE
c
c990   write(iw,3010)
c990   write(6,*) ' do you only want another time point '
c3010  format(' do you only want another time point? 1=yes')
c      read(ir,11) nn
c      if(nn.ne.1) goto 122
c
c      write(iw,3011)
c3011  format(' enter new time point f10.5')
c      read(ir,*) t
       goto 1980
c
 9998  continue
       do 351 i=1,nri
       write(iw,*) nas(i),'synapses were active for',
     c  ' random input type',i
 351   continue
       close(unit=ir)
       close(unit=iw)
       close(unit=is1)
c      close(unit=is2)
c      close(unit=is3)
       close(unit=is4)
c      close(unit=is5)
       close(unit=is6)
       close(unit=is7)
       close(unit=is8)
       close(unit=is9)
       goto 1
 9999  stop
       end
c ******************************************************
       SUBROUTINE UPSCON(nseg,nsc,neq,igs,vrs,icstyp,dNper,al,
     c   gbar,exm,exh,taum,tauh,aminf,ahinf,am,ah,ibrp,cai,
     c   g0,signum,alamda,RMM,RM,D,RI,NDETRM,timinc,iprt,iw,ir)
c
       IMPLICIT real*8 (A-H,O-Z)
c
       DIMENSION IGS(NSC,NSEG),VRS(nsc),icstyp(NSC),gbar(NSC,NSEG),
     C   exm(NSC),exh(NSC),taum(NSC,1500),tauh(NSC,1500),
     C   aminf(NSC,1500),ahinf(NSC,1500),am(NSC,NSEG),ah(NSC,NSEG),
     C   g0(NEQ),signum(NSEG),alamda(NSEG),RMM(NSEG),RM(NSEG),
     C   D(NSEG),RI(NSEG),cai(nseg),al(nseg),iprt(nseg)
c
       pi=4.d0*datan(1.d0)
       dlcar=dlog10(0.01d0)
       DO 1410 I=1,NSEG
          SIGNUM(I)=0.0d0
          CTOT=0.0d0
          DO 1420 J=1,NSC
             if(igs(j,i).eq.0) goto 1420
             if(icstyp(j).eq.1) cond=gbar(j,i)
             if(icstyp(j).eq.2.or.icstyp(j).eq.3) then
c               k=i-1
c               if(k.lt.1) k=neq
                k=iprt(i)
                gtab=(g0(k)+g0(i))*0.5d0*dNper
                it=idint(gtab)+ibrp
c               if(it.gt.1500) then
c                  write(iw,*) ' it = ',it
c                  it=1500
c               endif
                am(j,i)=aminf(j,it)+(am(j,i)-aminf(j,it))*
     c               dexp(-timinc/taum(j,it))
                if(exh(j).lt.1.d-4) then
                   cond=gbar(j,i)*(am(j,i)**exm(j))
                else
                   ah(j,i)=ahinf(j,it)+(ah(j,i)-ahinf(j,it))*
     c                dexp(-timinc/tauh(j,it))
                   cond=gbar(j,i)*(am(j,i)**exm(j))*
     c                       (ah(j,i)**exh(j))
                endif
c               write(6,*) ' uam,cond= ',am(j,i),cond
             endif
             if(icstyp(j).eq.3) then
c               k=i-1
c               if(k.lt.1) k=neq
                k=iprt(i)
                cacur=gbar(j,i)*(am(j,i)**exm(j))*(ah(j,i)**exh(j))*
     c             ((g0(k)+g0(i))*0.5d0-vrs(j))*
     c             pi*d(i)*al(i)
                cleak=gbar(6,i)* 0.027546d0*
     c             ((g0(k)+g0(i))*0.5d0-vrs(6))*
     c             pi*d(i)*al(i)
c                     current in micro-amps = 10-6 coul/sec
c               vol=pi*al(i)*(0.1d-4*d(i)-0.01d-8)
c                     volume=pi*l*0.25*(d**2-(d-.2d-4)**2) cm**3
                vol=pi*al(i)*d(i)*0.2d-4  
c                     durand vol=area*.2d-4
                dcacur=(cacur+cleak)/(2.d0*9.648456d4*vol)
c               dcacur=cacur*timinc/(2.d0*9.648456d4*vol)
c                     F=9.6484564e4 coul/mol  z=2
c                     num=10-6coul/sec*sec/10e3ms*timeinc
c                     denom=zF coul/mol*mol/10e6umol*vol cm3*l/10e3ml
c               cai(i)=cai(i)-(cai(i)-0.07d0)/9.d0-dcacur
                cai(i)=cai(i)-timinc*(cai(i)/9.d0+dcacur)
                if(cai(i).lt.0.07d0) cai(i)=0.07d0
c               cai(i)=cai(i)-timinc*((cai(i)-0.07d0)/9.d0+dcacur)
c               write(6,*) ' cacur,dcacur,cai ',cacur,dcacur,cai(i),i
             endif
c
             if(icstyp(j).eq.4) then
                gtab=(dlog10(cai(i))-dlcar)/4.d0*1500.d0     
                it=idint(gtab)+1
                if(it.lt.1) it=1
c               write(6,*) ' gtab,it ',gtab,it
                am(j,i)=aminf(j,it)+(am(j,i)-aminf(j,it))*
     c               dexp(-timinc/taum(j,it))
                cond=gbar(j,i)*(am(j,i)**exm(j))
c               write(6,*) ' ucam,cond= ',am(j,i),cond,it,aminf(j,it)
c               if(i.eq.1) write(6,8888) am(j,i),cond,cai(i),dcacur
c8888           format(4d14.5)
             endif
c
             CTOT=CTOT+cond   
             SIGNUM(I)=SIGNUM(I)+VRS(J)*cond
 1420        CONTINUE
          IF(NDETRM.NE.1) RMM(I)=1.0d0/(1.0d0/RM(I)+CTOT)
          IF(NDETRM.EQ.1) RMM(I)=1.0d0/CTOT
          ALAMDA(I)=dSQRT((RMM(I)*D(I))/(4.d0*RI(I)))
          SIGNUM(I)=SIGNUM(I)*RMM(I)
 1410     CONTINUE
C
       return
       end
c ********************************************************
       SUBROUTINE UPPCON(nc,nseg,neq,gp,gpr,ityp,t,tp,itp,tg2,
     c   arx,ar0,a1a,ar0s,alf,a1,a2a,a2,sa2b2,ars,a1bk,be,a2bk,
     c   x0,x0max,scc,igpc,vrp,xg,yg,g,block,amg,ipt,brp,
     c   tabpc,sopen,sblk,rr,timinc)
       IMPLICIT real*8(A-H,O-Z)
c
       DIMENSION gp(nseg),gpr(nseg),ityp(nc),tp(nc,300),
     C   ITP(NC),arx(nc),ar0(nc),a1a(nc),ar0s(nc),
     C   ALF(NC),A1(NC),A2A(NC),A2(NC),sa2b2(nc),ars(nc),
     C   A1BK(NC),be(nc),A2BK(NC),X0(NC),X0MAX(NC),SCC(NC),
     C   igpc(NC,nseg),vrp(nc),xg(nc,9),yg(nc,9),g(neq),
     c   sopen(nc,nseg),sblk(nc,nseg),
     C   block(nseg),ipt(9),tabpc(nc,1500000)
c
      external RAND
      real rr
c
       dt=timinc
       DO 298 I=1,NSEG
          GP(I)=0.d0
          GPR(I)=0.d0
 298      CONTINUE
       DO 300 I=1,NC
          if(ityp(i).eq.0) goto 300
          if(t.lt.tp(i,1)+1.d-8) goto 300
C
          if(ityp(i).eq.1.or.ityp(i).eq.2.or.ityp(i).eq.3) then
             tent=t-tp(i,1)
             ilook=idint(tent*100d0+1.00001d0)
             if(ilook.gt.1500000) then 
                ars(i) = 0
             else  
                ars(i)=tabpc(i,ilook)
             endif
          else if(ityp(i).eq.4.or.ityp(i).eq.14) then
             if(t.ge.tp(i,1).and.t.le.tp(i,2)) then
               ars(i)=x0(i)
             else
               ars(i)=0.d0
             endif
          else
             do 320 k=1,itp(i)
                tg1=t-tp(i,k)
                if(tg1.gt.-1.d-8) then
                    tg=tg1
                    ksav=k
                endif
 320         continue
C
             if(dabs(tg).lt.1.d-8) then
                tg2=tg
                if(itp(i).gt.1) tg=tp(i,ksav)-tp(i,ksav-1)
             endif
c
             ARX(I)=( (AR0(I)*A1A(I)+AR0S(I)*ALF(I))*dEXP(A1(I)*TG) -
     C            (AR0(I)*A2A(I)+AR0S(I)*ALF(I))*dEXP(A2(I)*TG) ) /
     C             SA2B2(I)
             ARS(I)=( (AR0S(I)*A1BK(I)+AR0(I)*BE(I))*dEXP(A1(I)*TG) -
     C             (AR0S(I)*A2BK(I)+AR0(I)*BE(I))*dEXP(A2(I)*TG) ) /
     C             SA2B2(I)
c
             if(dabs(tg2).lt.1.d-8) then
                tg2=10.0
                AR0S(I)=ARS(I)
                AR0(I)=ARX(I)+X0(I)*(x0max(i)-arx(i)-ars(i))/x0max(i)
c               if(ar0(i).gt.600) ar0(i)=600
             endif
          endif
c
          if(ityp(i).eq.3) then
             do 720 j=1,nseg
                if(igpc(i,j).eq.0) goto 720
                if(ars(i).gt.sopen(i,j)+sblk(i,j)+1.d-4) then
                   iadd=idint(ars(i)-sopen(i,j)-sblk(i,j)+1.d-4)
                   pblk=amg*0.61d0*dexp(-(g(j)+brp)/17.d0)*dt
                   punb=5.4d0*dexp((g(j)+brp)/47.d0)*dt
                   pblk=1.0-dexp(-pblk)
                   punb=1.0-dexp(-punb)
                   do 820 k=1,iadd
                     call RAND(rr)
                     rr=(pblk+punb)*rr
                     if(rr.ge.pblk) sopen(i,j)=sopen(i,j)+1.d0
                     if(rr.lt.pblk) sblk(i,j)=sblk(i,j)+1.d0
  820              continue
                endif
                if(ars(i).lt.sopen(i,j)+sblk(i,j)-1.d-4) then
                   iadd=idint(sopen(i,j)+sblk(i,j)-ars(i)+1.d-4)
                   do 830 k=1,iadd
                     if(sopen(i,j).gt.1d-4.and.sblk(i,j).gt.1.d-4)
     c                 then
                        popen=sopen(i,j)/(sopen(i,j)+sblk(i,j)+1.d-5)
                        call RAND(rr)
                        if(rr.lt.popen) sopen(i,j)=sopen(i,j)-1.d0
                        if(rr.ge.popen) sblk(i,j)=sblk(i,j)-1.d0
                     else if (sopen(i,j).gt.1.d-4.and.sblk(i,j)
     c                .lt.1.d-4) then
                        sopen(i,j)=sopen(i,j)-1.d0
                     else if (sopen(i,j).lt.1.d-4.and.sblk(i,j)
     c                .gt.1.d-4) then
                        sblk(i,j)=sblk(i,j)-1.d0
                     endif
  830              continue
                endif
c               if(ars(i).ne.0.) then
c               write(6,*) ' 3 ',ars(i),sopen(i,j),sblk(i,j)
c               endif
  720        continue
          endif
c
          GNOWT=ARS(I)*SCC(I)*1.0d-6
c
          IF(ITYP(I).lt.2.or.ityp(i).gt.4) THEN
             DO 340 J=1,NSEG
                GP(J)=GP(J)+GNOWT*dFLOAT(IGPC(I,J))
                GPR(J)=GPR(J)+GNOWT*dFLOAT(IGPC(I,J))*VRP(I)
 340         continue
             do 342 l=1,9
                xg(i,l)=gnowt*1.0d+6*dfloat(IGPC(i,ipt(l))) 
                YG(I,l)=-XG(I,l)*(G(ipt(l))-VRP(I))
 342         continue
          endif
          IF(ITYP(I).EQ.2.or.ityp(i).eq.4) THEN
             do 349 j=1,nseg
                BLOCK(j)=8.8d0*dEXP((G(J)+brp)/12.5d0)/
     c              (amg+8.8d0*dEXP((G(J)+brp)/12.5d0))
                GP(J)=GP(J)+GNOWT*dFLOAT(IGPC(I,J))*BLOCK(j)
                GPR(J)=GPR(J)+GNOWT*dFLOAT(IGPC(I,J))*VRP(I)*BLOCK(j)
 349         continue
           endif
           IF(ITYP(I).eq.3) then
             do 345 j=1,nseg
                if(igpc(i,j).ne.0) then
                   blk=0.d0
                   unb=0.d0
                   pblk=amg*0.61d0*dexp(-(g(j)+brp)/17.d0)*dt
                   punb=5.4d0*dexp((g(j)+brp)/47.d0)*dt 
                   pblk=1.0-dexp(-pblk) 
                   punb=1.0-dexp(-punb)
c               if(ars(i).ne.0.) then
c               write(6,*) ' 4 ',ars(i),sopen(i,j),sblk(i,j)
c               endif
                   do 400 k=1,idint(sopen(i,j)+1.d-4)
                      call RAND(rr)
                      if(rr.lt.pblk) blk=blk+1.d0
 400               continue
                   do 405 k=1,idint(sblk(i,j)+1.d-4)
                      call RAND(rr)
                      if(rr.lt.punb) unb=unb+1.d0
 405               continue
                   sopen(i,j)=sopen(i,j)-blk+unb
                   sblk(i,j)=sblk(i,j)-unb+blk
                   gnowt=sopen(i,j)*scc(i)*1.d-6
                   GP(J)=GP(J)+GNOWT
                   GPR(J)=GPR(J)+GNOWT*VRP(I)
c               if(ars(i).ne.0.) then
c               write(6,*) ' 5 ',ars(i),sopen(i,j),sblk(i,j)
c               endif
                endif
 345         continue
          endif
          if(ityp(i).gt.1.and.ityp(i).lt.5) then
             do 346 l=1,9
                xg(i,l)=gnowt*1.0d+6*dfloat(igpc(i,ipt(l)))
                if(ityp(i).ne.3) xg(i,l)=xg(i,l)*BLOCK(ipt(l))
                YG(I,l)=-XG(I,l)*(G(ipt(l))-VRP(I))
 346         continue
          ENDIF
 300   CONTINUE
c
       return
       end
C
c *****************************
       SUBROUTINE UPRCON(kks,fr,nsyn,probf,inf,tar,t,be,x0,x0max,
     c      arxr,ar0r,a1a,ar0sr,alf,a1,a2a,a2,sa2b2,a1bk,a2bk,
     C      gnowt,scc,gp,gpr,vrp,block,g,xg,yg,brp,nseg,neq,amg,nc,
     c      arsr,p,k,nri,timinc,irc,ityp,stt,ett,nas,ipt)
c
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      real random
       real rr
c
       DIMENSION gp(nseg),gpr(nseg),tar(nri,nseg),ARXR(nc,nseg),
     C  a1a(nc),ALF(nc),A1(nc),A2A(nc),A2(nc),ar0sr(nc,nseg),sa2b2(nc),
     C  A1BK(nc),be(nc),A2BK(nc),X0(nc),X0MAX(NC),
     c  SCC(nc),vrp(nc),g(neq),
     C  kks(3,990),tgr(990),arsr(nc,nseg),ar0r(nc,nseg),ityp(nc),
     c  block(nseg),fr(nri),nsyn(nri),inf(nri,nseg),irc(nc,nri),
     c  stt(nri),ett(nri),nas(nri),xg(nc,9),yg(nc,9),ipt(9)
       external RAND
c
       Pi=3.1415926535897932384626433832795
       do 604 m=1,nri
         if(t.lt.stt(m)) goto 617
         if(t.gt.ett(m)) goto 617
         probf=fr(m)*timinc/1000.d0
c        probf=fr(m)*timinc/(1000.d0*nsyn(m))
         refr=0.0
         do 602 i=1,nsyn(m)
           call rand(rr)
           if(rr.lt.probf) then
              IF(inf(m,i).EQ.1) THEN
                if (tar(m,i).lt.(t-refr)) then
                  tar(m,i)=t
                  nas(m)=nas(m)+1
                endif
              ELSE
                TAR(M,I)=T
                inf(m,i)=1
                nas(m)=nas(m)+1
              ENDIF
           endif
 602     continue
c
 617     continue
         do 603 i=1,nsyn(m)
           IF(inf(m,i).EQ.1) then
              tgr(i)=t-tar(m,i)
              do 609 j=1,nc
                 if(irc(j,m).eq.0) goto 609
                 if (dabs(tgr(i)).lt.1d-8) then
                     AR0Sr(j,i)=ARSr(j,i)
             AR0r(j,i)=ARXr(j,i)+X0(j)*(x0max(j)-arxr(j,i)-
     c                    arsr(j,i))/x0max(j)
                  endif
            ARXR(j,i)=( (AR0r(j,i)*A1A(j)+AR0Sr(j,i)*ALF(j))*EXP(A1(j)*
     C        TGR(i))-(AR0r(j,i)*A2A(j)+AR0Sr(j,i)*ALF(j))*EXP(A2(j)*
     C        TGR(i)) ) / SA2B2(j)
            ARSr(j,i)=( (AR0Sr(j,i)*A1BK(j)+AR0r(j,i)*BE(j))*EXP(A1(j)*
     C        TGR(i))-(AR0Sr(j,i)*A2BK(j)+AR0r(j,i)*BE(j))*EXP(A2(j)*
     C        TGR(i)) ) / SA2B2(j)
c
                  GNOWT=ARSr(j,i)*SCC(j)*1.0d-6
c
                  IF(ITYP(J).ne.2) THEN
                    GP(kks(m,i))=GP(kks(m,i))+GNOWT
                    GPR(kks(m,i))=GPR(kks(m,i))+GNOWT*VRP(j)
*
                  do  173  l=1,9
                     if(kks(m,i).ne.ipt(l)) goto 173
                     xg(j,l)=gnowt*1.0d+6
                     YG(J,l)=-XG(J,l)*(G(ipt(l))-VRP(J))
 173               continue
*
                  ENDIF
c
                  IF (ITYP(J).EQ.2) THEN
            BLOCK(kks(m,i))=8.8D0*DEXP((G(kks(m,i))+brp)/12.5D0)/
     c           (amg+8.8D0*DEXP((G(kks(m,i))+brp)/12.5D0))
            GP(kks(m,i))=GP(kks(m,i))+GNOWT*BLOCK(kks(m,i))
            GPR(kks(m,i))=GPR(kks(m,i))+GNOWT*VRP(J)*BLOCK(kks(m,i))
*
                  do  174  l=1,9
                     if(kks(m,i).ne.ipt(l)) goto 174
                     xg(j,l)=gnowt*1.0d+6*BLOCK(ipt(l))
                     YG(J,l)=-XG(J,l)*(G(ipt(l))-VRP(J))
 174               continue
*
                  ENDIF
609           continue
            ENDIF
603         continue
604     CONTINUE
        return
        END
c *****************************
        subroutine RAND(rr)
        INTEGER K,M,CONST1
        real rr,CONST2
        parameter (k=789,CONST1=2147483647, CONST2=.4656613e-9)
        SAVE
        DATA M/0/
c
        If (M.EQ.0) M=K
        m=m*65539
        if(m.lt.0) m=m+1+const1
        rr=m*const2
c
        return
        END
c ******************************************************
       SUBROUTINE GETDM(D,AL,RI,RM,CM,NSEG,iw,ir,isn,ish,idsp,
     c   JA,IA,IRR,ICC,IC,NEQ,MAXIEL,IEL,fac,ar,spta,irepl,
     c   inode,inmax,ndau,IPRT,ID,ISOMAD,IND)
       IMPLICIT real*8 (A-H,O-Z)
       DIMENSION D(Nseg),AL(Nseg),fac(nseg),ar(nseg),spta(nseg),
     C   RI(NSEG),RM(NSEG),CM(NSEG),irepl(nseg),
     c   JA(MAXIEL),IA(NEQ+1),IRR(NEQ),ICC(NEQ),IC(NEQ),
     c   inode(inmax),NDAU(nseg),IPRT(nseg),ID(4,nseg),ISOMAD(20)
       character*50 filen,filen2
       character*10 label
c
       idI=34
       ico=35
       PI=datan(1.0d0)*4.d0
       write(iw,*) ' enter filename for dimensions xxxxxx.dat '
       read(ir,10) filen
       write(iw,10) filen
       write(iw,*) ' enter filename for connectivity xxxxxx.dat '
       read(ir,10) filen2
       write(iw,10) filen2
 10    format(a50)
       open(unit=idI,file=filen,status='old')
       open(unit=ico,file=filen2,status='old')
c
       inode(1)=neq
       k=2
       write(iw,*) ' enter diam,len mult factors to get to um'
       write(iw,*) ' i.e., if d=12 means d=1.2 um enter .1'
       read(ir,*) fdi,fln
       do 100 i=1,nseg
          read(idI,*) ibr,isg,di,alen,iseg,inod,iloc,ire
c         read(idI,*) label,isg,iseg,alen,di
c9        format(10a,i4,i5,2f10.5)
c         irepl(iseg)=1
c         iloc=1
c         inod=0
          irepl(iseg)=ire
          spn=0.d0
          spa=0.d0
          spta(i)=0.d0
          if(iloc.eq.1) then
c             spn=242.01195d+2
             spn=1.34d+4
c     spn= spine # per cm   spa= average spine area
             spa=1.7794996d-8
          endif
          if(iloc.eq.2) then
c             spn=211.99237d+2
             spn=1.22d+4
             spa=1.7376424d-8
          endif
         if(iloc.eq.3) then
c             spn=189.92433d+2
             spn=1.16d+4
             spa=1.6372787d-8
          endif
          d(iseg)=di*1.0d-4*fdi
          al(iseg)=alen*1.0d-4*fln
          ar(iseg)=pi*d(iseg)*al(iseg)
          spta(iseg)=spa*spn*al(iseg)
          if(inod.eq.0.or.inod.eq.2) then
             inode(k)=iseg
             k=k+1
             if(k.gt.inmax) write(iw,*) ' inode max exceeded '
          endif
 100   continue
       IND=K-1
c
       DO 7000 I=1,NEQ
           IRR(I)=I
           ICC(I)=I
           IC(I)=I
 7000  CONTINUE
C
C                                         get connectivity
C
c       WRITE (Iw,*) ' ENTER CONNECTIVITY. '
c       WRITE (Iw,*) ' NOTE ONLY ONE SEGMENT CAN HAVE MORE THAN 3 '
c       WRITE (Iw,*) ' DAUGHTERS.  ENTER THIS ONE FIRST!!! ENTER '
c       WRITE (Iw,*) ' ONLY SEG#, PARENT#, AND #DAUGHTERS FOR THIS ONE '
c      READ (Ico,*) ISOMA,IPRT(ISOMA),NDAU(ISOMA)
C
C       WRITE (IU2,*) ' ENTER DAUGHTER SEGMENTS 5 PER LINE TIL DONE '
c      L=1
c      ILND=NDAU(ISOMA)/5
c      IREM=NDAU(ISOMA)-ILND
C
c      DO 150 I=1,ILND
c         READ (Ico,*) (ISOMAD(J),J=L,L+4)
c         L=L+5
c150   CONTINUE
C
c      IF (IREM.GT.0) READ (Ico,*) (ISOMAD(J),J=L,NDAU(ISOMA))
C
c      WRITE (Iw,*) ' NOW ENTER CONNECTIVITY FOR OTHER SEGMENTS '
c      WRITE (Iw,*) ' ENTER SEG#, PARENT#, #DAU, 1..#DAU DAUGHTERS '
C
c      DO 200 I=1,NSEG-1
       DO 200 I=1,NSEG
c         READ (Ico,*) label,isg,IS,IPRT(IS),NDAU(IS),
c    c       (ID(J,IS),J=1,NDAU(IS))
          READ (Ico,*) IS,IPRT(IS),NDAU(IS),(ID(J,IS),J=1,NDAU(IS))
          WRITE (iw,*) IS,IPRT(IS),NDAU(IS),
     c       (ID(J,IS),J=1,NDAU(IS))
c 8       format(10a,i5,i5,i5,i5,i5,i6)
 200   CONTINUE
c
       write(iw,*) ' enter default rm,ri, and cm values '
       read(ir,*) drm,dri,dcm
       do 260 i=1,nseg
          fac(i)=1.d0
          RM(I)=drm  
          RI(I)=dri  
          CM(I)=dcm  
 260   continue
c
       write(iw,*) ' enter first spine neck & head seg numbers '
       read(ir,*) isn,ish
       idsp=ish-isn
       ilasts=isn-1
       do 250 i=1,ilasts
          jj=id(2,i)
          if(jj.ge.isn.and.jj.lt.ish) then
             sta=spta(i)-irepl(jj)*(ar(jj)+ar(jj+idsp))
             if(sta.lt.0.d0) then
                write(iw,*) ' irepl too big ',jj,irepl(jj)
c               sta=spta(i)
                sta=0.d0
             endif
c            sta=0.d0
          else
             sta=spta(i)
          endif
          fac(i)=(sta+ar(i))/ar(i)
          RM(I)=drm/fac(i)
          CM(I)=dcm*fac(i)
 250   continue
c
c      do 290 i=1,nseg
c         write(iw,291) i,ar(i),spta(i),rm(i),cm(i),irepl(i)
c291      format(i4,2d12.4,2f8.4,i4)
c290   continue
c
C
C        form matrices--node, ia, ir   and get r and el
C
       M=1
       DO 300 I=1,NSEG
          IF(M+3.GT.MAXIEL) WRITE(IW,*) ' MAXIEL EXCEEDED '
          IA(I)=M
          JA(M)=IPRT(I)
          JA(M+1)=I
c         IF(I.EQ.ISOMA) THEN
c            DO 304 J=1,NDAU(I)
c               JA(M+1+J)=ISOMAD(J)
c304         CONTINUE
c         ELSE
             DO 305 J=1,NDAU(I)
                JA(M+1+J)=ID(J,I)
 305         CONTINUE
c         ENDIF
          M=M+2+NDAU(I)
 300   CONTINUE
       IA(NEQ)=M
       IA(NEQ+1)=M+2
       JA(M)=NEQ
       JA(M+1)=1
       IEL=M+1
       WRITE(iw,*) ' iel= ',iel,' PREORDER 1=NO '
       READ(ir,11) IPRE
       write(iw,11) IPRE
 11    FORMAT(I2)
       IF(IPRE.EQ.1) GOTO 9999
C
C   PREORDER ROWS
C
        DO 205 I=1,NEQ
 205       ICC(I)=0
        DO 210 K=1,NEQ
           KDEG=IA(K+1) - IA(K)
           IF (KDEG .EQ. 0) KDEG = KDEG + 1
           IC(K) = ICC(KDEG)
           ICC(KDEG) = K
 210       CONTINUE
        I=0
        DO 230 J= 1,NEQ
           IF (ICC(J) .EQ. 0) GOTO 230
           K = ICC(J)
 220       I = I+1
           IRR(I) = K
           K = IC(K)
           IF (K .GT. 0) GOTO 220
 230       CONTINUE
        DO 240 I=1,NEQ
           ICC(I)=I
           IC(I) = I
 240       CONTINUE
C
 9999  CONTINUE
       close(unit=idI)
       close(unit=ico)
       RETURN
       END
c  ************************************
       SUBROUTINE CHGSPC(D,AL,RI,RM,CM,ar,fac,spta,irepl,id,
     c   NSEG,iw,ir,PI,isn,ish,idsp)
       IMPLICIT real*8(A-H,O-Z)
       DIMENSION D(NSEG),AL(NSEG),RI(NSEG),RM(NSEG),CM(NSEG),
     c   fac(nseg),ar(nseg),spta(nseg),irepl(nseg),id(4,nseg)
 122   write(iw,12)
  12   format(' change length of any section 1=yes')
       read(ir,11) nn
  11   format(i2)
       if(nn.ne.1) goto 222
       write(iw,15)
  15   format(' enter START, END SEC, AND new length in um ')
       read(ir,*) NN1,nn,alnew
       DO 100 I=NN1,NN
          al(I)=alnew*1.0d-4
          ar(i)=pi*d(i)*al(i)
          if(i.lt.isn) then
             spn=0.d0
             spa=0.d0
             if(iloc.eq.1) then
                spn=242.01195d+2
                spa=1.7794996d-8
             endif
             if(iloc.eq.2) then
                spn=211.99237d+2
                spa=1.7376424d-8
             endif
             if(iloc.eq.3) then
                spn=189.92433d+2
                spa=1.6372787d-8
             endif
             spta(i)=spa*spn*al(iseg)
             rm(i)=rm(i)*fac(i)
             cm(i)=cm(i)/fac(i)
             fac(i)=spta(i)/ar(i)+1.d0
             rm(i)=rm(i)/fac(i)
             cm(i)=cm(i)*fac(i)
          endif
 100   CONTINUE
       goto 122
c
 222   CONTINUE
 110   write(iw,22)
  22   format(' change diam of any section 1=yes')
       read(ir,11) nn
       if(nn.ne.1) goto 14
       write(iw,25)
  25   format(' enter START,END section and new diam in um ')
       read(ir,*) NN1,nn,dnew
       DO 101 I=NN1,NN
          d(I)=dnew*1.0d-4
          ar(i)=pi*d(i)*al(i)
          if(i.lt.isn) then
             rm(i)=rm(i)*fac(i)
             cm(i)=cm(i)/fac(i)
             fac(i)=spta(i)/ar(i)+1.d0
             rm(i)=rm(i)/fac(i)
             cm(i)=cm(i)*fac(i)
          endif
  101  CONTINUE
       goto 110
C
 14    continue
C
 94    CONTINUE
       WRITE(iw,45)
 45    FORMAT('  CHANGE RM,RI,CM DEFAULTS 1=YES')
       READ(ir,11) NN
       WRITE(iW,11) NN
       IF (NN.NE.1) GOTO 188
       WRITE(iw,770)
 770   FORMAT(' ENTER NEW DEFAULT RM,RI,CM ')
       READ(ir,*) RMM,RII,CMM
       WRITE(iW,*) RMM,RII,CMM
       DO 18 I=1,nseg
          RI(I)=RII
          RM(I)=RMM/fac(i)
          CM(I)=CMM*fac(i)
 18    CONTINUE
c
 188   WRITE (iw,19)
 19    FORMAT(' DO YOU WANT TO CHANGE RI 1=YES')
       READ(ir,11) NN
       IF (NN.NE.1) GOTO 629
       WRITE(iw,321)
 321   FORMAT(' ENTER START, END SECTION AND NEW VALUE I2,F12.5')
       READ(ir,*) NN,NN1,RII
       do 20 i=NN,NN1
          RI(i)=RII
 20    continue
       GOTO 188
c
 629   CONTINUE
 488   WRITE (iw,49)
 49    FORMAT(' DO YOU WANT TO CHANGE RM 1=YES')
       READ(ir,11) NN
       IF (NN.NE.1) GOTO 449
       WRITE(iw,321)
       READ(ir,*) NN,NN1,RMM
       do 40 i=NN,NN1
          RM(i)=RMM/fac(i)
 40    continue
       GOTO 488
c
 449   CONTINUE
 588   WRITE (iw,59)
 59    FORMAT(' DO YOU WANT TO CHANGE CM 1=YES')
       READ(ir,11) NN
       IF (NN.NE.1) GOTO 549
       WRITE(iw,321)
       READ(ir,*) NN,NN1,CMM
       do 50 i=NN,NN1
          CM(i)=CMM*fac(i)
 50    continue
       GOTO 588
c
 549   CONTINUE
c
 1020  write(iw,*) ' change replication value 1=yes '
       read(ir,11) nn
       if(nn.ne.1) goto 1900
       write(iw,*) ' enter start segment end seg and repl value '
       read(ir,*) ist,iend,ire
       do 1100 j=ist,iend
          irepl(j)=ire
 1100  continue
       goto 1020
 1900  continue
c
       do 250 i=1,isn-1
          jj=id(2,i)
          if(jj.ge.isn.and.jj.lt.ish) then
             sta=spta(i)-irepl(jj)*(ar(jj)+ar(jj+idsp))
             if(sta.lt.0.d0) then
                write(iw,*) ' irepl too big ',jj,irepl(jj)
c               sta=spta(i)
                sta=0.d0
             endif
c            sta=0.d0
          else
             sta=spta(i)
          endif
          rm(i)=rm(i)*fac(i)
          cm(i)=cm(i)/fac(i)
          fac(i)=(sta+ar(i))/ar(i)
          RM(I)=rm(i)/fac(i)
          CM(I)=cm(i)*fac(i)
 250   continue
c
       kk=360
       if(kk.gt.nseg) kk=nseg
c      do 290 i=1,kk 
c         write(iw,291) i,ar(i),spta(i),rm(i),cm(i),irepl(i)
c291      format(i4,2d12.4,2f8.4,i4)
c290   continue
c
       totmar=0.d0
       do 2000 i=1,nseg
          res=4.d0*al(i)*1.0d-4/(pi*d(i)*d(i))
          rmarea=pi*d(i)*al(i)*1.0d+8*dfloat(irepl(i))
          if(i.ge.ish) rmarea=rmarea*dfloat(irepl(i-idsp))
          totmar=totmar+rmarea
c          write(iw,2001) i,res,rmarea
 2001     format(' seg ',i4,' ri ',d18.6,' rm ',d18.6)
 2000     continue
       write(iw,2002) totmar
 2002  format(' totmar ',d18.6)
c
       RETURN
       END
C ***************************************************
       SUBROUTINE PRTEL(AL,ALAMDA,RM,NSEG,IPRT,INODE,IND,IW,IR)
       IMPLICIT real*8 (A-H,O-Z)
       DIMENSION AL(NSEG),ALAMDA(NSEG),EL(990),TL(990),TEL(990),
     C   AAL(990),TTL(990),RM(NSEG),IPRT(NSEG),INODE(IND)
c
       DO 100 I=1,NSEG
          EL(I)=AL(I)/ALAMDA(I)
 100   CONTINUE
       TL(1)=AL(1)
       TEL(1)=EL(1)
       DO 200 I=2,NSEG
          IF(IPRT(I).EQ.1) THEN
             TL(I)=AL(I)
             TEL(I)=EL(I)
          ELSE
             TL(I)=AL(I)+TL(IPRT(I))
             TEL(I)=EL(I)+TEL(IPRT(I))
          ENDIF
 200   CONTINUE
c
       DO 9033 I=1,NSEG
           AAL(I)=AL(I)*1.0d+4
           TTL(I)=TL(I)*1.0d+4
 9033  CONTINUE
c
       write(iw,*) ' to node points only. 1=yes '
       read(ir,11) nn
 11    format(i2)
       write(iw,*) ' write to file 1=yes '
       read(ir,11) nnn
       if(nnn.eq.1) open(file='ELDAT',unit=63,status='new')
       if(nn.eq.1) then
          DO 296 j=2,IND
             i=inode(j)
             RMM=RM(I)*1.0d+3
             WRITE(iw,9001) I,AAL(I),TTL(I),EL(I),TEL(I),RMM
             if(nnn.eq.1) WRITE(63,9002) I,AAL(I),
     c         TTL(I),EL(I),TEL(I),RMM
 9001        FORMAT(' SEG ',I4,' LEN ',F9.4,' TLEN ',F9.4,' EL ',
     C         F8.5,' TEL ',F8.5,' RM ',F10.4)
 9002        format(i4,2f9.2,2f8.5,f10.3)
 296      CONTINUE
       else
          DO 299 I=1,nseg
             RMM=RM(I)*1.0d+3
              WRITE(iw,9001) I,AAL(I),TTL(I),EL(I),TEL(I),RMM
             if(nnn.eq.1) WRITE(63,9002) I,AAL(I),
     c         TTL(I),EL(I),TEL(I),RMM
 299      CONTINUE
       endif
       if(nnn.eq.1) close(unit=63)
c
       RETURN
       ENd
C *****************************************************
       SUBROUTINE CURIP(H,h2,h3,CUR,cur2,cur3,
     c   NVC,NIP,TEXP,NCURLO,CSTR,CEND,CEND2,iw,ir)
       IMPLICIT real*8(A-H,O-Z)
       WRITE(iw,77)
 77    FORMAT(' DO YOU WANT VOLTAGE CLAMP 1=YES')
       READ(ir,11) NVC
 11    FORMAT(I2)
       WRITE (iw,110)
 110   FORMAT(' ENTER cur(OR VOLTAGE IF VC)  ')
       READ (ir,*) cur,cur2,cur3
       WRITE (iW,*) cur,cur2,cur3
       IF(NVC.EQ.1) GOTO 1980
       WRITE(iw,901)
 901   FORMAT(' ENTER 1=DELTA FUN, 2=CONSTANT,3=A2TE-AT ')
       READ (ir,11) NIP
       WRITE(iW,11) NIP
       if(nip.eq.2) then
          write(iw,*) ' enter start, end times '
          read(ir,*) cstr,cend,cend2
       endif
       IF(NIP.NE.3) GOTO 339
       WRITE(iw,910)
 910   FORMAT(' ENTER alpha I.E. 10  ')
       READ(ir,*) TEXP
 339   WRITE (iw,130)
 130   FORMAT('  WHERE INPUT--ENTER SEGMENT # 0=SOMA')
       READ(ir,12) NCURLO
       WRITE(iW,12) NCURLO
 12    FORMAT(I4)
 138   H=.001d0*CUR
       h2=.001d0*cur2
       h3=.001d0*cur3
 1980  CONTINUE
       RETURN
       END
C *****************************************************
       SUBROUTINE CONDIS(NCOND,NSEG,NSC,igs,vrs,icstyp,NperMv,
     c   gbar,exm,exh,taum,tauh,aminf,ahinf,brp,cai,iw,ir)
       IMPLICIT real*8(A-H,O-Z)
       DIMENSION IGS(NSC,NSEG),VRS(nsc),icstyp(NSC),gbar(NSC,NSEG),
     C   exm(NSC),exh(NSC),taum(NSC,1500),tauh(NSC,1500),
     C   aminf(NSC,1500),ahinf(NSC,1500),cai(nseg)
       NCOND=1
C
       DO 100 I=1,NSC
          WRITE(IW,*) ' DO YOU WANT TO CHANGE CONDUCTANCE',I,
     C      ' ? 1=YES '
          READ(IR,11) NN
 11       format(i2)
          WRITE(IW,11) NN
          IF(NN.NE.1) GOTO 100
          write(iw,*) ' enter conductance type, 1=constant, '
          write(iw,*) ' 2=voltage dependent (HH), 3=Ca (HH volt dep)'
c         write(iw,*) ' 4=Ca dependent K, 5=A current 6=KC    '
          read(ir,11) icstyp(i)
          WRITE(IW,*) ' ENTER REVERSAL POTENTIAL '
          READ(IR,*) VRS(I)
          vrs(i)=vrs(i)-brp
          WRITE(IW,*) VRS(I)
c
          if(icstyp(i).eq.1) then
             exm(i)=0.d0
             exh(i)=0.d0
          endif
c
          if(icstyp(i).eq.2.or.icstyp(i).eq.3) then
             write(iw,*) ' enter temperature.  Default is 6.3 C '
             read(ir,*) temp
             q10=3.d0**((temp-6.3d0)/10.d0)
 31          write(iw,*) ' which ion? 1=Na, 2=K, 3=Ca '
             read(ir,*) iontyp
             if(iontyp.le.0.or.iontyp.ge.5) goto 31
             if(iontyp.eq.1) then
                write(iw,*)
     c            ' do you want 1) HH 2)Y and D 3) Traub or 4)other '
                write(iw,*) ' enter 1 or 2 or 3 or 4 '
                read(ir,*) idef
                if(idef.eq.1) then
                   exm(i)=3.d0
                   exh(i)=1.d0
c                  write(iw,32) exm(i),exh(i)
c32                format(' default Na m and h exponents are ',2f4.0)
c                  write(iw,*) '  v-d params [a+bV]/[c+exp((d+V)/f)] '
c                  write(iw,*) '  m: alpha a=-18.3 b=-0.3 c=-1.0 ',
c    c              ' d=61.0 f=-05.0 '
c                  write(iw,*) '  m: beta  a=9.900 b=0.3 c=-1.0 ',
c    c              ' d= 33.0 f=05.0 '
c                  write(iw,*) '  h: alpha a=0.23 b=0.0 c=0.0 ',
c    c              ' d=83.0 f=20.0 '
c                  write(iw,*) '  h: beta a=3.33 b=0.0 c=1.0 ',
c    c              ' d= 30.5 f=-10.0 '
                   cmaa=-4.5d0
                   cmab=-.1d0
                   cmac=-1.0d0
                   cmad=45.d0
                   cmaf=-10.d0
                   cmba=4.d0
                   cmbb=0.d0
                   cmbc=0.d0
                   cmbd=70.d0
                   cmbf=18.d0
                   chaa=0.07d0
                   chab=0.d0
                   chac=0.d0
                   chad=70.d0
                   chaf=20.d0
                   chba=1.d0
                   chbb=0.d0
                   chbc=1.d0
                   chbd=40.d0
                   chbf=-10.d0
                endif
                if(idef.eq.2) then
                   write(iw,*) ' enter Na shift '
                   read(ir,*) sh
                   exm(i)=3.d0
                   exh(i)=1.d0
c                  cmaa=-18.3d0
c                  cmaa=-15.3d0
                   cmab=-0.3d0
                   cmaa=cmab*(70.d0-9.d0-sh)
                   cmac=-1.0d0
                   cmad=61.0d0-sh
                   cmaf=-5.0d0
c                  cmba=9.900d0
c                  cmba=6.900d0
                   cmbb=0.3d0
                   cmba=cmbb*(70.d0-37.d0-sh)
                   cmbc=-1.0d0
                   cmbd= 33.d0-sh
                   cmbf=05.d0
                   chaa=0.23d0
                   chab=0.0d0
                   chac=0.0d0
                   chad=83.d0-sh
                   chaf=20.d0
                   chba=3.33d0
                   chbb=0.0d0
                   chbc=1.0d0
                   chbd= 30.5d0-sh
                   chbf=-10.0d0
                endif
                if(idef.eq.3) then
                   exm(i)=2.d0
                   exh(i)=1.d0
                   cmaa=-18.208d0
                   cmab=-.32d0
                   cmac=-1.d0
                   cmad=56.9d0
                   cmaf=-4.d0
                   cmba=8.372d0
                   cmbb=0.28d0
                   cmbc=-1.d0
                   cmbd=29.9d0
                   cmbf=5.d0
                   chaa=0.128d0
                   chab=0.d0
                   chac=0.d0
                   chad=53.d0
                   chaf=18.d0
                   chba=4.d0
                   chbb=0.d0
                   chbc=1.d0
                   chbd=30.d0
                   chbf=-5.d0
                endif
                if(idef.eq.4) then
                   write(iw,*) ' enter m and h exponents '
                   read(ir,*) exm(i), exh(i)
                   write(iw,*) ' enter m alpha parameters: a,b,c,d,f '
                   read(ir,*) cmaa,cmab,cmac,cmad,cmaf
                   write(iw,*) ' enter m beta  parameters: a,b,c,d,f '
                   read(ir,*) cmba,cmbb,cmbc,cmbd,cmbf
                   write(iw,*) ' enter h alpha parameters: a,b,c,d,f '
                   read(ir,*) chaa,chab,chac,chad,chaf
                   write(iw,*) ' enter h beta  parameters: a,b,c,d,f '
                   read(ir,*) chba,chbb,chbc,chbd,chbf
                endif
C
                NperMv=10
                Ntab=150*NperMv
                V=-100.d0
                dV=1.d0/dfloat(NperMv)
                do 200 k=1,Ntab
                   alfm=(cmaa+cmab*V)/(cmac+dexp((cmad+V)/cmaf))*q10 
                   betm=(cmba+cmbb*V)/(cmbc+dexp((cmbd+V)/cmbf))*q10 
                   alfh=(chaa+chab*V)/(chac+dexp((chad+V)/chaf))*q10 
                   beth=(chba+chbb*V)/(chbc+dexp((chbd+V)/chbf))*q10 
                   taum(i,k)=1.d0/(alfm+betm)
                   tauh(i,k)=1.d0/(alfh+beth)
                   aminf(i,k)=alfm*taum(i,k)
                   ahinf(i,k)=alfh*tauh(i,k)
                   V=V+dV
 200            continue
c
C               compute NA tables and initial cond
C
             endif
c
             if(iontyp.eq.2) then
                write(iw,*)
     c              ' do you want 1) HH 2)Y and D 3) Traub or 4)other '
                write(iw,*) ' enter 1 or 2 or 3 or 4 '
                read(ir,*) idef
c               write(iw,*) ' defaults K values are n**4  '
                exm(i)=4.d0
                exh(i)=0.d0
                if(idef.eq.1) then
c                  write(iw,*) '  v-d params [a+bV]/[c+exp((d+V)/f)] '
c                  write(iw,*) '  n: alpha a=-2.73 b=-0.07 c=-1.0 ',
c    c               ' d= 39.0 f=-06.0 '
c                  write(iw,*) '  n: beta  a=.264 b=0.0 c=0.0 ', 
c    c               ' d=64.0 f=40.0 '
                   cnaa=-0.6d0
                   cnab=-0.01d0
                   cnac=-1.d0
                   cnad=60.d0
                   cnaf=-10.d0
                   cnba=.125d0
                   cnbb=0.d0
                   cnbc=0.d0
                   cnbd=70.d0
                   cnbf=80.d0
                endif
                if(idef.eq.2) then
                   exm(i)=4.d0
                   write(iw,*) ' enter K shift '
                   read(ir,*) sh
c                  cnaa=-2.73d0
c                  cnaa=-2.03d0
                   cnab=-0.07d0
                   cnaa=cnab*(70.d0-31.d0-sh)
                   cnac=-1.0d0
                   cnad= 39.d0-sh
                   cnaf=-06.0d0
                   cnba=0.264d0
                   cnbb=0.d0
                   cnbc=0.d0
                   cnbd=64.d0-sh
                   cnbf=40.d0
                endif
                if(idef.eq.3) then
                   exm(i)=1.d0
                   cnaa=-0.5584
                   cnab=-0.016
                   cnac=-1.d0
                   cnad=34.9d0
                   cnaf=-5.d0
                   cnba=0.25d0
                   cnbb=0.d0
                   cnbc=0.d0
                   cnbd=50.d0
                   cnbf=40.d0
                endif
                if(idef.eq.4) then
                   write(iw,*) ' enter n exponent '
                   read(ir,*) exm(i)
                   write(iw,*) ' enter n alpha parameters: a,b,c,d,f '
                   read(ir,*) cnaa,cnab,cnac,cnad,cnaf
                   write(iw,*) ' enter n beta  parameters: a,b,c,d,f '
                   read(ir,*) cnba,cnbb,cnbc,cnbd,cnbf
                endif
c
C
                V=-100.d0
                dV=1.d0/dfloat(NperMv)
                do 300 k=1,Ntab
                   alfn=(cnaa+cnab*V)/(cnac+dexp((cnad+V)/cnaf))*q10 
                   betn=(cnba+cnbb*V)/(cnbc+dexp((cnbd+V)/cnbf))*q10 
                   taum(i,k)=1.d0/(alfn+betn)
                   tauh(i,k)=1.d0
                   aminf(i,k)=alfn*taum(i,k)
                   ahinf(i,k)=1.d0 
                   V=V+dV
 300            continue
c
C               compute K tables and initial cond
C
             endif
c
             if(iontyp.eq.3) then
                write(iw,*) ' defaults Ca values are s**2 r**1 '
                write(iw,*)
     c            ' do you want 2) Y and D 3) Traub or 4) other '
                write(iw,*) ' enter 2 or 3 or 4 '
                read(ir,*) idef
                if(idef.eq.2) then
c                  write(iw,*) '  v-d params [a+bV]/[c+exp((d+V)/f)] '
c                  write(iw,*) '  m: alpha a=-1.4105 b=-.031 c=-1.0 ',
c    c              ' d= 45.5 f=-3.00 '
c                  write(iw,*) '  m: beta  a=0.93 b=.031 c=-1.0 ',
c    c              ' d= 30.0 f=3.00 '
c                  write(iw,*) '  h: alpha a=3d-6 b=0.0 c=0.0 ',
c    c              ' d= 58.5 f=17.5 '
c                  write(iw,*) '  h: beta a=3d-4 b=0.0 c=1.0 ',
c    c              ' d= 50.5 f=-14.0 '
                   exm(i)=2.d0
                   exh(i)=1.d0
                   write(iw,*) ' enter Ca voltage shift '
                   read(ir,*) sh
c                  cmaa=-1.4105d0
                   cmab=-.031d0
                   cmaa=cmab*(70.d0-24.5d0-sh)
                   cmac=-1.d0
                   cmad= 45.5d0-sh
                   cmaf=-3.0d0  
c                  cmba=0.93d0
                   cmbb=0.031d0
                   cmba=cmbb*(70.d0-40.d0-sh)
                   cmbc=-1.d0
                   cmbd= 30.d0-sh
                   cmbf=3.0d0
                   chaa=3d-6  
                   chab=0.d0
                   chac=0.d0
                   chad= 58.5d0-sh
                   chaf=17.5d0
                   chba=3d-4  
                   chbb=0.d0
                   chbc=1.d0
                   chbd= 50.5d0-sh
                   chbf=-14.d0
                endif
                if(idef.eq.3) then
                   exm(i)=2.d0
                   exh(i)=1.d0
                   cmaa=1.6d0
                   cmab=0.d0
                   cmac=1.d0
                   cmad=5.d0
		   cmaf=13.89
                   cmba=.378d0
                   cmbb=.02d0
                   cmbc=-1.d0
                   cmbd=19.9d0
                   cmbf=5.d0
                   chaa=.005d0
                   chab=0.d0
                   chac=0.d0
                   chad=70.d0
                   chaf=20.d0
                endif
                if(idef.eq.4) then
                   write(iw,*) ' enter m and h exponents '
                   read(ir,*) exm(i), exh(i)
                   write(iw,*) ' enter m alpha parameters: a,b,c,d,f '
                   read(ir,*) cmaa,cmab,cmac,cmad,cmaf
                   write(iw,*) ' enter m beta  parameters: a,b,c,d,f '
                   read(ir,*) cmba,cmbb,cmbc,cmbd,cmbf
                   write(iw,*) ' enter h alpha parameters: a,b,c,d,f '
                   read(ir,*) chaa,chab,chac,chad,chaf
                   write(iw,*) ' enter h beta  parameters: a,b,c,d,f '
                   read(ir,*) chba,chbb,chbc,chbd,chbf
                endif
C
                V=-100.d0
                dV=1.d0/dfloat(NperMv)
                do 400 k=1,Ntab
                   alfm=(cmaa+cmab*V)/(cmac+dexp((cmad+V)/cmaf))*q10 
                   betm=(cmba+cmbb*V)/(cmbc+dexp((cmbd+V)/cmbf))*q10 
                   if(idef.eq.3) then
                      if(v.lt.brp) then
                         alfh=0.005d0
                         beth=0.d0
                      else
                         alfh=(chaa+chab*V)/(chac+dexp((chad+V)/chaf))*q10 
                         beth=0.005d0-alfh
                      endif 
                   else
                      alfh=(chaa+chab*V)/(chac+dexp((chad+V)/chaf))*q10 
                      beth=(chba+chbb*V)/(chbc+dexp((chbd+V)/chbf))*q10 
                   endif 
            
                   taum(i,k)=1.d0/(alfm+betm)
                   tauh(i,k)=1.d0/(alfh+beth)
                   aminf(i,k)=alfm*taum(i,k)
                   ahinf(i,k)=alfh*tauh(i,k)
                   V=V+dV
 400            continue
c
C               compute CA tables and initial cond
C
             endif
c
          endif
c
          if(icstyp(i).eq.4) then
             write(iw,*) ' defaults K values are q**2  '
c            write(iw,*) '  v-d params [a/[exp((c+d*log(Ca))/f)] '
c            write(iw,*) '  q: alpha a=0.0041 b=0.0 c=4.48 ',
c    c         ' d=10.0 f=-4.50 '
c            write(iw,*) '  q: beta  a=0.01 b=0.0 c=36.4 ', 
c    c         ' d=10.0 f=35.0 '
             write(iw,*) ' do you want 2) Y and D 3) Traub or 4) other '
             write(iw,*) ' enter 2 or 3 or 4 '
             read(ir,*) idef
             if(idef.eq.2) then
                exm(i)=2.d0
                exh(i)=0.d0
                write(iw,*) ' enter ca shift '
                read(ir,*) sh
                cnaa=0.0041d0
                cnab=0.0d0
                cnac=4.48d0
                cnad=10.d0+sh
                cnaf=-4.50d0
                cnba=0.01d0
                cnbb=0.d0
                cnbc=36.4d0
                cnbd=10.d0+sh
                cnbf=35.d0
             endif
             if(idef.eq.3) then
                exm(i)=1.d0
                exh(i)=0.d0
             endif
             if(idef.eq.4) then
                write(iw,*) ' enter q exponent '
                read(ir,*) exm(i)
                write(iw,*) ' enter q alpha parameters: a,b,c,d,f '
                read(ir,*) cnaa,cnab,cnac,cnad,cnaf
                write(iw,*) ' enter q beta  parameters: a,b,c,d,f '
                read(ir,*) cnba,cnbb,cnbc,cnbd,cnbf
             endif
c
C
             ca=0.01d0
             cal=dlog10(ca)
             dca=4.d0/dfloat(Ntab)
             do 500 k=1,Ntab
c               alfn=cnaa/(dexp((cnac+cnad*cal)/cnaf))*q10 
c               betn=cnba/(dexp((cnbc+cnbd*cal)/cnbf))*q10 
                if(idef.eq.2) then
                   alfn=cnaa/(dexp((cnac+cnad*cal)/cnaf))
                   betn=cnba/(dexp((cnbc+cnbd*cal)/cnbf)) 
                endif
                if(idef.eq.3) then
                   alfn=0.2d-4*(10.0**cal)
                   if(alfn.gt.0.01) alfn=0.01
                   betn=0.001d0
                endif
                taum(i,k)=1.d0/(alfn+betn)
                tauh(i,k)=1.d0
                aminf(i,k)=alfn*taum(i,k)
                ahinf(i,k)=1.d0 
                cal=cal+dca
 500         continue
c
             do 600 k=1,nseg
                cai(k)=0.07d0
 600         continue
             write(iw,*) ' initial Ca = .07 uM, OK? 1=no '
             read(ir,11) nn
             if(nn.eq.1) then
 611            write(iw,*) ' enter new Ca, start, end segments '
                read(ir,*) canew,ist,iend
                do 620 k=ist,iend
                   cai(k)=canew
 620            continue
                write(iw,*) ' any more 1=YES '
                read(ir,11) nn1
                if(nn1.eq.1) goto 611
             endif
c
C               compute K tables and initial cond
C
          endif
c
             if(icstyp(j).eq.5) then
	        write(iw,*) ' enter 3 for Traub, 4 for other '
                read(ir,*) idef
                if(idef.eq.3) then
                   exm(i)=1.d0
                   exh(i)=1.d0
                   cmaa=1.138d0
                   cmab=-.02d0
                   cmac=-1.d0
                   cmad=56.9d0
                   cmaf=-10.d0
                   cmba=0.52325d0
                   cmbb=.0175d0
                   cmbc=-1.d0
                   cmbd=29.9d0
                   cmbf=10.d0
                   chaa=0.0016d0
                   chab=0.d0
                   chac=0.d0
                   chad=83.d0
                   chaf=18.d0
                   chba=0.05d0
                   chbb=0.d0
                   chbc=1.d0
                   chbd=59.9d0
                   chbf=-5.d0
                endif
                if(idef.eq.4) then
                   write(iw,*) ' enter m and h exponents '
                   read(ir,*) exm(i), exh(i)
                   write(iw,*) ' enter m alpha parameters: a,b,c,d,f '
                   read(ir,*) cmaa,cmab,cmac,cmad,cmaf
                   write(iw,*) ' enter m beta  parameters: a,b,c,d,f '
                   read(ir,*) cmba,cmbb,cmbc,cmbd,cmbf
                   write(iw,*) ' enter h alpha parameters: a,b,c,d,f '
                   read(ir,*) chaa,chab,chac,chad,chaf
                   write(iw,*) ' enter h beta  parameters: a,b,c,d,f '
                   read(ir,*) chba,chbb,chbc,chbd,chbf
                endif
C
                NperMv=10
                Ntab=150*NperMv
                V=-100.d0
                dV=1.d0/dfloat(NperMv)
                do 700 k=1,Ntab
                   alfm=(cmaa+cmab*V)/(cmac+dexp((cmad+V)/cmaf))*q10 
                   betm=(cmba+cmbb*V)/(cmbc+dexp((cmbd+V)/cmbf))*q10 
                   alfh=(chaa+chab*V)/(chac+dexp((chad+V)/chaf))*q10 
                   beth=(chba+chbb*V)/(chbc+dexp((chbd+V)/chbf))*q10 
                   taum(i,k)=1.d0/(alfm+betm)
                   tauh(i,k)=1.d0/(alfh+beth)
                   aminf(i,k)=alfm*taum(i,k)
                   ahinf(i,k)=alfh*tauh(i,k)
                   V=V+dV
 700            continue
c
         endif
c
         if(icstyp(j).eq.6) then
C
             exm(i)=1.d0
             exh(i)=0.d0
c
                NperMv=10
                Ntab=150*NperMv
                V=-100.d0
                dV=1.d0/dfloat(NperMv)
                do 800 k=1,Ntab
                   if(v.lt.-20.d0) then
                      alfm=dexp((v+60.d0)/11.d0 - (v+64.5d0)/27)/18.975
                      betm=2.d0*dexp(-(v+64.5)/27)-alfm
                   else
                      alfm=2.d0*dexp(-(v+64.5)/27)
                      betm=0.d0
                   endif
                   taum(i,k)=1.d0/(alfm+betm)
                   tauh(i,k)=1.d0/(alfh+beth)
                   aminf(i,k)=alfm*taum(i,k)
                   ahinf(i,k)=alfh*tauh(i,k)
                   V=V+dV
 800            continue
c
         endif
c
c
 3        WRITE(IW,*) ' ENTER START SEG, END SEG '
          read(ir,*) ist,iend
          WRITE(IW,*) IST,IEND
          write(iw,*) ' enter cond in mS/cm2 '
          read(ir,*) gb
          DO 110 J=IST,IEND
             gbar(i,j)=gb
             iGs(I,J) = 1
 110      CONTINUE
          WRITE(IW,*) ' ANY MORE? 1=YES '
          READ(IR,11) NN
          WRITE(IW,11) NN
          IF(NN.EQ.1) GOTO 3
c
 100   CONTINUE
C
       return
       end
C *****************************************************
       SUBROUTINE CONDIP(NCOND,IGPC,VRP,NSEG,NC,tp,brp,
     c   ITP,ITYP,ALF,BE,A1,A2,A1A,A2A,A1BK,A2BK,SA2B2,X0,X0MAX,
     c   SCC,tabpc,iw,ir)
       IMPLICIT real*8(A-H,O-Z)
       DIMENSION IGPC(NC,NSEG),VRP(nc),tp(NC,300),
     C   ITP(NC),ITYP(NC),ALF(NC),BE(NC),A1(NC),A2(NC),A1A(NC),
     C   A2A(NC),A1BK(NC),A2BK(NC),SA2B2(NC),X0(NC),X0MAX(NC),
     C   SCC(NC),tabpc(nc,1500000),x(18)
       character*50 file1
       NCOND=1
C
       DO 200 I=1,NC
          WRITE(IW,*) ' DO YOU WANT TO CHANGE CONDUCTANCE',I,
     C      ' ? 1=YES '
          READ(IR,11) NN
 11       format(i2)
          WRITE(IW,11) NN
          IF(NN.NE.1) GOTO 200
          WRITE(IW,*) ' ENTER REVERSAL POTENTIAL, TYPE, #INPUTS '
          READ(IR,*) VRP(I),ITYP(I),ITP(I)
          vrp(i)=vrp(i)-brp
          WRITE(IW,*) VRP(I),ITYP(I),ITP(I)
          if(ityp(i).eq.1.or.ityp(i).eq.2.or.ityp(i).eq.3) then
             k=1
             write(iw,*) ' enter lookup file name'
             read(ir,10) file1
 10          format(a50)
             write(iw,*) ' enter number of columns, column to use '
             read(ir,*) ncol,iuse
             open(file=file1,unit=61,status='old')
 66          read(61,*,end=77) (x(nlkup),nlkup=1,ncol)
             tabpc(i,k)=x(iuse)
             k=k+1
             goto 66
 77          continue
             close(unit=61)
          else if(ityp(i).eq.4.or.ityp(i).eq.14) then
             write(iw,*) ' enter the conductance value for constant g '
             read(ir,*) x0(I) 
          else 
             WRITE(IW,*) ' ENTER ALPHA,BETA,AKM1 '
             READ(IR,*) ALF(I),BE(I),AKM1
             WRITE(IW,*) ALF(I),BE(I),AKM1
             WRITE(IW,*) ' ENTER X0,X0MAX '
             READ(IR,*) X0(I),X0MAX(I)
             WRITE(IW,*) X0(I),X0MAX(I)
c
             aa=0.5*(alf(i)+be(i)+akm1)
             bb=alf(i)*akm1
             sa2b=dsqrt(aa*aa-bb)
             a1(i)=-aa+sa2b
             a2(i)=-aa-sa2b
             a1a(i)=a1(i)+alf(i)
             a2a(i)=a2(i)+alf(i)
             a1bk(i)=a1(i)+be(i)+akm1
             a2bk(i)=a2(i)+be(i)+akm1
             sa2b2(i)=2.d0*sa2b
         endif
c
c
          WRITE(IW,*) ' ENTER ',ITP(I),' TIME POINTS FOR INPUTS '
          READ(IR,*) (TP(I,J),J=1,ITP(I))
c         write(Iw,*) (TP(I,J),J=1,ITP(I))
          WRITE(IW,*) ' ENTER SINGLE CHANNEL COND IN pS '
          READ(IR,*) SCC(I)
          WRITE(IW,*) SCC(I)
c         SCAC(I)=YFAC(I)*0.472d-6
 4        WRITE(IW,*) ' ENTER START SEG, END SEG '
          READ(IR,*) IST,IEND
          WRITE(IW,*) IST,IEND
          DO 210 J=IST,IEND
             IGPC(I,J) = 1
 210      CONTINUE
          WRITE(IW,*) ' ANY MORE? 1=YES '
          READ(IR,11) NN
          WRITE(IW,11) NN
          IF(NN.EQ.1) GOTO 4
 200   CONTINUE
c      do 300 i=1,nc-1
c         if(ityp(i).eq.1.and.ityp(i+1).ne.2)
c    c     write(iw,*) ' bad conductance types, type 2 must follow',
c    c      ' type 1 '
c300   continue
C
 9999  CONTINUE
       RETURN
       END
C  ***************************************************
        SUBROUTINE getrin(ir,iw,nri,stt,ett,fr,nc,igpc,kks,nsyn,irc,
     c                     nseg)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION  irc(nc,nri),stt(nri),ett(nri),fr(nri),kks(nri,nseg),
     c        igpc(NC,NSEG),nsyn(nri)
c
        WRITE(IW,*) 'HOW MANY DIFFERENT RANDOM INPUTS?'
        READ(IR,*) NRI
        WRITE(IW,*) NRI
        DO 701 K=1,NRI
         KS=0
         WRITE(iw,*) 'RANDOM INPUT NUMBER',K
         WRITE(IW,*) 'ENTER START TIME, END TIME, FREQUENCY'
         READ(IR,*) STT(K),ETT(K),FR(K)
         WRITE(iw,*) STT(K),ETT(K),FR(K)
  44     continue
         WRITE(iw,*) ' ENTER START SEG, END SEG '
         READ(ir,*) IST,IEND
         WRITE(iw,*) IST,IEND
            DO 607 J=IST,IEND
            ird=0
               DO 608 I=1,NC
                  IF (IGPC(I,J).eq.1) ird=1
 608           CONTINUE
               if (ird.eq.1) goto 607
               KS=KS+1
               kks(K,KS)=J
 607        CONTINUE
          WRITE(iw,*) ' ANY MORE? 1=YES '
          READ(ir,11) NN
 11       format(I2)
          IF(NN.EQ.1) GOTO 44
          NSYN(K) = KS
          WRITE(iw,*) 'type, number of segments',K,NSYN(K)
 43       WRITE(iw,*) 'indicate CONDUCTANCE NUMBER YOU WANT'
          READ(ir,*)  i
          write(iw,*) i
          irc(i,k)=1
          WRITE(iw,*) ' ANY MORE? 1=YES '
          READ(ir,11) NN
          IF(NN.EQ.1) GOTO 43
 701      CONTINUE
        RETURN
        END
C  ***************************************************
       subroutine linv(n,V,m)
C        numerical inversion of laplace transforms
C
C        ref:  harald stehfest,  cacm vol 13 #1 jan 1970
C
C        if a laplace transform p(s) is given in the form of a
C         real function, linv prodices an approximate value fa
C         of the inverse f(.) at t.
C
C        n must be even
C        at the first call m may be any integer not = n
C
C     REVISED 8/3/83 TO DO JUST THE INITIALIZATION.  THE ACTUAL
C     INVERSION IS DONE IN THE ROOT MODULE.
C
       implicit real*8(a-h,o-z)
       dimension v(N),gA(24),h(12)
       integer sn
       if(m.eq.n) goto 1
       gA(1)=1.0d0
       nh=n/2
       do 2 i=1,n
   2      gA(i+1)=gA(i)*dfloat(i)
       h(1)=2.0d0/gA(nh)
       do 3 i=2,nh
  3      h(i)=dfloat(i)**dfloat(nh)*gA(2*i+1)/(gA(nh-i+1)*
     c      gA(i+1)*gA(i))
       iarg=nh-(nh/2)*2
       if(iarg)9,10,9
  9    sn=2*isign(1,iarg)-1
       goto 11
 10    sn=-1
 11    do 4 i=1,n
          v(i)=0.0d0
          linf=(i+1)/2
          lsup=i
          if(i.gt.nh) lsup=nh
          do 5 k=linf,lsup
             imk=i-k+1
             k2m1=k+k-i+1
  5          v(i)=v(i)+h(k)/(gA(imk)*gA(k2m1))
          v(i)=dfloat(sn)*v(i)
          sn=-sn
  4    continue
       m=n
  1    f1=0.0d0
       return
       end
C ***************************************************
       subroutine getpqr(s,Gg,H,alamda,al,RMM,CM,r,cur,
     C   NVC,NIP,TEXP,NCOND,NSEG,NEQ,
     C   G0,vvso,IICOND,NCURLO,NCIP,iw,
     C   SIGNUM,GP,GPR,irepl,
     C   IA,JA,IRR,IC,ICC,IEL,n3,ISW)
       IMPLICIT real*8 (A-H,O-Z)
       dimension alamda(NSEG),al(NSEG),r(NSEG),CM(NSEG),
     C   RMM(NSEG),GG(NEQ),G0(NEQ),SIGNUM(NSEG),GP(NSEG),
     C   GPR(NSEG),irepl(nseg),
     C   Z(990),SI(990),TA(990),GPCC(990),GPRC(990),
     C   JA(IEL),RTEMP(2191),A(3990),ITEMP(3094),
     C   B(991),IRR(NEQ),ICC(NEQ),IC(NEQ),IA(NEQ+1),
     C   TC(990),C(990),RTW(990),SIG(990),BB(990),
     C   CZS(990),CZT(990)
       EXTERNAL NSPIV
c
       MAX=1290
c
c      MAX is the maximum number of off diagonal elements to store
c      RTEMP size = N + MAX where N is number of equations
c      ITEMP size = 2*N + MAX + 2
c
c      write(6,*) ' in getpqr '
       DO 58 I=1,NSEG
          rtw(I)=dsqrt(RMM(I)*CM(I)*s+1.0d0)
          TC(I)=AL(I)/ALAMDA(I)*RTw(I)
          Z(I)=(RMM(I)*CM(I))/(1.0d0+RMM(I)*CM(I)*S)
          TA(I)=dTANH(TC(I))
          SI(I)=dSINH(TC(I))
          sig(i)=SIGNUM(I)/(s*rtw(i)*rtw(i))
          GPCC(I)=GP(I)*R(I)*ALAMDA(I)/RTW(I)
          GPRC(I)=GPR(I)*R(I)*ALAMDA(I)/RTW(I)
  58   CONTINUE
       DO 1401 J=1,NSEG
          B(J)=0.d0
          BB(J)=-SIG(J)*(1.d0/SI(J)-1.d0/TA(J))
 1401  CONTINUE
       B(NEQ)=0.d0
C
c       call getc(C,rtw,alamda,r,nseg,IA,JA,IEL)
C ********************************************************
c      SUBROUTINE getc(C,rtw,alamda,r,nseg,IA,JA,IEL)
c      IMPLICIT real*8(A-H,O-Z)
c      DIMENSION RTW(NSEG),ALAMDA(NSEG),R(NSEG),
c     C   IA(NSEG+2),JA(IEL),C(NSEG)
C
       c(1)=1.d0
       DO 300 I=2,NSEG
          IP=JA(IA(I))
          c(I)=(RTw(I)/RTw(IP))*(ALAMDA(IP)/ALAMDA(I))*(R(IP)/R(I))
 300   CONTINUE
C
c       RETURN
c       END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       CALL GETA(A,SI,TA,C,NSEG,IEL,IA,JA,irepl)
C ********************************************************
c      SUBROUTINE GETA(A,S,T,C,NSEG,IEL,IA,JA,irepl)
c      IMPLICIT real*8(A-H,O-Z)
c      DIMENSION A(IEL),T(NSEG),S(NSEG),C(NSEG),irepl(nseg),
c     C   IA(NSEG+2),JA(IEL)
C
      NEQ=NSEG+1
      DO 100 I=1,NSEG
         IST=IA(I)
         IEND=IA(I+1)-1
         ist1=ist+1
         JS=JA(ist1)
         A(IST)=-1.0d0/SI(JS)
         su=0.d0
         KK=1
         DO 150 K=IST+2,IEND
            IDA=JA(IST1+KK)
            su=su+C(IDA)/TA(IDA)*dfloat(irepl(ida))
            A(IST1+KK)=-C(IDA)/SI(IDA)*dfloat(irepl(ida))
            KK=KK+1
 150     CONTINUE
         a(ist1)=1.0d0/tA(js)+su
 100  CONTINUE
      A(IA(NEQ))=1.0d0/TA(1)
      A(IA(NEQ)+1)=-1.0d0/SI(1)
c      RETURN
c      End
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
       IF (NCOND.EQ.1) THEN
c          CALL BCOND(B,BB,C,NSEG,IEL,IA,JA,irepl)
C ********************************************************
c      SUBROUTINE BCOND(B,BB,C,NSEG,IEL,IA,JA,irepl)
c      IMPLICIT real*8(A-H,O-Z)
c      DIMENSION B(NSEG+1),BB(NSEG),C(NSEG),irepl(nseg),
c     C   IA(NSEG+2),JA(IEL)
C
c      write(6,*) ' in bcond '
      NEQ=NSEG+1
      DO 400 I=1,NSEG
         IST=IA(I)
         IEND=IA(I+1)-1
         KK=1
         B(I)=BB(I)
         DO 450 K=IST+2,IEND
            IDA=JA(IST+1+KK)
            B(I)=B(I)+C(IDA)*BB(IDA)*dfloat(irepl(ida))
            KK=KK+1
 450     CONTINUE
 400  CONTINUE
      B(NEQ)=bB(1)
c      RETURN
c      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          DO 2200 J=1,NSEG
             B(J)=B(J)+GPRC(J)/S
             IAJP1=IA(J)+1
             A(IAJP1)=A(IAJP1)+GPCC(J)
c             if(n3.eq.1) then
c             if(gprc(j).gt.1.d-8) write(iw,*) ' gprc(j) ',j,gprc(j),
c     c          gpr(j)
c             if(gpcc(j).gt.1.d-8) write(iw,*) ' gpcc(j) ',j,gpcc(j),
c     c          gp(j)
c             endif
 2200     CONTINUE
c     added 3/31/93 cccccccccccccccccccccccccccc
          b(neq)=b(neq)+gprc(1)/s
          a(ia(neq))=a(ia(neq))+gpcc(1)
ccccccccccccccccccccccccccccccccccccc
       ENDIF
C
       IF (IICOND.EQ.1) then
c          CALL SETIC (SI,TA,C,G0,Z,TC,B,NSEG,NEQ,IA,JA,IEL,irepl)
C ************************************************************
c       SUBROUTINE SETIC (S,T,C,G0,Z,TC,B,NSEG,NEQ,IA,JA,IEL,irepl)
c       IMPLICIT real*8 (A-H,O-Z)
c       DIMENSION S(NSEG),T(NSEG),C(NSEG),Z(NSEG),TC(NSEG),
c     C   IA(NEQ+1),JA(IEL),G0(NEQ),B(NEQ),irepl(nseg),
c     C   CZS(990),CZT(990)
       DO 500 I=1,NSEG
          CZS(I)=Z(I)*(1.0d0/SI(I)-1.0d0/TC(I))
          CZT(I)=Z(I)*(1.0d0/TA(I)-1.0d0/TC(I))
  500     CONTINUE
       DO 501 I=1,NSEG
          IST=IA(I)
          IEND=IA(I+1)-1
          JP=JA(IST)
          B(I) = B(I) -G0(JP)*CZS(I) + G0(I)*CZT(I)
          KK=1
          DO 550 K=IST+2,IEND
             ist1k=ist+1+kk
             IDA=JA(IST1K)
             B(I) = B(I) + G0(I)*C(IDA)*CZT(IDA)*dfloat(irepl(ida)) -
     C           G0(IDA)*C(IDA)*CZS(IDA)*dfloat(irepl(ida))
             KK=KK+1
 550      CONTINUE
 501   CONTINUE
       B(NEQ) = B(NEQ) +G0(NEQ)*CZT(1) -G0(1)*CZS(1)
c       RETURN
c       END
       endif
C
       IF(NVC.NE.1) GOTO 88
          B(1)=B(1)+(CUR/SI(1))/S
          A(1)=0.d0
ccccccccccccccccccccccccccccccccccccc 3/31/93
          A(IEL)=0.d0
cold      B(NEQ)=(CUR/TA(1))/S
          vvso=(rtw(1)/(r(1)*alamda(1)))*(a(iel-1)*cur/s-b(neq))
          a(iel-1)=1.d0
          b(neq)=cur/s
ccccccccccccccccccccccccccccccccccccc
          GOTO 89
C
 88    IF(NCIP.NE.1) GOTO 89
       IF (NIP.EQ.1) GOTO 901
       IF (NIP.EQ.2) GOTO 902
       IF (NIP.EQ.3) GOTO 903
       IF (NIP.LT.1.OR.NIP.GT.3) GOTO 89
 901   IF(NCURLO.NE.0) B(NCURLO)=B(NCURLO)+(alamda(NCURLO)/
     C     RTW(NCURLO))*(r(NCURLO)*h)
       IF(NCURLO.EQ.0) B(NEQ)=B(NEQ)+R(1)*H*ALAMDA(1)/RTW(1)
       GOTO 89
 902   IF (NCURLO.NE.0) B(NCURLO)=B(NCURLO)+(alamda(NCURLO)/
     C     RTW(NCURLO))*(r(NCURLO)*h)/s
       IF (NCURLO.EQ.0) B(NEQ)=B(NEQ)+R(1)*H*ALAMDA(1)/RTW(1)/S
       GOTO 89
 903   IF (NCURLO.NE.0) B(NCURLO)=B(NCURLO)+(alamda(NCURLO)/
     C     RTW(NCURLO))*(r(NCURLO)*h)
C    C     *TEXP*dEXP(1.0d0)/((S+TEXP)*(S+TEXP))
     C     *TEXP*TEXP/((S+TEXP)*(S+TEXP))
       IF (NCURLO.EQ.0) B(NEQ)=B(NEQ)+R(1)*H*ALAMDA(1)/RTW(1)
C    C     *TEXP*dEXP(1.0d0)/((S+TEXP)*(S+TEXP))
     C         *TEXP*TEXP/((S+TEXP)*(S+TEXP))
       GOTO 89
C
C
 89     CONTINUE
c       write(6,*) ' entering nspiv '
        CALL NSPIV (NEQ,IEL,IA,JA,A,B,MAX,IRR,ICC,IC,
     C       GG,ITEMP,RTEMP,ierR)
        if (ierR.LT.0.OR.ISW.NE.1) then     
        write(iw,4200) ierR,ISW
        ISW=1
c       write(6,4200) ierR
 4200   format(' ieRr=  ',2i8)
        endif
 4202   CONTINUE
ccccccccccccccccccccccccccccc3/31/93
       if(nvc.eq.1) then
          vvso=vvso-rtw(1)/(r(1)*alamda(1))*gg(1)/si(1)
       endif
cccccccccccccccccccccccccccccc
       return
       END
C  **********************************************
       SUBROUTINE NSPIV (N,IEL,IA,JA,A,B,MAX,RR,CC,IC,X,ITEMP,RTEMP
     C   ,IERR)
       real*8 A(IEL),B(N),X(N),RTEMP(N+MAX)
       INTEGER IA(N+1),JA(IEL),RR(N),CC(N),IC(N),ITEMP(2*N+MAX+2)
       INTEGER IU,JU,U,Y,P
       Y=1
       U=Y+N
       P=1
       IU=P+N+1
       JU=IU+N+1
       CALL NSPIV1 (N,IEL,IA,JA,A,B,MAX,RR,CC,IC,X,RTEMP(Y),ITEMP(P),
     C              ITEMP(IU),ITEMP(JU),RTEMP(U),IERR)
       RETURN
       END
C
C
       SUBROUTINE NSPIV1 (N,IEL,IA,JA,A,B,MAX,RR,CC,IC,X,Y,P,IU,JU,U
     C     ,IERR)
       real*8 A(IEL),B(N),U(MAX),X(N),Y(N)
       real*8 DK,LKI,ONE,XPV,XPVMAX,YK,ZERO
       INTEGER CC(N),IA(N+1),IC(N),IU(N+1),JA(IEL),JU(MAX),P(N+1),RR(N)
       INTEGER CK,PK,PPK,PV,V,VI,VJ,VK
       IF (N.EQ.0) GOTO 1001
       ONE=1.0
       ZERO=0.0
       DO 10 J=1,N
          X(J)=ZERO
 10       CONTINUE
       IU(1)=1
       JUPTR=0
       DO 170 K=1,N
          P(N+1)=N+1
          VK=RR(K)
          JMIN=IA(VK)
          JMAX=IA(VK+1)-1
          IF (JMIN.GT.JMAX) GOTO 1002
          J=JMAX
 20          JAJ=JA(J)
             VJ=IC(JAJ)
             X(VJ)=A(J)
             PPK=N+1
 30          PK=PPK
             PPK=P(PK)
             IF(PPK-VJ) 30,1003,40
 40          P(VJ)=PPK
             P(PK)=VJ
             J=J-1
             IF(J.GE.JMIN) GOTO 20
          VI=N+1
          YK=B(VK)
 50       VI=P(VI)
          IF(VI.GE.K) GOTO 110
          LKI=-X(VI)
          X(VI)=ZERO
          YK=YK+LKI*Y(VI)
          PPK=VI
          JMIN=IU(VI)
          JMAX=IU(VI+1)-1
          IF(JMIN.GT.JMAX) GOTO 50
          DO 100 J=JMIN,JMAX
             JUJ=JU(J)
             VJ=IC(JUJ)
             IF (X(VJ).NE.ZERO) GOTO 90
             IF (VJ-PPK) 60,90,70
 60          PPK=VI
 70          PK=PPK
             PPK=P(PK)
             IF(PPK-VJ) 70,90,80
 80          P(VJ)=PPK
             P(PK)=VJ
             PPK=VJ
 90          X(VJ)=X(VJ)+LKI*U(J)
 100         CONTINUE
          GOTO 50
 110      IF (VI.GT.N) GOTO 1004
          XPVMAX=dABS(X(VI))
          MAXC=VI
          NZCNT=0
          PV=VI
 120         V=PV
             PV=P(PV)
             IF(PV.GT.N) GOTO 130
             NZCNT=NZCNT+1
             XPV=dABS(X(PV))
             IF (XPV.LE.XPVMAX) GOTO 120
             XPVMAX=XPV
             MAXC=PV
             MAXCL=V
             GOTO 120
 130      IF (XPVMAX.EQ.ZERO) GOTO 1004
          IF (VI.EQ.K) GOTO 140
          IF (VI.EQ.MAXC) GOTO 140
          P(MAXCL)=P(MAXC)
          GOTO 150
 140      VI=P(VI)
 150      DK=ONE/X(MAXC)
          X(MAXC)=X(K)
          I=CC(K)
          CC(K)=CC(MAXC)
          CC(MAXC)=I
          CK=CC(K)
          IC(CK)=K
          IC(I)=MAXC
          X(K)=ZERO
          Y(K)=YK*DK
          IU(K+1)=IU(K)+NZCNT
          IF (IU(K+1).GT.MAX+1) GOTO 1005
          IF(VI.GT.N) GOTO 170
          J=VI
 160         JUPTR=JUPTR+1
             JU(JUPTR)=CC(J)
             U(JUPTR)=X(J)*DK
             X(J)=ZERO
             J=P(J)
             IF (J.LE.N) GOTO 160
 170      CONTINUE
       K=N
       DO 200 I=1,N
          YK=Y(K)
          JMIN=IU(K)
          JMAX=IU(K+1)-1
          IF (JMIN.GT.JMAX) GOTO 190
          DO 180 J=JMIN,JMAX
             JUJ=JU(J)
             JUJ=IC(JUJ)
             YK=YK-U(J)*Y(JUJ)
 180         CONTINUE
 190      Y(K)=YK
          CK=CC(K)
          X(CK)=YK
          K=K-1
 200      CONTINUE
       IERR=IU(N+1)-IU(1)
       RETURN
 1001  IERR=0
       RETURN
 1002  IERR=-K
       RETURN
 1003  IERR=-(N+K)
       RETURN
 1004  IERR=-(2*N+K)
       RETURN
 1005  IERR=-(3*N+K)
       RETURN
       END
