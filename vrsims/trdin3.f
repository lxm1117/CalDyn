       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION X(10),xmax(10),txmax(10),hz(20),f0(3),cca(3),
     c   ccamax(3),tcamax(3),sumf(3),ispn(20),brp(20),xcpre(2)
       CHARACTER*20 FILE1(20),file2(20),filint
c      data file1/'ns0241.ha1','ns0241.ha2','ns0241.ha3',
c    c	 'ns0216.ha1','ns0216.ha2','ns0216.ha3','vsv010.sav'/
c      data file2/'os0241.ha1','os0241.ha2','os0241.ha3',
c    c	 'os0216.ha1','os0216.ha2','os0216.ha3','vsv010.sav'/
c      data hz/400,400,400,400,400,400,200/
c      data ispn/114,114,114,114,114,114,005/
c
       iw=6
       ir=40
       open(file='trddat.dat',unit=40,status='old')
c      open(file='junk.dat',unit=49,status='unknown')
c
       write(iw,*) ' enter number of files '
       read(ir,*) inf
       write(iw,*) ' enter ',inf,' input file names, one per line '
       do 8 i=1,inf
          read(ir,10) file1(i)
 10       format (a20)
 8     continue
       write(iw,*) ' enter ',inf,' output file names, one per line '
       do 9 i=1,inf
          read(ir,10) file2(i)
 9     continue
       write(iw,*) ' enter int file name '
       read(ir,10) filint
c
       write(iw,*) ' enter ',inf,' hz values separated by commas '
       read(ir,*) (hz(i),i=1,inf)
       write(iw,*) ' enter ',inf,' spine #s,separated by commas '
       read(ir,*) (ispn(i),i=1,inf)
       write(iw,*) ' enter ',inf,' base rest voltage sep by commas '
       read(ir,*) (brp(i),i=1,inf)
c
c
	  anai=10.d0
	  anao=145.d0
	  aki=140.d0
	  ako=2.5d0
	  cai=.00002d0
	  cao=2.0d0
	  rtof=26.73d0
	  f=96.48456d0
          actk=.75d0
          actna=.75d0
          actca=.19d0
          anai=anai*actna
          anao=anao*actna
          aki=aki*actk
          ako=ako*actk
          cai=cai*actca
          cao=cao*actca
c
       vk= rtof*dlog(ako/aki)
       vna= rtof*dlog(anao/anai)
       vca= rtof/2.d0*dlog(cao/cai)
c      write(6,*) ' vk,vna,vca,rtof= ',vk,vna,vca,rtof
c      write(6,*) ' nai,noo,ki,ko,cai,cao= ',
c    c     anai,anao,aki,ako,cai,cao
c
       open(file=filint,unit=33,status='unknown')
       nd=9
       do 900 l=1,inf
	  OPEN (FILE=FILE1(l),UNIT=31,status='old')
	  OPEN (FILE=FILE2(l),UNIT=34,status='unknown')
	  do 5 j=1,nd
	     xmax(j)=0.d0
	     txmax(j)=0.d0
 5	  continue
          do 6 j=1,3
 	     sumf(j)=0.d0
             ccamax(j)=0.d0
             tcamax(j)=0.d0
             f0(j)=0.d0
 6        continue
	  T0=0.D0
C
          xcpre(1)=x(2)
          xcpre(2)=x(6)
 1	  READ(31,*,END=9999) T1,(X(I),I=1,ND)
c
          do 400 k=1,2
             ii=4*k-2
c                save previous voltage, the current was computed
c                using the previous voltage.  This becomes important
c                when voltage crosses synaptic reversal potential
             v=xcpre(k)+brp(l)
c            v=x(ii)+brp(l)
c
             de1=1.d0-dexp(-v/rtof)
             de2=1.d0-dexp(-2.d0*v/rtof)
c
	     c=f/rtof*v/de1*(anai-anao*dexp(-v/rtof))
c    c	       (1.d0-dexp(-v/rtof))
c
	     d=f/rtof*v/de1*(aki-ako*dexp(-v/rtof))
c    c	       (1.d0-dexp(-v/rtof))
c
	     e=4.d0*f/rtof*v/de2*(cai-cao*dexp(-2.d0*v/rtof))
c    c	       (1.d0-dexp(-2.d0*v/rtof))
c
	     pna=x(ii+3)/(c+1.0d0*d+10.6d0*e)
	     pk=pna*1.0d0
	     pca=10.6d0*pna
	     cna=c*pna
	     ck=d*pk
	     cca(k)=e*pca
c            if(cca(k).lt.0) write(6,*) 'c,d,e= ',c,d,e,de1,de2
c            sum=c+d+10.6*e
c            if(k.eq.1) write(49,49) t1,v,x(ii+3),cna,ck,cca(1),
c    c          c,d,e,sum
c 49         format(f8.3,9e12.5)
c
c                why was this next statement here?  There must
c                have been a reason.
c 	     IF(cca(k).LE.1.0D-4) cca(k)=0.D0
	     if(cca(k).gt.ccamax(k)) then
	        ccamax(k)=cca(k)
	        tcamax(k)=t1
	     endif
	     SUMF(k)=SUMF(k)+0.5D0*(cca(k)+F0(k))*(T1-T0)
	     F0(k)=cca(k)
             xcpre(k)=x(ii)
c
 400     continue
c
         gn1=x(3)*1.d3
         gnn1=x(4)*1.d3
         gn2=x(7)*1.d3
         gnn2=x(8)*1.d3
         v1=x(2)+brp(l)
         v2=x(6)+brp(l)
	 write(34,23,err=3333) t1,x(1),v1,gn1,gnn1,cca(1),
     c      v2,gn2,gnn2,cca(2)
 23	 format(f7.2,f9.3,2(f9.3,2f8.2,e12.4))
c
	  do 20 j=1,nd
	     if(x(j).ge.xmax(j)) then
		xmax(j)=x(j)
		txmax(j)=t1
	     endif
 20	  continue

	  T0=T1
	  GOTO 1
C
 9999	  continue
c
	  write(33,109) ispn(l),HZ(l)
 109	  format(' THERE WERE ',I4,' SPINE INPUTS AT ',F6.0,' HZ ')
	  WRITE(33,110) SUMF(1),SUMF(2)
 110	  FORMAT(' SUMF= ',2F15.8)
	  write(33,120) xmax(1),xmax(2),xmax(3),xmax(4),ccamax(1),
     c      xmax(6),xmax(7),xmax(8),ccamax(2)
	  write(33,121) txmax(1),txmax(2),txmax(3),txmax(4),
     c      tcamax(1),
     c      txmax(6),txmax(7),txmax(8),tcamax(2)
 120	  format('  xmax= ',f7.2,2(f7.2,2f7.4,f7.4))
 121	  format(' txmax= ',f7.1,2(f7.1,2f7.1,f7.1))
C
	  CLOSE(UNIT=31)
	  CLOSE(UNIT=34)
c
 900   continue
 3333  write(6,*) t1,x(1)
c
       STOP
       END
