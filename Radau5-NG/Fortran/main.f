      program radau
      implicit double precision (a-h,o-z)
      dimension y(1000),RPAR(3),work(1000000),iwork(1000000),rtol(1)
      dimension atol(1)
      external F,jactbb
      n=1000


      itol=0
c
      ijac=1
c
      mljac=1
      mujac=1
      imas=0
      mlmas=n
      mumas=0
      iout=0
      nwork=1000
      do i=1,20
         work(i)=0.0
         iwork(i)=0
      enddo

      
      lwork=100000
      liwork=100000
      rtol(1)=1.d-5
      atol(1)=1.d-5
      do i=1,10
         do j=1,n/2
            y(j)=0.9999999
         enddo
         do j=n/2+1,n
            y(j)=0.0
         enddo
         x=0
         h=1.0d0
         xend=500.d0
         IWORK(1)=0
         iout=0
         call radau5(n,F,x,y,xend,h,rtol,atol,itol,jactbb,
     &        IJAC,MLJAC,MUJAC,
     &        MAS ,IMAS,MLMAS,MUMAS,
     &        SOLOUT,IOUT,
     &        WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
c         write(*,*) "nstep=",IWORK(16),"njac=",IWORK(15)
c         write(*,*) "naccpt=",IWORK(17),"nrejct=",IWORK(18)
c         write(*,*) "ndec=",IWORK(19)
      enddo
      write(*,*) "nstep=",IWORK(16),"njac=",IWORK(15)
      write(*,*) "naccpt=",IWORK(17),"nrejct=",IWORK(18)
      write(*,*) "ndec=",IWORK(19)
c      write(*,*) "y= ",y
      write(*,*) "idid=",idid," x=",x
      write(*,*)"fin."
      end
      subroutine solout()
      end
