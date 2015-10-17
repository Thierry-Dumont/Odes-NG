      program radau
      implicit double precision (a-h,o-z)
      dimension y(3),RPAR(3),work(1000),iwork(1000),rtol(1),atol(1)
      external F
      n=3


      itol=0
      ijac=0
      mljac=n
      mujac=0
      imas=0
      mlmas=n
      mumas=0
      iout=0
      nwork=1000
      do i=1,20
         work(i)=0.0
         iwork(i)=0
      enddo

      
      lwork=1000
      liwork=1000
      rtol(1)=1.d-5
      atol(1)=1.d-5
      do i=1,1000
         y(1)=1.0d0
         y(2)=2.0d0
         y(3)=3.0d0
         x=0
         h=1.0d0
         xend=500.d0
         IWORK(1)=0
         call radau5(n,F,x,y,xend,h,rtol,atol,itol,jac,
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
      write(*,*) "y= ",y
      write(*,*) "idid=",idid," x=",x
      write(*,*)"fin."
      end
      subroutine solout()
      end
