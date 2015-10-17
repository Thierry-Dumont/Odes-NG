      subroutine jactbb(N,X,Y,DFY,LDFY,RPAR,IPAR)
      implicit none
      integer n,ipar,ldfy,mujac,i
      double precision Y(n),X,DFY(LDFY,N),rpar(*),h,eps,eh2
      h=100./(n-1)
      eps=0.001d0
      eh2=eps/(h*h)
      mujac=1
c     DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
c      u*(2.-3.*u);
      DFY(MUJAC+1,1)=-2.d0*eh2+ y(1)*(2.d0-3.d0*y(1))
      DFY(1-2+MUJAC+1,2) =2.d0*eh2
      do i=2,n-1
         DFY(MUJAC+1,i) =-2.d0*eh2+y(i)*(2.d0-3.d0*y(i))
         DFY(1+MUJAC+1,i-1)=eh2
         DFY(MUJAC,i+1)=eh2
      enddo
      DFY(MUJAC+1,n)=-2.d0*eh2+ y(n)*(2.d0-3.d0*y(n))
      DFY(1+MUJAC+1,n-1)=2.d0*eh2
      end
