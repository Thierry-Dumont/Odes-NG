  
         subroutine F(N,X,Y,R,RPAR,IPAR)
         implicit none
         double precision Y(3),R(3),X,RPAR,h,eps,eh2
         integer N,IPAR,i
         double precision function f
         h=100./(n-1)
         eps=0.001;
         eh2=eps/(h*h);
         r(1)=2*eh2*(y(2)-y(1))+y(1)*y(1)*(1.d0-y(1))
         do i=2,n-1
             r(i)=eh2*(y(i-1)-2*y(i)+y(i+1))+y(i)*y(i)*(1.d0-y(i))
         enddo
         r(n)=2*eh2*(y(n-1)-y(n))++y(n)*y(n)*(1.d0-y(n))
         end
