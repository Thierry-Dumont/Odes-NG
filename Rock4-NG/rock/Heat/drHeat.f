c * * * * * * * * * * * * * * * * * * * * * * * * *
c    Driver for ROCK4 (or ROCK2) at Brusselator-2dim problem
c * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    This driver shows how to use ROCK4 (or ROCK2). It solves a
c    system of ODEs resulting from the 2-dimensional space 
c    discretization of the Brusselator equations (u=u(x,y,t),v=v(x,y,t)):
c     
c    u_t=1+u^2*v-4.4*u+0.1*(u_{xx}+u_{yy})+f(x,y,t)
c    v_t=3.4*u-u^2*v+0.1*(v_{xx}+v_{yy})     for t>=0, 0<= x <= 1, 0<= y <= 1
c
c    with initial conditions 
c
c    u(x,y,0)=22*y*(1-y)^{3/2}  v(x,y,0)=27*x*(1-x)^{3/2}
c
c    and periodic boundary conditions
c
c    u(x+1,y,t)=u(x,y,t),  v(x,y+1,t)=v(x,y,t).
c
c    The function f is defined by (inhomogeneity)
c
c    f(x,y,t)=5 if (x-0.3)^2+(y-0.6)^2<= 0.1^2 and t>=1.1
c            =0 else
c
c    We discretize the space variables with
c    x_i=i/(N+1), y_i=i/(N+1) for i=0,1,...,N,
c    with N=128. We obtain a system of 32768
c    equations. The spectral radius of the Jacobian 
c    can be estimated with the Gershgorin theorem   
c    (13200 is an estimation for it). Thus we 
c    provide an external function RHO, giving 
c    the spectral radius of the Jacobian matrix.
c    As output point we choose t_out=1.5 
c
c--------------------------------------------------------
c ----- to integrate with rock2.f ----- 
c     include 'rock2.f'
c  
      include 'rock4.f'     
      implicit double precision (a-h,o-z)
      parameter (neqn=5000)
c--------------------------------------------------------
c      Work is of length 7*neqn because the spectral radius 
c      is computed externally (otherwise should be 8*neqn).
c      If integrating with ROCK2 define work(4*neqn).
c-------------------------------------------------------- 
c ----- to integrate with rock2.f 
c     dimension y(neqn),work(4*neqn) 
c 
      dimension y(neqn),work(7*neqn)
      integer iwork(12),idid
      external fbrus
c --- common parameters for the problem -----
      common/trans/anu,uh2,rspec
c ----- file for solution -----
       open(8,file='sol.out')
       rewind 8
c 
       anu=0.01
       h=1.d0/(neqn-1)
       uh2=anu/(h*h)
       rspec=2.d0*uh2
c--------------------------------------------------------
c     Initialize iwork: 
c      iwork(1)=1  RHO returns an upper bound for the spectral radius.
c      iwork(2)=1  The Jacobian is constant (RHO is called once).
c      iwork(3)=0  Return and solution at tend.
c      iwork(4)=0  Atol and rtol are scalars.
c--------------------------------------------------------
      iwork(1)=1
      iwork(2)=1
      iwork(3)=0
      iwork(4)=0
c ----- initial and end point of integration -----
      t=0.0d0
      tend=0.1
c ----- initial values -----
      do i=1,neqn
         if(i.le.neqn/2) then
            y(i)=1.
         else
            y(i)=0.
         endif
      enddo
c ----- required tolerance -----
      rtol=1.d-5
      atol=rtol
c ----- initial step size -----
      h=0.1
c ----- integration -----
      write(6,*) 'Integration of the 2-dim Brusselator problem' 
c ----- to integrate with rock2.f     
c      call rock2(neqn,t,tend,h,y,fbrus,atol,rtol,work,
c     &           iwork,idid) 
c
c ----- call of the subroutine rock4 -----
      call rock4(neqn,t,tend,h,y,fbrus,atol,rtol,work,
     &           iwork,idid)    
c ----- print solution -----
      do j=1,neqn
        write (8,*) y(j)
      end do
c ----- print statistics -----
      write(6,*) 'Solution is tabulated in file sol.out'
      write(6,*) 'The value of IDID is',idid
      write(6,*) 'Max estimation of the spectral radius=',iwork(11)
      write(6,*) 'Min estimation of the spectral radius=',iwork(12)
      write(6,*) 'Max number of stages used=',iwork(10)
      write(6,*) 'Number of f eval. for the spectr. radius=',iwork(9)
      write (6,91) iwork(5),iwork(6),iwork(7),iwork(8)
 91   format(' Number of f evaluations=',i10,' steps=',i10,
     &        ' accpt=',i10,' rejct=',i10) 
c--------------------------------------------------------
c     End of main program
c--------------------------------------------------------
      end      
c--------------------------------------------------------
c     The subroutine RHO gives an estimation of the spectral 
c     radius of the Jacobian matrix of the problem. This
c     is a bound for the whole interval and thus RHO is called
c     once.
c--------------------------------------------------------
      double precision function rho(neqn,t,y)
      implicit double precision (a-h,o-z)
      common/trans/anu,uh2,rspec
      rho=rspec
      write(*,*) "RHO",rho
      return
      end 
c--------------------------------------------------------
c     The subroutine FBRUS compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine fbrus(neqn,x,y,f)
c ----- brusselator with diffusion in 2 dim. space -----
      implicit double precision (a-h,o-z)
      dimension y(neqn),f(neqn)
      common/trans/anu,uh2,rspec
       f(1)=uh2*(y(2)-y(1))
      do i=2,neqn-1
         f(i)=uh2*(y(i-1)-2.0*y(i)+y(i+1))
      enddo
      f(neqn)=uh2*(y(neqn-1)-y(neqn))
c      write(*,*) uh2
c      write(*,*)(f(i),i=1,neqn)
c      read(*,*) i
      return
      end  
      
      


