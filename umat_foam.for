!************************************************************************
!
! User material subroutine (UMAT) for the large-deformation behavior of 
!  compressible, elastomeric foams, using a logarithmic-strain-based
!  hyperelastic model. This UMAT is not for use in plane stress or in any 
!  other situation in which there are more strain terms than stress terms.
!
! Xiuqi Li, March 2019
!
!************************************************************************
! Usage:
!************************************************************************
!
!     Material Properties Vector (*user material, constants = 14)
!     --------------------------------------------------------------
!     G0      = props(1)  ! Shear modulus
!     B       = props(2)  ! Bulk modulus
!     Jmin    = props(3)  ! Locking value of J in f-fuction
!     C1      = props(4)  ! C_1 in L-function
!     K10     = props(5)  ! K1_0 in X-function
!     delta_K = props(6)  ! delta_K in X-function
!     X1      = props(7)  ! X1' in X-function
!     X2      = props(8)  ! X2' in X-function
!     C0      = props(9)  ! C_0 in L-function
!     p       = props(10) ! power p in L-function
!     q       = props(11) ! power q in L-function
!     C2      = props(12) ! C_2 in f-function
!     C3      = props(13) ! C_3 in f-function
!     r       = props(14) ! power r in f-function
!
!     State Variables (*depvar 3) (for visualization purposes only)
!     --------------------------------------------------------------
!     statev(1) = K1 ---- ln(J), range(-inf, inf)
!     statev(2) = K2 ---- amount of distortional deformation, range[0, inf)
!     statev(3) = K3 ---- deformation mode, range[-1, 1], 
!                             -1 is compression, 0 is shear, 1 is tension
!
!************************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     +     rpl,ddsddt,drplde,drpldt,
     +     stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     +     ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     +     celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      !
      include 'aba_param.inc'
      !
      dimension stress(ntens),statev(nstatv),
     +     ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     +     stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     +     props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
      !
      character*8 cmname
      !
      ! Variables defined and used in the UMAT
      !
      integer i,j,stat
      !
      real*8 F_tau(3,3),T_tau(3,3),K1,K2,K3,psi,SpTanMod(3,3,3,3)
      !
      real*8 zero,one,two,three,fourth,third,half,six,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +     six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0)


      ! Do nothing if a "dummy" step
      !
      if(dtime.eq.zero) return


      ! Get the deformation gradient at the end of the increment. The
      !  subscript tau denotes the time at the end of the increment.
      !
      do i=1,3
	do j=1,3
	  F_tau(i,j) = dfgrd1(i,j)
	end do
      end do
      
      
      ! At the start of an Abaqus calculation, the state
      !  variables are passed into UMAT with zero values.
      !  Initialize the state variables. At this point,
      !  the time total_time and step_time both have a value
      !  equal to zero and the step counter, kstep, is 
      !  equal to 1.
      !
      if ((time(1).eq.zero).and.(kstep.eq.1)) then
        statev(1) = zero
        statev(2) = zero
        statev(3) = zero
      end if


      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !
      ! Perform the constitutive update
      !
      call foam(props,nprops,dtime,F_tau,
     +          T_tau,K1,K2,K3,SpTanMod,psi,stat)
      !
      if (stat.eq.0) then
         pnewdt = 0.5
         return
      endif
      !
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     
     
      ! Update the state variables
      !
      statev(1) = K1
      statev(2) = K2
      statev(3) = K3


      ! Update the stress
      !	
      stress = zero
      !
      do i=1,ndi
        stress(i) = T_tau(i,i)
      end do
      !
      if (nshr.ne.0) then
        stress(ndi+1) = T_tau(1,2)
        if (nshr.ne.1) then
          stress(ndi+2) = T_tau(1,3)
          if (nshr.ne.2) then
            stress(ndi+3) = T_tau(2,3)
          endif
        endif
      endif


      ! Update the tangents
      !	
      ddsdde = zero
      !
      do i=1,ndi
        do j=1,ndi
          ddsdde(i,j) = SpTanMod(i,i,j,j)
        end do
      end do
      !
      if (nshr.ne.0) then
        do i=1,ndi
          ddsdde(i,ndi+1) = SpTanMod(i,i,1,2)
          ddsdde(ndi+1,i) = SpTanMod(1,2,i,i)
        end do
        ddsdde(ndi+1,ndi+1) = SpTanMod(1,2,1,2)
        if (nshr.ne.1) then
          do i=1,ndi
            ddsdde(i,ndi+2) = SpTanMod(i,i,1,3)
            ddsdde(ndi+2,i) = SpTanMod(1,3,i,i)
          end do
          ddsdde(ndi+2,ndi+2) = SpTanMod(1,3,1,3)
          ddsdde(ndi+1,ndi+2) = SpTanMod(1,2,1,3)
          ddsdde(ndi+2,ndi+1) = SpTanMod(1,3,1,2)
          if (nshr.ne.2) then
            do i=1,ndi
              ddsdde(i,ndi+3) = SpTanMod(i,i,2,3)
              ddsdde(ndi+3,i) = SpTanMod(2,3,i,i)
            end do
            ddsdde(ndi+3,ndi+3) = SpTanMod(2,3,2,3)
            ddsdde(ndi+1,ndi+3) = SpTanMod(1,2,2,3)
            ddsdde(ndi+3,ndi+1) = SpTanMod(2,3,1,2)
            ddsdde(ndi+2,ndi+3) = SpTanMod(1,3,2,3)
            ddsdde(ndi+3,ndi+2) = SpTanMod(2,3,1,3)
          endif
        endif
      endif


      ! Update the free energy
      !	
      sse = psi


      return
      end subroutine umat

!************************************************************************
!     Material subroutines
!************************************************************************

      subroutine foam(
     +        ! Inputs
     +        props,nprops,dtime,F_tau,
     +        ! Outputs
     +        T_tau,K1,K2,K3,SpTanMod,psi,stat)

      implicit none
      !
      integer i,j,k,l,m,w,s,t,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),T_tau(3,3),K1,K2,K3,
     +  SpTanMod(3,3,3,3),psi,Iden(3,3),G0,B,Jmin,C1,K10,delta_K,X1,X2,
     +  C0,p,q,C2,C3,r,detF,Rot(3,3),U_tau(3,3),E_tau(3,3),B_tau(3,3),
     +  dev_E(3,3),N(3,3),Y(3,3),f,dfdK1,dfdK1_2,X,dXdK1,dXdK1_2,Lf,
     +  dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,dpsidK1,dpsidK2,dpsidK3,dpsidK1_2,
     +  dpsidK2_2,dpsidK2dK1,dpsidK2dK3,DTkDE(3,3,3,3),DlnBDB(3,3,3,3)
      !
      real*8 zero,one,two,three,fourth,third,half,six,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +     six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0)
 

      ! Identity matrix
      !
      call onem(Iden)
      
      
      ! Obtain material properties
      !
      G0 = props(1) ! Shear modulus
      B = props(2) ! Bulk modulus
      Jmin = props(3) ! Locking value of J in f-fuction
      C1 = props(4) ! C_1 in L-function
      K10 = props(5) ! K1_0 in X-function
      delta_K = props(6) ! delta_K in X-function
      X1 = props(7) ! X1' in X-function
      X2 = props(8) ! X2' in X-function
      C0 = props(9) ! C_0 in L-function
      p = props(10) ! power p in L-function
      q = props(11) ! power q in L-function
      C2 = props(12) ! C_2 in f-function
      C3 = props(13)! C_3 in f-function
      r = props(14) ! power r in f-function
      
      
      ! Calulate kinematic quantities
      !
      call mdet(F_tau,detF)
      call skinem(F_tau,Rot,U_tau,E_tau,stat)
      E_tau = matmul(matmul(Rot,E_tau),transpose(Rot))
      B_tau = matmul(F_tau,transpose(F_tau))
      
      
      ! Calulate the invariants K1, K2, K3 
      !
      K1 = E_tau(1,1) + E_tau(2,2) + E_tau(3,3)	  
      dev_E = E_tau - third*K1*Iden
      K2 = dsqrt(sum(dev_E*dev_E))
      if (K2.eq.zero) then
        ! if K2 == 0, assign an arbitrary value to N, K3, and Y
        N = zero
        K3 = zero
        Y = zero
      else
        ! else compute N, K3, and Y
        N = dev_E/K2
        call mdet(N,K3)
        K3 = K3*three*dsqrt(six)  
        Y  = three*dsqrt(six)*matmul(N,N) 
     +              - dsqrt(six)*Iden - three*K3*N
      end if 


      ! Compute the free energy function and its derivatives
      !
      call f_fun(f,dfdK1,dfdK1_2,K1,Jmin,C2,C3,r)
      call X_fun(X,dXdK1,dXdK1_2,K1,X1,X2,K10,delta_K)
      call L_fun(Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C0,C1,p,q)
      dpsidK1 = G0*dXdK1*K2**two + B*dfdK1
      dpsidK2 = G0*(two*X*K2 + dLdK2)
      dpsidK3 = G0*dLdK3
      dpsidK1_2 = G0*dXdK1_2*K2**two + B*dfdK1_2
      dpsidK2_2 = G0*(two*X + dLdK2dK2)
      dpsidK2dK1 = G0*two*dXdK1*K2
      dpsidK2dK3 = G0*dLdK2dK3
      psi = G0*(X*(K2**two) + Lf) + B*f
      
      
      ! Compute the Cauchy stress and the material contribution
      !   to the tangent (derivative of Kirchhoff stress wrt log strain)
      !
      if (K2.eq.zero) then
        !
        T_tau = (one/detF)*(dPSIdK1*Iden) 
        DTkDE = zero
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DTkDE(i,j,k,l) = DTkDE(i,j,k,l)
     +              + two*B*dfdK1_2*Iden(i,j)*Iden(k,l)
     +              + two*G0*(half*Iden(i,k)*Iden(j,l) 
     +                         + half*Iden(i,l)*Iden(j,k)
     +                            - third*Iden(i,j)*Iden(k,l))
              end do
            end do
          end do
        end do
        !
      else		 		  
        !
        T_tau = (one/detF)*(dpsidK1*Iden + dpsidK2*N + dpsidK3/K2*Y)
        DTkDE = zero
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DTkDE(i,j,k,l) = DTkDE(i,j,k,l) +
     +  (dev_E(k,l)*Iden(i,j) + Iden(k,l)*dev_E(i,j))*
     +  (-two*dsqrt(six)/K2**three*dpsidK3 + one/K2*dpsidK2dK1)+
     +  (dev_E(k,l)*Y(i,j) + Y(k,l)*dev_E(i,j))*
     +  (one/K2**two*dpsidK2dK3 - three/K2**three*dpsidK3)+
     +  dev_E(k,l)*dev_E(i,j)*
     +  (-one/K2**three*dpsidK2 + one/K2**two*dpsidK2_2 -
     +  three*K3/K2**four*dpsidK3)+
     +  (Iden(i,k)*dev_E(l,j) + Iden(j,l)*dev_E(i,k) + 
     +  Iden(i,l)*dev_E(k,j)+Iden(k,j)*dev_E(i,l))
     +  *three/two*dsqrt(six)/K2**three*dpsidK3  
     +  +(half*Iden(i,k)*Iden(j,l) + half*Iden(i,l)*Iden(j,k) - 
     +  third*Iden(i,j)*Iden(k,l))*
     +  (one/K2*dpsidK2 - three*K3/K2**two*dpsidK3)+
     +  Iden(i,j)*Iden(k,l)*dpsidK1_2
              end do
            end do
          end do
        end do
        !
      end if
      
      
      ! Compute the fourth order tensor dlnBDB
      !
      call dlnxdx(B_tau,DlnBDB)	
      !
      ! Compute the tangent required by Abaqus 	 
      !
      SpTanMod = zero
      do i=1,3
        do j=1,3
          do k=1,3 
            do l=1,3
              do m=1,3
                do w=1,3
                  do t=1,3
                    do s=1,3					
                      SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +  fourth/detF*(B_tau(t,l)*Iden(k,s)+B_tau(s,l)*Iden(t,k)+
     +	B_tau(t,k)*Iden(l,s)+B_tau(s,k)*Iden(t,l))
     +  *DTkDE(i,j,m,w)*DlnBDB(m,w,t,s)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
	     
	   	 
      return
      end subroutine foam 

!****************************************************************************

      subroutine f_fun(f,dfdK1,dfdK1_2,K1,Jmin,C2,C3,r)
      !
      ! This subroutine calculates the f-function and its derivatives
      !
      implicit none
      !
      real*8 f,dfdK1,dfdK1_2,K1,Jmin,B,C2,C3,r,J
      !
      real*8 zero,one,two,three,fourth,third,half,six,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +     six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0)
     
     
      J = dexp(K1)
      
      
      ! Compute f
      !
      f = (dexp(C2*(K1))-C2*K1 - one)/(C2**two) + (C3/(r-one))*
     +           (-J**(one-r) + (J-Jmin)**(one-r)*(one-Jmin)**r + Jmin)
      
      
      ! Compute first derivative of f wrt K1  
      !
      dfdK1 = (dexp(C2*(K1))-one)/C2 + 
     + C3*J*(J**(-r)-(one-Jmin)**r/(J-Jmin)**r)
     
     
      ! Compute second derivative of f wrt K1
      !
      dfdK1_2 = dexp(C2*K1) + C3*(J**(-r)-(one-Jmin)**(r)
     + /(J-Jmin)**r)*dexp(K1)
     + + C3*(-r*J**(-r-one) + r*(one-Jmin)**(r)/(J-Jmin)**(r+one))
     +  *dexp(K1)*dexp(K1)


      return
      end subroutine f_fun
      
!****************************************************************************

      subroutine X_fun(X,dXdK1,dXdK1_2,K1,X1,X2,K10,delta_K)
      !
      ! This subroutine calculates the X-function and its derivatives
      !
      implicit none
      !
      real*8 K1,X,dXdK1,X1,X2,delta_K,K10,dXdK1_2
      !
      real*8 zero,one,two,three,fourth,third,half,pi
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0, pi=3.141592653d0)


      ! Compute X		
      X =  (X1+X2)/two*K1 + (X1-X2)/two*delta_K*
     +	     dlog(dcosh((K1-K10)/delta_K)/dcosh(K10/delta_K)) + one


      ! Compute first derivative of X wrt K1
      !
      dXdK1 = (X1+X2)/two + (X1-X2)/two*dtanh((K1 - K10)/delta_K)
      
      
      ! Compute second derivative of X wrt K1
      !
      dXdK1_2 = (X1-X2)/two/dcosh((K10-K1)/delta_K)**two/delta_K	


      return
      end subroutine X_fun
      
!****************************************************************************

      subroutine L_fun(Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C0,C1,p,q)
      !
      ! This subroutine calculates the L-function and its derivatives
      !
      implicit none
      !
      real*8 Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C1,C0,p,q
      !
      real*8 zero,one,two,three,fourth,third,half,pi
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0, pi=3.141592653d0)



      ! Compute L
      !
      Lf = C0*K2**p + C1*(one+K3)*K2**q


      ! Compute first derivative of L wrt K2	
      !
      dLdK2 = C0*p*K2**(p-one) + C1*(one+K3)*q*K2**(q-one)
     
     
      ! Compute first derivative of L wrt K3
      !
      dLdK3 = C1*K2**q 
      
      
      ! Compute second derivative of L wrt K2
      !
      dLdK2dK2 = C0*p*(p-one)*K2**(p-two) + 
     +           C1*(one+K3)*q*(q-one)*K2**(q-two)
      
      
      ! Compute mixed derivative of L wrt K2 and K3
      !
      dLdK2dK3 = C1*q*K2**(q-one)
      
      
      return
      end subroutine L_fun
      
!****************************************************************************
!     The next subroutine calculates various kinematical quantities
!      associated with the deformation gradient
!****************************************************************************

      subroutine skinem(F,R,U,E,istat)
      !
      ! This subroutine performs the right polar decomposition
      !  F = RU of the deformation gradient F into a rotation
      !  R and the right stretch tensor U.  The logarithmic 
      !  strain E = ln(U) is also returned.
      !
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	U(3,3):       the right stretch tensor; output
      !	Uinv(3,3):    the inverse of U
      !	C(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Ueigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of U
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
      implicit none
      !
      integer istat
      !
      real*8 F(3,3),C(3,3),omega(3),Ueigval(3),eigvec(3,3),
     +  U(3,3),E(3,3),Uinv(3,3),R(3,3),detF
     

      !	Store the identity matrix in R, U, and Uinv
      !
      call onem(R)
      call onem(U)
      call onem(Uinv)
      

      ! Store the zero matrix in E
      !
      E = 0.d0
      

      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the 
      !  time increment.
      !
      call mdet(F,detF)
      if (detF.le.0.d0) then
        write(*,'(/5X,A/)') '--problem in kinematics-- the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
      

      ! Calculate the right Cauchy-Green tensor C
      !
      C = matmul(transpose(F),F)
      
 
      ! Calculate the eigenvalues and eigenvectors of C
      !
      call spectral(C,omega,eigvec,istat)
      

      ! Calculate the principal values of U and E
      !
      Ueigval(1) = dsqrt(omega(1))
      Ueigval(2) = dsqrt(omega(2))
      Ueigval(3) = dsqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
      !
      E(1,1) = dlog(Ueigval(1))
      E(2,2) = dlog(Ueigval(2))
      E(3,3) = dlog(Ueigval(3))
      

      ! Calculate the complete tensors U and E
      !
      U = matmul(matmul(eigvec,U),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
      

      ! Calculate Uinv
      !
      call matInv3D(U,Uinv,detF,istat)
      

      ! calculate R
      !
      R = matmul(F,Uinv)
      

      return
      end subroutine skinem

!****************************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!****************************************************************************

      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)


      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      call eigsrt(D,V,3,np)
	

      return
      end subroutine spectral
	
!****************************************************************************

      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
     
      
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
	B(ip) = A(ip,ip)
	D(ip) = B(ip)
	Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the 
              !  off-diagonal element is small.
              !
	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +            .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
	          T =A(ip,iq)/H
	        else
	          theta = 0.5d0*H/A(ip,iq)
	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
	          if (theta.lt.0.d0) T = -T
	        end if
	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
	        S = T*C
	        tau = S/(1.d0 + C)
	        H = T*A(ip,iq)
	        Z(ip) = Z(ip) - H
	        Z(iq) = Z(iq) + H
	        D(ip) = D(ip) - H
	        D(iq) = D(iq) + H
	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
		!		
	        do j=1,ip-1
	          G = A(j,ip)
	          H = A(j,iq)
	          A(j,ip) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations P < J < Q
                !
	        do j=ip+1,iq-1
	          G = A(ip,j)
	          H = A(j,iq)
	          A(ip,j) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations Q < J <= N
                !
	        do j=iq+1,n
                  G = A(ip,j)
	          H = A(iq,j)
	          A(ip,j) = G - S*(H + G*tau)
	          A(iq,j) = H + S*(G - H*tau)
	        end do
	        do j = 1,n
	          G = V(j,ip)
	          H = V(j,iq)
	          V(j,ip) = G - S*(H + G*tau)
	          V(j,iq) = H + S*(G - H*tau)
	        end do
	        nrot = nrot + 1
              end if
	    end do
	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
	  do ip=1,n
	    B(ip) = B(ip) + Z(ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  end do
	end do


      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
      

      return
      end subroutine jacobi
	
!****************************************************************************

      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  jacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      

      do i=1,n-1
	k = i
	P = D(i)
	do j=i+1,n
	  if (D(j).ge.P) then
	    k = j
	    P = D(j)
	  end if
	end do
	if (k.ne.i) then
	  D(k) = D(i)
	  D(i) = P
	  do j=1,n
	    P = V(j,i)
	    V(j,i) = V(j,k)
	    V(j,k) = P
	  end do
  	end if
      end do
      

      return
      end subroutine eigsrt

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: SUBROUTINE matInv3:'
        write(*,*) 'WARNING: DET of MAT=',DET_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************

      subroutine dlnxdx(X,DYDX)
      !
      ! This subroutine calculates the derivative of the logarithm
      ! of a symmetric tensor with respect to that tensor
      !
      implicit none
      !
      integer i,j,k,l
      !
      real*8 X(3,3),DYDX(3,3,3,3),Iden(3,3),Iden4(3,3,3,3),eigval(3),
     +  eigvec(3,3),ehat1(3),ehat2(3),ehat3(3),E1(3,3),E2(3,3),E3(3,3),
     +  y(3),DX2DX(3,3,3,3),s1,s2,s3,s4,s5,s6
      !
      real*8 zero,one,two,half,three,third,small
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     +     three=3.d0,third=1.d0/3.d0,small=1.d-12)
      
      
      ! Initialize
      !
      ! Second order identity tensor
      !
      call onem(Iden)
      !
      ! Fourth order symmetric identity tensor
      !
      Iden4 = zero
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              Iden4(i,j,k,l) = half*(Iden(i,k)*Iden(j,l) + 
     +                                Iden(i,l)*Iden(j,k))
            end do
          end do
        end do
      end do
      
      
      ! Calculate the eigenvalues and eigenvectors of X
      !
      call spectral(X,eigval,eigvec)
      !
      ! Extract the eigenvectors
      !
      do i=1,3
        ehat1(i) =  eigvec(i,1)
        ehat2(i) =  eigvec(i,2)
        ehat3(i) =  eigvec(i,3)
      end do
      !
      ! Assemble the eigenprojections
      !
      do i=1,3
        do j=1,3
	  E1(i,j) = ehat1(i)*ehat1(j)
	  E2(i,j) = ehat2(i)*ehat2(j)
	  E3(i,j) = ehat3(i)*ehat3(j)
	end do
      end do
      
      
      ! Calculate the eigenvalues of Y = ln(X)
      !
      y = dlog(eigval)
      
      
      ! Calculate the derivative of X^2 with respect to X
      !
      DX2DX = zero
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              DX2DX(i,j,k,l) = half*(Iden(i,k)*X(j,l) + 
     +                                Iden(i,l)*X(j,k) + 
     +                                X(i,k)*Iden(j,l) + 
     +                                X(i,l)*Iden(j,k))
            end do
          end do
        end do
      end do
         
            
      ! Calculate DYDX
      !
      DYDX = zero
      if (dabs(eigval(1)-eigval(3)).le.small) then
        !
        ! Three repeated eigenvalues
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = (one/eigval(1))*Iden4(i,j,k,l)
              end do
            end do
          end do
        end do
        !
      elseif (dabs(eigval(2)-eigval(3)).le.small) then
        !
        ! The eigenvalues 2 and 3 are repeated. Eigenvalue 1 is distinct.
        !
        s1 = (y(1) - y(2))/((eigval(1)-eigval(2))**two) - 
     +                  (one/eigval(2))/(eigval(1)-eigval(2))
        s2 = two*eigval(2)*(y(1)-y(2))/((eigval(1)-eigval(2))**two) - 
     +     (one/eigval(2))*(eigval(1)+eigval(2))/(eigval(1)-eigval(2))
        s3 = two*(y(1)-y(2))/((eigval(1)-eigval(2))**three) - 
     +   ((one/eigval(1))+(one/eigval(2)))/((eigval(1)-eigval(2))**two)
        s4 = eigval(2)*s3
        s5 = s4
        s6 = (eigval(2)**two)*s3
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = s1*DX2DX(i,j,k,l) - s2*Iden4(i,j,k,l) - 
     +                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) + 
     +                     s5*Iden(i,j)*X(k,l) - s6*Iden(i,j)*Iden(k,l)
              end do
            end do
          end do
        end do
        !
      elseif (dabs(eigval(1)-eigval(2)).le.small) then
        !
        ! The eigenvalues 1 and 2 are repeated. Eigenvalue 3 is distinct.
        !
        s1 = (y(3) - y(2))/((eigval(3)-eigval(2))**two) - 
     +                  (one/eigval(2))/(eigval(3)-eigval(2))
        s2 = two*eigval(2)*(y(3)-y(2))/((eigval(3)-eigval(2))**two) - 
     +     (one/eigval(2))*(eigval(3)+eigval(2))/(eigval(3)-eigval(2))
        s3 = two*(y(3)-y(2))/((eigval(3)-eigval(2))**three) - 
     +   ((one/eigval(3))+(one/eigval(2)))/((eigval(3)-eigval(2))**two)
        s4 = eigval(2)*s3
        s5 = s4
        s6 = (eigval(2)**two)*s3
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = s1*DX2DX(i,j,k,l) - s2*Iden4(i,j,k,l) - 
     +                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) + 
     +                     s5*Iden(i,j)*X(k,l) - s6*Iden(i,j)*Iden(k,l)
              end do
            end do
          end do
        end do
        !
      else
        !
        ! Eigenvalues are distinct.
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = (y(1)/((eigval(1)-eigval(2))*
     +                                 (eigval(1)-eigval(3))))*
     +         (DX2DX(i,j,k,l) - (eigval(2)+eigval(3))*Iden4(i,j,k,l) - 
     +  ((eigval(1)-eigval(2))+(eigval(1)-eigval(3)))*E1(i,j)*E1(k,l) - 
     +       (eigval(2)-eigval(3))*(E2(i,j)*E2(k,l)-E3(i,j)*E3(k,l))) + 
     +                        (one/eigval(1))*E1(i,j)*E1(k,l) +
     +                          (y(2)/((eigval(2)-eigval(1))*
     +                                 (eigval(2)-eigval(3))))*
     +         (DX2DX(i,j,k,l) - (eigval(1)+eigval(3))*Iden4(i,j,k,l) - 
     +  ((eigval(2)-eigval(1))+(eigval(2)-eigval(3)))*E2(i,j)*E2(k,l) - 
     +       (eigval(1)-eigval(3))*(E1(i,j)*E1(k,l)-E3(i,j)*E3(k,l))) + 
     +                        (one/eigval(2))*E2(i,j)*E2(k,l) +
     +                          (y(3)/((eigval(3)-eigval(1))*
     +                                 (eigval(3)-eigval(2))))*
     +         (DX2DX(i,j,k,l) - (eigval(1)+eigval(2))*Iden4(i,j,k,l) - 
     +  ((eigval(3)-eigval(1))+(eigval(3)-eigval(2)))*E3(i,j)*E3(k,l) - 
     +       (eigval(1)-eigval(2))*(E1(i,j)*E1(k,l)-E2(i,j)*E2(k,l))) + 
     +                        (one/eigval(3))*E3(i,j)*E3(k,l)
              end do
            end do
          end do
        end do
        !
      end if
      
      return
      end subroutine dlnxdx

!****************************************************************************

      subroutine mprod4(A,B,C)
      !
      ! This subroutine calculates the product of two fourth order tensors,
      ! A and B, and places the product in C:
      !             C_{ijkl} = A_{ijmn}B_{mnkl}.
      !
      implicit none
      !
      integer i,j,k,l,m,n
      !
      real*8 A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
      
      
      C = 0.d0
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              do m=1,3
                do n=1,3
                  C(i,j,k,l) = C(i,j,k,l) + A(i,j,m,n)*B(m,n,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
      
      
      return
      end subroutine mprod4

!****************************************************************************