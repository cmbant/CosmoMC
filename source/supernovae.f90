!--------------------------------------------------------------------!
! This can be used to replace the supernovae module in cosmomc.      !
!--------------------------------------------------------------------!
!                                                                    !
! Original module modified by Yun Wang (assisted by Pia Mukherjee)   ! 
! to implement flux-averaging for SNe with correlated errors.        !
! SN data files from the SNLS data arXiv are supplied too.           !
!                                                                    !
!                                                                    !
! For details of the flux-averaging method, see arXiv:1109.3172,     !
! Yun Wang, Chia-Hsun Chuang, Pia Mukherjee, "A Comparative Study    !
! of Dark Energy Constraints from Current Observational Data".       !
! See also http://www.nhn.ou.edu/~wang/SNcode/index.html             !
!           							     !
! -------------------------------------------------------------------!
!                                                                    !
! This module uses the compilation of SNe by Conley et al. (2011)    !
! [Conley, A., et al., 2011, Astrophys.J.Suppl., 192, 1]             !
!                                                                    !
!--------------------------------------------------------------------!
!								     !
! Note: Marginalize over the SN nuisance parameters 		     !
!	(alpha, beta, Msn) by varying them.                          !
!                                                                    !	
! To modify the code for arbitrary H(z): replace function "invEz(z)" !
! with H0/H(z) you wish to study [H0 is the Hubble constant].        !
!                                                                    !
!--------------------------------------------------------------------!
! Yun Wang, September 16, 2011                                       !
!--------------------------------------------------------------------!
! Yun Wang, Oct 18, 2011: added clarification on  (alpha, beta, Msn) !                                !
!--------------------------------------------------------------------!

module snovae
use cmbtypes
use MatrixUtils
implicit none

 integer, parameter :: SN_num = 472
 double precision, parameter :: Pi_num = 3.14159265359D0 
 double precision :: SN_z(SN_num), SN_moduli(SN_num)
 double precision :: SN_dz(SN_num), SN_dmoduli(SN_num)
 double precision :: SN_stretch(SN_num), SN_dstretch(SN_num)
 double precision :: SN_color(SN_num), SN_dcolor(SN_num)
 double precision :: SN_cov_ms(SN_num), SN_cov_mc(SN_num), SN_cov_sc(SN_num)
 double precision :: SN_dii0(SN_num)
 double precision :: SN_Ninv(SN_num,SN_Num)
 double precision :: SN_sumninv,dum
 double precision :: alpha,beta,Msn,pz
 double precision :: v0(SN_num,SN_num),v0a(SN_num,SN_num),v0b(SN_num,SN_num)
 double precision :: va(SN_num,SN_num),vb(SN_num,SN_num),vab(SN_num,SN_num)
 real :: Omegak,Omegam,w0,wa
 integer :: nmin
 logical, parameter :: SN_marg = .True.
 logical, parameter :: SN_syscovmat = .True.  !! Use covariance matrix with or without systematics
 

contains


 subroutine SN_init
   use settings
   character (LEN=12):: name
   integer i, set
   real :: siginta(0:3),zmin,zmax

   siginta(0)=0.068
   siginta(1)=0.113
   siginta(2)=0.082
   siginta(3)=0.099

   alpha=1.2
   beta=3.64
   Msn=23.94

   if (Feedback > 0) write (*,*) 'Reading: snls3 supernovae data'
   call OpenTxtFile(trim(DataDir)//'snls3_combined.dat',tmp_file_unit)
   do i=1,  sn_num
      read(tmp_file_unit, *) name, SN_z(i), dum, SN_dz(i), SN_moduli(i), &
           SN_dmoduli(i), SN_stretch(i), SN_dstretch(i), SN_color(i), &
           SN_dcolor(i), dum, dum, SN_cov_ms(i), SN_cov_mc(i), SN_cov_sc(i), set

   pz=0.0005 ! 150km/s
   SN_dz(i)=sqrt(SN_dz(i)**2+pz*pz) ! adding peculiar velocity residuals
   SN_dii0(i)=SN_dmoduli(i)**2+siginta(set)**2

   if (SN_z(i).lt.zmin) zmin=SN_z(i)
   if (SN_z(i).gt.zmax) zmax=SN_z(i)	 

   end do
   close(tmp_file_unit)

   call OpenTxtFile(trim(DataDir)//'snls3_v0_covmatrix.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum
   do i=1, sn_num
      read (tmp_file_unit,*) v0(i,1:sn_num)
   end do
   close (tmp_file_unit)


   call OpenTxtFile(trim(DataDir)//'snls3_v0a_covmatrix.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum
   do i=1, sn_num
      read (tmp_file_unit,*) v0a(i,1:sn_num)
   end do
   close (tmp_file_unit)

   call OpenTxtFile(trim(DataDir)//'snls3_v0b_covmatrix.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum
   do i=1, sn_num
      read (tmp_file_unit,*) v0b(i,1:sn_num)
   end do
   close (tmp_file_unit)

   call OpenTxtFile(trim(DataDir)//'snls3_va_covmatrix.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum
   do i=1, sn_num
      read (tmp_file_unit,*) va(i,1:sn_num)
   end do
   close (tmp_file_unit)

   call OpenTxtFile(trim(DataDir)//'snls3_vb_covmatrix.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum
   do i=1, sn_num
      read (tmp_file_unit,*) vb(i,1:sn_num)
   end do
   close (tmp_file_unit)

   call OpenTxtFile(trim(DataDir)//'snls3_vab_covmatrix.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum
   do i=1, sn_num
      read (tmp_file_unit,*) vab(i,1:sn_num)
   end do
   close (tmp_file_unit)

   
  end subroutine SN_init


 function SN_LnLike(CMB)
   use camb
  !Assume this is called just after CAMB with the correct model  use camb
  implicit none
  type(CMBParams) CMB
  logical, save :: do_SN_init = .true.
  integer, parameter :: nq=20
  integer, parameter :: np=1000

  real SN_LnLike
  integer i,j,nmin,nbin,snnum,q,count,counta(np),indx(np),i1,i2
  integer SNa_nzlow(np),SNA_nzhigh(np)
  real z
  real diffs(SN_num), chisq
  real  zlowa(np), zhigha(np), dzhighmin(np), dzlowmin(np),za(np),zabin(np)
  real dzlow,dzhigh,zlowprev(np),zhighprev(np),H0dLa(np),H0dl,xsn(np),ysn(np),drdz
  real lumd(np),fluxd(np),rcov0,rcov0a(np),Dii(np),sigmz2,H0dld
  real suml,zmean,facz
  real faczmean, sum0,fluxbin(np),y0a(np)
  real dzfluxbin,rcov0bina(np),invcovmatF(np,np)
  real eps,kappa,Gamma,del
  real, allocatable :: nn0(:,:)
  double precision covmat(SN_num,SN_num),covmatF(np,np),sum

  write(*,*) 'in SN_LnLike'

  if (do_SN_init) then 
     call SN_init
     do_SN_init = .false.
  end if


  Omegak=CMB%omk
  Omegam=CMB%omb+CMB%omdm
  w0=CMB%w
  wa=CMB%nufrac

! test values 
!  Omegak=0.
!  Omegam=0.21
!  w0=-1.
!  wa=0.

  nmin=10
  eps=1.e-6
  kappa=sqrt(abs(Omegak))
  dzfluxbin=0.07

! converting to flux
     do i=1, SN_num
        z= SN_z(i)

        call intg(0.,z,invEz,eps,Gamma,del)  
        if (Omegak.lt.0.0) then ! closed universe
           rcov0 = sin(kappa*Gamma)/kappa
        else if (Omegak.gt.0.) then ! open universe
           rcov0 = sinh(kappa*Gamma)/kappa
        else 
           rcov0=Gamma	         
        end if
        H0dL=rcov0*(1.0+z) !H0*dL	   
        rcov0a(i)=rcov0

        if (Omegak.lt.0.0) then ! closed universe
           drdz = cos(kappa*Gamma)*invEz(z)
        else if (Omegak.gt.0.) then ! open universe
           drdz = cosh(kappa*Gamma)*invEz(z)
        else 
           drdz=invEz(z) !flat universe      
        end if
	      	      
        sigmz2=((5./log(10.))*(1./(1.+z)+drdz/rcov0))**2 * &
             SN_dz(i)**2      !c=3.e5km/s   
	      
        Dii(i) =SN_dii0(i)+ (alpha*SN_dstretch(i))**2+ (beta*SN_dcolor(i))**2 +&  
             2.*alpha*SN_cov_ms(i)- 2.* beta*SN_cov_mc(i) - &
             2.*alpha*beta*SN_cov_sc(i) + sigmz2  

        H0dLd= log(10.)*(SN_moduli(i)-(Msn-&
             alpha*(SN_stretch(i)-1.0)+beta*SN_color(i)) )/5.
        H0dld=ex(H0dld)
        fluxd(i)=1./H0dLd**2
        lumd(i)=fluxd(i)*H0dL**2
     end do

! sum up snls3 cov matrix
     do i=1,SN_num
        do j=i,SN_num            
           covmat(i,j)=v0(i,j)+alpha**2*va(i,j) + beta**2*vb(i,j) + &
                2.*alpha*v0a(i,j) - 2.*beta*v0b(i,j) - 2.*alpha*beta*vab(i,j)
           if (i.eq.j) then
              covmat(i,j)=covmat(i,j)+Dii(i)
           end if
           covmat(j,i)=covmat(i,j)
        end do
     end do

! binning
     za(1)=0.0
     do i=1,nq
        za(i+1)=za(1)+ dzfluxbin*real(i)
        zlowprev(i+1)=za(i+2)
     end do

     nbin=0       
     snnum=0   
        	  
     do q=1,nq
        sumL=0.
        zmean=0.
        count=0
        
        dzlowmin(q)=0.1
        dzhighmin(q)=0.1     	    
        
        do j=1,SN_num
           z=SN_z(j)
           
           if (z.gt.za(q).and.z.le.za(q+1)) then
              count=count+1
              snnum=snnum+1
              rcov0=rcov0a(j)    
              facz= ((1.+z)*rcov0)**2 
              sumL=sumL+ fluxd(j)*facz
              zmean=zmean+z
              dzlow=abs(z-za(q)) 	
              if (dzlow.le.dzlowmin(q)) then
                 zlowa(q)=z	    
                 SNa_nzlow(q)=snnum
                 dzlowmin(q)=dzlow	
                 zlowprev(q)=z	    
              end if
              dzhigh=abs(z-za(q+1)) 	
              if (dzhigh.le.dzhighmin(q).and.z.ge.zhighprev(q)) then
                 zhigha(q)=z	    
                 SNa_nzhigh(q)=snnum
                 dzhighmin(q)=dzhigh
                 zhighprev(q)=z		    
              end if
           end if
        end do
        
        SNa_nzlow(1)=1
        if (count.gt.0) then
           if (q.gt.1.and.SNa_nzlow(q).gt.(SNa_nzhigh(q-1)+1).and.&
                SNa_nzhigh(q-1).gt.0) then 
              SNa_nzlow(q)=SNa_nzhigh(q-1)+1
           end if
           counta(q)=count
        end if
        SNa_nzlow(1)=1
        if (count.gt.0) then
           if (q.gt.1.and.SNa_nzlow(q).gt.(SNa_nzhigh(q-1)+1).and.&
                SNa_nzhigh(q-1).gt.0) then 
              SNa_nzlow(q)=SNa_nzhigh(q-1)+1
           end if
           counta(q)=count
        end if

        if (count.gt.0) then
           nbin=nbin+1
           zmean=zmean/real(count)
           call intg(0.,zmean,invEz,eps,Gamma,del)  
           if (Omegak.lt.0.0) then ! closed universe
              rcov0 = sin(kappa*Gamma)/kappa
           else if (Omegak.gt.0.) then ! open universe
              rcov0 = sinh(kappa*Gamma)/kappa
           else 
              rcov0=Gamma	         
           end if
           faczmean= ((1.+zmean)*rcov0)**2  		   
           sum0= sumL/( faczmean*dble(count) )
           rcov0bina(nbin)=rcov0      
           zabin(nbin)=zmean
           fluxbin(nbin)=sum0
        end if
             
     end do	!end binning

! compute covariance matrix for correlated errors in SNe:	

     do i=1,nbin
        xsn(i)=zabin(i)
        z=xsn(i)
        rcov0=rcov0bina(i)
        H0dLa(i)=(1.+z)*rcov0	
        y0a(i)=1./H0dLa(i)**2 ! predicted "flux"		  
        ysn(i)=fluxbin(i)
     end do
	
     do i=1,nbin
        do j=1,nbin
           sum=0.	  	    
           do i1=SNa_nzlow(i), SNa_nzhigh(i) !1, counta(j)
	      do i2=SNa_nzlow(j), SNa_nzhigh(j)
                 sum=sum+lumd(i1)*lumd(i2)*covmat(i1,i2)
	      end do
           end do
           
           covmatF(i,j)=sum*(log(10.)/(2.5*H0dLa(i)*H0dLa(j)))**2 /&
                (counta(i)*counta(j))
        end do
     end do
	
     allocate(nn0(nbin,nbin))


     do i=1,nbin
        nn0(i,1:nbin)=covmatF(i,1:nbin)
     end do

! invert it:	
     call Matrix_Inverse(nn0)
     
     chisq=0.
     do i=1,nbin
        do j=1,nbin
           chisq=chisq+(ysn(i)-y0a(i))*nn0(i,j)*(ysn(j)-y0a(j))
        end do
     end do


     SN_LnLike = chisq/2

     deallocate(nn0)

     if(Feedback>0) write(*,*) 'SN_LnLike from flux averaged SNLS3 SNe', SN_LnLike

 end function SN_LnLike

 subroutine intg (a,b,f,eps,ans,del)
   implicit real (a-h,o-z)
   integer nmin,n,m,k,i,j
   external f
   
   dimension t(20,20)	!	matrix of approximations
   c = 0.50 * (b + a)	!	center
   d = 0.50 * (b - a)	!	width
   n = 2			!	level
   m = 1			!	counts something or other
   e = 1.0		!	something like a step size
   t(1,1) = 0.0		!	(cf. HP comment above)
   t(1,2) = 2.0 * d * f(c)
   t(2,1) = 0.750 * t(1,2)
   
  do
     n = n + 1
     m = m * 2
     e = e * 0.50
     s = 0.0
     do j = 2, m, 2
        y = real(j-1) * e
        x = 0.50*y*(3.0 - y**2)
        s = s + (1.0 - y**2) * (f(c-d*x) + f(c+d*x))
     end do
     t(n,1) = 1.50*s*d*e + 0.50*t(n-1,1)
   
     p = 1.0
     do k = 1, n-1
        p = p * 4.0
        i = n + 1 - k
        t(i-1,k+1) = t(i,k) + (t(i,k) - t(i-1,k))/(p-1.0)
     end do
   
     ans = t(1,n)
     del = abs(t(1,n)-t(2,n-1))
!     if (n .lt. nmin) go to 1000
     if (n .ge. 20) return
     if (abs(del) .le. eps*abs(ans)) return
  end do
         
  return
 end subroutine intg

!------------------------------------------------------------------

	function invEz(z)
	implicit none
	real sum0, invEz, OmegaX, Xz,z,w
							
	w=1.+z

	OmegaX=1.-Omegam-Omegak
	Xz= ex(3.*(1.+w0+wa)*log(w)+ 3.*wa*(1./w-1.))	!w0+wa(1-a)
	sum0=Omegam*w**3+ Omegak*w**2+ OmegaX*Xz

	invEz = 1.0/sqrt(sum0)
	
	return
	end function invEz

	function ex(x)
	implicit none
	real x, ex
	
	if (x.lt.-15.) then
	  ex=0.
	else
	  ex=exp(x)
	end if      

	return
	end function ex

end module snovae
