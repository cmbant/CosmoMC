!Estimate the covariance matrix by making steps on a grid in nd
! slb 5 aug 2004

! Overwrites whatever is in propose_matrix with the newly estimated matrix
! (In future could add/replace some entries with any better info on the covmat
! which is already known.)
! Uses the steps from the params.ini file to choose grid

module EstCovmatModule
 use Random
 use ParamDef
 use CalcLike
 use settings
 use Matrixutils
 implicit none

 real, dimension(:,:,:,:), allocatable :: Lgrid
 ! Store likelihood values in LGrid to save on computations
 ! This isn't really used much at the moment though..
 ! Ideally would have an array with n dimensions but too hard to code.
 ! This method has a lot of redundancy, but its the best I can think of.
 ! For each pair of parameters i,j we need to consider 9 likelihood values:
 ! all combinations of stepping forward, backward and staying still in each parameter
 ! These are stored in Lgrid(i,j,istep,jstep)
 ! where istep is 1, 2 or 3 for backwards step, staying still and forwards step
 ! similarly jstep.
 ! This means there are the following redundancies:
 ! Lgrid(i,j,2,2) is the center of the grid
 ! So Lgrid(:,:,2,2) all have the same value. So only use and fill Lgrid(1,1,2,2).
 ! Lgrid(i,j,1,2) means that parameter j has stayed still (jstep=2)
 ! so Lgrid(i,:,1,2) all have the same value. So only use and fill Lgrid(i,1,1,2)
 ! Similarly Lgrid(:,j,2,1) all have the same value. So only use and fill Lgrid(1,j,2,1)
 ! Similarly for Lgrid(i,:,3,2). So only use and fill Lgrid(i,1,3,2).
 ! Similarly for Lgrid(:,j,2,3). So only use and fill Lgrid(1,j,2,3).
 ! Finally Lgrid(i,j,istep,jstep)= Lgrid(j,i,jstep,istep). 
 ! So in practice only fill and use Lgrid for j>i 

  
contains


 function EstCovmat(Params, bestdll_in, sucess)
   integer :: ntries = 10 ! Maximum number of attempts to get a covmat with all +ve evals
   integer itries ! Number of tries so far
   integer, intent(out) :: sucess
   real parsteps(num_params_used), bestdll
   real, intent(in) :: bestdll_in  ! the dloglike value used to determine grid spacing eg. try 4
   Type(ParamSet) Params, CenterParams
   !double precision Gaussian1
   real Hess(num_params_used,num_params_used), U(num_params_used,num_params_used)
   real D(num_params_used),EstCovmat(num_params_used,num_params_used)

   itries=0
   sucess=0
   bestdll = bestdll_in


   if (Feedback>1) write(*,*)
   if (Feedback>1) write(*,*) 'Estimating the covariance matrix by finding likelihood values on a grid'

   ! Set up the grid position
   CenterParams = Params 
   allocate(LGrid(num_params_used,num_params_used,3,3))
   Lgrid=0.0 ! Clear it to be on the safe side
   parsteps=Scales%PWidth(params_used)
   call PlaceGrid(CenterParams,parsteps,bestdll)
  
   ! Loop over attempts to find the covariance matrix
   do while (sucess.eq.0)
      itries = itries + 1
      call GetHess(CenterParams,parsteps,Hess)

      U=Hess
      call Matrix_Diagonalize(U, D, num_params_used)  
      if (any(D <= 0)) then
        write(*,*) 
        write(*,*) '!!! -ve evals, so spreading the grid out by a factor of ~2'
        bestdll=bestdll*2 ! try using a more spread out grid 
        call PlaceGrid(CenterParams,parsteps,bestdll)
      else
             sucess=1
      end if

      ! just give up if we have overrun the maximum no. of tries specified above
      if (itries > ntries) then
             write(*,*) 'EstCovmat.f90: Couldn''t find covmat with all positive evals after ', &
                         ntries,' tries - check starting ranges'
         stop
      end if

   end do

   ! uncomment this bit to update the propose_matrix..
    if (sucess>0) then
      U=Hess
      call Matrix_Inverse(U) ! gives the covariance matrix
      EstCovmat(1:num_params_used,1:num_params_used)=U
    end if

   call AcceptReject(.true.,CenterParams%Info,Params%Info)

   deallocate(LGrid)

 end function EstCovmat

 
 subroutine PlaceGrid(CenterParams,parsteps,bestdll)
   ! Takes the input CenterParams and parsteps values and improves on them
   ! While doing so, fills the likelhood grid LGrid
   ! Leaving LGrid(:,:,1:3,2) and LGrid(:,:<2,1:3) set correctly for the returned 
   ! CenterParams and parsteps values
   Type(ParamSet) P, CenterParams
   real parsteps(num_params_used),bestdll
   integer maxstepsizetries

   if (Feedback>1) write(*,*) 'Improving the grid spacings'
   P = CenterParams
   Lgrid(1,1,2,2) = GetLogLike(CenterParams)
   call AcceptReject(.true., P%Info, CenterParams%Info)
   if (Lgrid(1,1,2,2).eq.LogZero) then
         write(*,*) 'ERROR: Trial parameters excluded by prior or bad likelihood'
         write(*,*) 'Try starting further away from problem regions?'
         call DoAbort
   end if

   ! Decide how big to make the grid
   maxstepsizetries=20
   call GetStepsForDChisq1(CenterParams,maxstepsizetries,parsteps,bestdll)
   ! See if this is a good place to center it
   !if (.any.istepsizetries>maxstepsizetries) then it failed...

 end subroutine PlaceGrid


 subroutine GetStepsForDChisq1(CenterParams,maxstepsizetries,parsteps,bestdll)
 !!New version using binary division, AL 22Aug04
   ! Stays centered on CenterParams
   ! But for each paramter, iterates for a maximum of maxstepsizetries to try to find
   ! parsteps such that Delta chisq is of order 1 (when stepping only in that parameter)
   ! Inputs: CenterParams, LCenter, maxstepsizetries
   ! Outputs: parsteps
   Type(ParamSet) CenterParams
   integer, intent(in) :: maxstepsizetries
   real parsteps(num_params_used), bestdll
   integer i, ii, tries
   Type(ParamSet) StepParams
   real LCenter, LStep,dloglike
   real step_min, step_max
   integer asign

   ! The center of the whole grid 
   LCenter = Lgrid(1,1,2,2) ! see notes at top of this file on Lgrid indices
 
   !! Estimate some good step sizes to use by interating until dchisq~1
   do i=1, num_params_used
      ii= params_used(i) 
      
      step_max = 0.
      step_min = 0.
      asign = 1
      if (.not. CenterParams%P(ii) == Scales%PMin(ii) .and. &
          CenterParams%P(ii) > Scales%PMax(ii) - Scales%PWidth(ii)/100) asign = -1

      ! Find step sizes that give a half sensible Delta chisq eg. 0.5 to 10
      do tries = 1, maxstepsizetries

         StepParams=CenterParams ! reset back to original starting point
         StepParams%P(ii)=CenterParams%P(ii) + parsteps(i)*asign
         LStep=GetLogLike(StepParams)
         call AcceptReject(.true., StepParams%Info, CenterParams%Info) 

         dloglike=abs(LStep-LCenter)
         if (dloglike > 3*bestdll) then
           step_max = parsteps(i)
           parsteps(i) = (parsteps(i) + step_min) / 2 
         elseif (dloglike < bestdll/3) then
           step_min = parsteps(i)
           if (step_max==0.) then
            parsteps(i) = parsteps(i) * 2 
           else
            parsteps(i) = (parsteps(i) + step_max) / 2 
           end if
         else
          exit
         end if       
     
         ! just give up if we have overrun the maximum no. of tries specified above
         if (tries == maxstepsizetries) then
            write(*,*) 'Couldn''t find a good stepsize for parameter ',ii,' with Delta chisq ~1 but continuing anyway'
         end if
      end do
      if (Feedback>1) write(*,*) ' Decided on stepsize ',parsteps(i),' for parameter ',ii

   end do

 end subroutine GetStepsForDChisq1


   
 subroutine GetHess(CenterParams,parsteps,Hess)
   ! Get the curvature matrix (Hess) for the grid defined by CenterParams
   ! and parsteps. 
   ! Return the answer in Hess.
   ! Input only: CenterParams, parsteps
   ! Ouput only: Hess
   Type(ParamSet) CenterParams,StepParams
   real parsteps(num_params_used)
   real Hess(num_params_used,num_params_used)
   real Steps(3),wii,wjj
   integer i,j,ii,jj,istep,jstep
       
   Steps(1)=-1 ! Some bits of the code assume these values are -1, 0 1
   Steps(2)=0  ! so can't go changing them without thinking carefully.
   Steps(3)=1  ! The only point of the array is to make loops possible.

   ! Find likelihood values on a grid to estimate curvature matrix (Hess)
   if (Feedback>1) write(*,*) ' Finding curvature matrix ...'
   do i=1, num_params_used
      ii= params_used(i) 
      if (Feedback>1) write(*,*)
      if (Feedback>1) write(*,*) ' Parameter ',i

      wii=parsteps(i)
      do istep=1,3
        jstep=2
        StepParams=CenterParams ! rest back to original starting point
        StepParams%P(ii)=CenterParams%P(ii) + Steps(istep)*wii
        Lgrid(i,1,istep,jstep) = GetLogLike(StepParams)
        call AcceptReject(.true., StepParams%Info, CenterParams%Info)
        if (Lgrid(i,1,istep,jstep).eq.LogZero) then
          write(*,*) 'ERROR: Trial parameters hit prior or error in function evaluation'
          write(*,*) 'Try starting further away from problem regions?'
          call DoAbort
        end if
      end do
      Hess(i,i)=1/(wii*wii) *(Lgrid(i,1,1,2) + Lgrid(i,1,3,2) - 2*Lgrid(i,1,2,2))

       do j=(i+1), num_params_used 
         jj=params_used(j)
         if (Feedback>1) write(*,*) ' Parameter pair ',i,' , ',j
         wjj=parsteps(j)

         ! Complete overkill here, make the whole matrix for debugging purposes
         do istep=1,3
             do jstep= 1,3,2  ! ie just do jstep=1 and jstep=3
               StepParams=CenterParams ! rest back to original starting point
               StepParams%P(jj)=CenterParams%P(jj) + Steps(jstep)*wjj
               StepParams%P(ii)=CenterParams%P(ii) + Steps(istep)*wii
               Lgrid(i,j,istep,jstep) = GetLogLike(StepParams)
               call AcceptReject(.true., StepParams%Info, CenterParams%Info)
               if (Lgrid(i,j,istep,jstep).eq.LogZero) then
                write(*,*) 'ERROR: Trial parameters hit prior or error in function evaluation'
                     write(*,*) 'Try starting further away from problem regions?'
                call DoAbort
               end if

            end do
         end do
!                write(4,*) CenterParams%P(ii),CenterParams%P(jj),wii,wjj,Lgrid

        Hess(j,j)=1/(wjj*wjj) * (Lgrid(i,j,2,1) + Lgrid(i,j,2,3) - 2*Lgrid(i,1,2,2))
        Hess(i,j)=1/(wii*wjj*4)*(Lgrid(i,j,1,1) + Lgrid(i,j,3,3) - Lgrid(i,j,1,3) - Lgrid(i,j,3,1))
        Hess(j,i)=Hess(i,j)

      end do
   end do
               
 end subroutine GetHess

end module EstCovmatModule
