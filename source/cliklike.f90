    module cliklike
    use clik
    use CosmologyTypes
    use CosmoTheory
    use settings
    use Likelihood_Cosmology
    implicit none

    logical :: use_clik = .false.

    integer, parameter :: dp = kind(1.d0)

    type, extends(TCMBLikelihood) :: ClikLikelihood
        type(clik_object) :: clikid
        integer(kind=4),dimension(6) :: clik_has_cl, clik_lmax
        integer :: clik_n,clik_ncl,clik_nnuis
        integer, allocatable :: clik_index(:,:)
        character (len=256), dimension(:), pointer :: names
    contains
    procedure :: LogLike => clik_LnLike
    procedure :: clik_likeinit
    procedure :: do_checks
    end type ClikLikelihood

    type, extends(ClikLikelihood) :: ClikLensingLikelihood
        integer(kind=4), dimension(7) :: lensing_lmaxs
    contains
    procedure :: LogLike => clik_lensing_LnLike
    procedure :: clik_likeinit => clik_lensing_likeinit
    end type ClikLensingLikelihood

    private

    public :: clik_readParams, use_clik

    contains

    subroutine clik_readParams(LikeList,Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) Ini
    character (LEN=:), allocatable ::  name,fname,params
    integer i
    Class(ClikLikelihood), pointer :: like
    logical(KIND=4) is_lensing
    character(LEN=:), allocatable :: tag

    do i=1, Ini%Count
        name = Ini%Items(i)%P%Name
        if (StringStarts(name,'clik_data_')) then
            fname = Ini%ReadFileName(name, NotFoundFail = .false.)
            if (fname=='') cycle
            if (MpiRank==0 .and. feedback > 0) &
                print*,'Using clik with likelihood file ',trim(fname)
            fname = fname // ' ' !just be safe to avoid clik bug
            call clik_try_lensing(is_lensing, fname)
            if (is_lensing) then
                allocate(ClikLensingLikelihood::like)
            else
                allocate(ClikLikelihood::like)
            end if
            call LikeList%Add(Like)
            Like%needs_powerspectra =.true.
            Like%LikelihoodType = 'CMB'
            tag = trim(name(len('clik_data_')+1:))
            Like%Tag = tag
            Like%name= File%ExtractName(fname, no_ext = .true.)
            allocate(Like%clik_index(2,6), source=0)
            Like%clik_index(:,1) = [1,1]
            Like%clik_index(:,2) = [2,2]
            Like%clik_index(:,3) = [3,3]
            Like%clik_index(:,4) = [2,1]
            Like%clik_index(:,5) = [3,1]
            Like%clik_index(:,6) = [3,2]

            params = Ini%ReadFileName('clik_params_'//tag, NotFoundFail = .false.)
            if (params/='') call Like%loadParamNames(params)
            like%speed = Ini%Read_Int('clik_speed_'//tag, 0)
            call Like%clik_likeinit(fname)
        end if
    end do

    end subroutine clik_readParams

    subroutine clik_likeinit(this, fname)
    class (ClikLikelihood) :: this
    character(LEN=*), intent(in) :: fname
    character (len=2),dimension(6) :: clnames
    integer i

    if (Feedback > 1 .and. MPIRank==0) Print*,'Initialising clik...'
    call clik_init(this%clikid,fname)
    call clik_get_has_cl(this%clikid,this%clik_has_cl)
    call clik_get_lmax(this%clikid,this%clik_lmax)

    !output Cls used
    clnames(1)='TT'
    clnames(2)='EE'
    clnames(3)='BB'
    clnames(4)='TE'
    clnames(5)='TB'
    clnames(6)='EB'
    if (Feedback > 1 .and. MPIRank==0) print*,'Likelihood uses the following Cls:'
    do i=1,6
        if (this%clik_has_cl(i) .eq. 1) then
            print*,'  ',trim(clnames(i)),' from l=0 to l=',this%clik_lmax(i)
        end if
    end do

    this%clik_ncl = sum(this%clik_lmax) + 6

    allocate(this%cl_lmax(CL_B,CL_B), source=0)
    this%cl_lmax(CL_T,CL_T) = this%clik_lmax(1)
    this%cl_lmax(CL_E,CL_T) = this%clik_lmax(4)
    this%cl_lmax(CL_E,CL_E) = this%clik_lmax(2)
    this%cl_lmax(CL_B,CL_B) = this%clik_lmax(3)
    where (this%cl_lmax<0)
        this%cl_lmax=0
    end where

    this%clik_nnuis = clik_get_extra_parameter_names(this%clikid,this%names)
    call this%do_checks()

    !tidying up

    this%clik_n = this%clik_ncl + this%clik_nnuis

    end subroutine clik_likeinit

    real(mcp) function clik_lnlike(this, CMB, Theory, DataParams)
    Class(ClikLikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    integer :: i,j ,l
    real(dp) clik_cl_and_pars(this%clik_n)

    !set C_l and parameter vector to zero initially
    clik_cl_and_pars = 0.d0

    j = 1

    !TB and EB assumed to be zero
    !If your model predicts otherwise, this function will need to be updated
    do i=1,4
        associate(Cls=> Theory%Cls(this%clik_index(1,i),this%clik_index(2,i)))
            do l=0,this%clik_lmax(i)
                !skip C_0 and C_1
                if (l >= 2) then
                    clik_cl_and_pars(j) = CLs%CL(L)/real(l*(l+1),mcp)*twopi
                end if
                j = j+1
            end do
        end associate
    end do
    j =j + sum(this%clik_lmax(5:6)+1) !zero arrays of TB, EB

    !Appending nuisance parameters
    if (this%clik_nnuis > 0) then
        do i=1,this%clik_nnuis
            clik_cl_and_pars(j) = DataParams(i)
            j = j+1
        end do
    end if

    !Get - ln this needed by CosmoMC
    clik_lnlike = -1.d0*clik_compute(this%clikid,clik_cl_and_pars)

    if (Feedback>1) Print*,trim(this%name)//' lnlike = ',clik_lnlike

    end function clik_lnlike

    subroutine do_checks(this)
    class (ClikLikelihood) :: this
    integer i

    !Safeguard
    if (this%clik_nnuis/= this%nuisance_params%nnames) call MpiStop(FormatString( &
        'clik_nnuis (%u) has different number of nuisance parameters than .paramnames  (%u)',this%clik_nnuis,this%nuisance_params%nnames))
    if (this%clik_nnuis /= 0 .and. MPIRank==0 .and. Feedback>0) then
        Print*,'Clik will run with the following nuisance parameters:'
        do i=1,this%clik_nnuis
            Print*,trim(this%names(i))
        end do
    end if
    if (this%clik_nnuis >0) deallocate(this%names)

    end subroutine do_checks

    subroutine clik_lensing_likeinit(this, fname)
    class (ClikLensingLikelihood) :: this
    character(LEN=*), intent(in) :: fname

    if (Feedback > 1) Print*,'Initialising clik lensing...'
    call clik_lensing_init(this%clikid,fname)
    !    call clik_lensing_get_lmax(this%clikid,this%lensing_lmax)

    ! fill the lmaxs array. It contains (in that order) lmax for cl_phiphi, cl_TT, cl_EE, cl_BB, cl_TE, cl_TB, cl_EB
    ! -1 means that a particular spectrum is not used.
    ! lmaxs(1) will never be -1, but all the other can be is the likelihood is set not to include renormalization
    call clik_lensing_get_lmaxs(this%clikid,this%lensing_lmaxs)

    this%clik_nnuis = clik_lensing_get_extra_parameter_names(this%clikid,this%names)
    allocate(this%cl_lmax(CL_phi,CL_phi), source=0)
    this%cl_lmax(CL_T,CL_T) = this%lensing_lmaxs(2)
    this%cl_lmax(CL_E,CL_T) = this%lensing_lmaxs(5)
    this%cl_lmax(CL_E,CL_E) = this%lensing_lmaxs(3)
    this%cl_lmax(CL_B,CL_B) = this%lensing_lmaxs(4)
    this%cl_lmax(CL_Phi,CL_Phi) = this%lensing_lmaxs(1)
    where (this%cl_lmax<0)
        this%cl_lmax=0
    end where
    call this%do_checks()
    print *,'lensing lmax: ', this%lensing_lmaxs
    this%clik_n = sum(this%lensing_lmaxs(1:7)+1) + this%clik_nnuis

    end subroutine clik_lensing_likeinit

    real(mcp) function clik_lensing_lnlike(this, CMB, Theory, DataParams)
    Class(ClikLensingLikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    integer :: i,j,l
    real(dp) clik_cl_and_pars(this%clik_n)


    !set C_l and parameter vector to zero initially
    clik_cl_and_pars = 0.d0

    j = 1
    associate(Cls=> Theory%Cls(CL_Phi,CL_Phi))
        do l=0,this%lensing_lmaxs(1)
            !skip C_0 and C_1
            if (l >= 2) then
                clik_cl_and_pars(j) = CLs%Cl(L)/real(l*(l+1),mcp)**2*twopi
            end if
            j = j+1
        end do
    end associate

    !TB and EB assumed to be zero
    !If your model predicts otherwise, this function will need to be updated
    do i=1,4
        associate(Cls=> Theory%Cls(this%clik_index(1,i),this%clik_index(2,i)))
            do l=0,this%lensing_lmaxs(i+1)
                !skip C_0 and C_1
                if (l >= 2) then
                    clik_cl_and_pars(j) = CLs%Cl(L)/real(l*(l+1),mcp)*twopi
                end if
                j = j+1
            end do
        end associate
    end do

    do i=1,this%clik_nnuis
        clik_cl_and_pars(j) = DataParams(i)
        j = j+1
    end do

    !Get - ln this needed by CosmoMC
    clik_lensing_lnlike = -1.d0*clik_lensing_compute(this%clikid,clik_cl_and_pars)

    if (Feedback>1) Print*,trim(this%name)//' lnlike = ',clik_lensing_lnlike

    end function clik_lensing_lnlike

    end module cliklike
