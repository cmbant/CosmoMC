    module cliklike
    use clik
    use cmbtypes
    use CosmoTheory
    use settings
    use Likelihood
    implicit none

    logical :: use_clik = .false.

    integer, parameter :: dp = kind(1.d0)

    type, extends(CosmologyLikelihood) :: ClikLikelihood
        type(clik_object) :: clikid
        integer(kind=4),dimension(6) :: clik_has_cl, clik_lmax
        integer :: clik_n,clik_ncl,clik_nnuis
        character (len=256), dimension(:), pointer :: names
    contains
    procedure :: LogLike => clik_LnLike
    procedure :: clik_likeinit
    procedure :: do_checks
    end type ClikLikelihood

    type, extends(ClikLikelihood) :: ClikLensingLikelihood
        integer(kind=4) lensing_lmax
    contains
    procedure :: LogLike => clik_lensing_LnLike
    procedure :: clik_likeinit => clik_lensing_likeinit
    end type ClikLensingLikelihood

    private
    integer, dimension(4) :: mapped_index

    public :: clik_readParams, use_clik

    contains

    subroutine clik_readParams(LikeList,Ini)
    class(LikelihoodList) :: LikeList
    class(TSettingIni) Ini
    character (LEN=Ini_max_string_len) ::  name,fname,params
    integer i
    Class(ClikLikelihood), pointer :: like
    logical(KIND=4) is_lensing

    do i=1, Ini%Count
        name = Ini%Items(i)%P%Name
        if (name(1:10)=='clik_data_') then
            fname = Ini%ReadFileName(name, NotFoundFail = .false.)
            if (fname=='') cycle
            if (MpiRank==0 .and. feedback > 0) &
            print*,'Using clik with likelihood file ',trim(fname)
            call clik_try_lensing(is_lensing, fname)
            if (is_lensing) then
                allocate(ClikLensingLikelihood::like)
            else
                allocate(ClikLikelihood::like)
            end if
            call LikeList%Add(Like)
            Like%needs_powerspectra =.true.
            Like%LikelihoodType = 'CMB'
            Like%name= ExtractFileName(fname)
            !                Like%version = CAMSpec_like_version
            call StringReplace('clik_data_','clik_params_',name)
            params = Ini%ReadFileName(name, NotFoundFail = .false.)
            if (params/='') call Like%loadParamNames(params)
            call StringReplace('clik_params_','clik_speed_',name)
            like%speed = Ini%Read_Int(name, 0)
            call Like%clik_likeinit(fname)
        end if
    end do

    !Mapping CosmoMC's power spectrum indices to clik's
    mapped_index(1) = 1
    mapped_index(2) = 3
    mapped_index(3) = 4
    mapped_index(4) = 2

    end subroutine clik_readParams

    subroutine clik_likeinit(like, fname)
    class (ClikLikelihood) :: like
    character(LEN=*), intent(in) :: fname
    character (len=2),dimension(6) :: clnames
    integer i

    if (Feedback > 1 .and. MPIRank==0) Print*,'Initialising clik...'
    call clik_init(like%clikid,fname)
    call clik_get_has_cl(like%clikid,like%clik_has_cl)
    call clik_get_lmax(like%clikid,like%clik_lmax)

    !output Cls used
    clnames(1)='TT'
    clnames(2)='EE'
    clnames(3)='BB'
    clnames(4)='TE'
    clnames(5)='TB'
    clnames(6)='EB'
    if (Feedback > 1 .and. MPIRank==0) print*,'Likelihood uses the following Cls:'
    do i=1,6
        if (like%clik_has_cl(i) .eq. 1) then
            print*,'  ',trim(clnames(i)),' from l=0 to l=',like%clik_lmax(i)
        end if
    end do

    like%clik_ncl = sum(like%clik_lmax) + 6

    like%clik_nnuis = clik_get_extra_parameter_names(like%clikid,like%names)
    call like%do_checks()

    !tidying up

    like%clik_n = like%clik_ncl + like%clik_nnuis

    end subroutine clik_likeinit

    real(mcp) function clik_lnlike(like, CMB, Theory, DataParams)
    Class(ClikLikelihood) :: like
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions) Theory
    real(mcp) DataParams(:)
    integer :: i,j ,l
    real(mcp) acl(lmax,num_cls_tot)
    real(dp) clik_cl_and_pars(like%clik_n)

    call Theory%ClsFromTheoryData(acl)

    !set C_l and parameter vector to zero initially
    clik_cl_and_pars = 0.d0

    j = 1

    !TB and EB assumed to be zero
    !If your model predicts otherwise, this function will need to be updated
    do i=1,4
        do l=0,like%clik_lmax(i)
            !skip C_0 and C_1
            if (l >= 2) then
                clik_cl_and_pars(j) = acl(l,mapped_index(i))
            end if
            j = j+1
        end do
    end do

    !Appending nuisance parameters
    !Not pretty. Oh well.
    if (like%clik_nnuis > 0) then
        do i=1,like%clik_nnuis
            !silvia--------------------------------------------
            if(i .eq. 6 .and. INDEX(trim(like%name),'plik').ne. 0)then
                clik_cl_and_pars(j) = DataParams(i)-2
            else
                clik_cl_and_pars(j) = DataParams(i)
            end if
            !------------------------------------------------
            j = j+1
        end do
    end if

    !Get - ln like needed by CosmoMC
    clik_lnlike = -1.d0*clik_compute(like%clikid,clik_cl_and_pars)

    if (Feedback>1) Print*,trim(like%name)//' lnlike = ',clik_lnlike

    end function clik_lnlike

    subroutine do_checks(like)
    class (ClikLikelihood) :: like
    integer i

    !Safeguard
    if ((lmax < maxval(like%clik_lmax)+500) .and. (lmax < 4500)) then
        print*,'lmax too low: it should at least be set to',min(4500,(maxval(like%clik_lmax)+500))
        call MPIstop
    end if

    if (like%clik_nnuis/= like%nuisance_params%nnames) &
    call MpiStop('clik_nnuis has different number of nuisance parameters than .paramnames')
    if (like%clik_nnuis /= 0 .and. MPIRank==0 .and. Feedback>0) then
        Print*,'Clik will run with the following nuisance parameters:'
        do i=1,like%clik_nnuis
            Print*,trim(like%names(i))
        end do
    end if
    if (like%clik_nnuis >0) deallocate(like%names)

    end subroutine do_checks

    subroutine clik_lensing_likeinit(like, fname)
    class (ClikLensingLikelihood) :: like
    character(LEN=*), intent(in) :: fname
    integer i

    if (Feedback > 1) Print*,'Initialising clik lensing...'
    call clik_lensing_init(like%clikid,fname)
    call clik_lensing_get_lmax(like%clikid,like%lensing_lmax)

    like%clik_nnuis = clik_lensing_get_extra_parameter_names(like%clikid,like%names)
    call like%do_checks()

    like%clik_n = 2*(like%lensing_lmax+1) + like%clik_nnuis

    end subroutine clik_lensing_likeinit

    real(mcp) function clik_lensing_lnlike(like, CMB, Theory, DataParams)
    Class(ClikLensingLikelihood) :: like
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions) Theory
    real(mcp) DataParams(:)
    integer :: i,j ,l
    real(mcp) acl(lmax,num_cls_tot)
    real(dp) clik_cl_and_pars(like%clik_n)

    call Theory%ClsFromTheoryData(acl)

    !set C_l and parameter vector to zero initially
    clik_cl_and_pars = 0.d0

    j = 1
    do l=0,like%lensing_lmax
        !skip C_0 and C_1
        if (l >= 2) then
            clik_cl_and_pars(j) = acl(l,num_cls+1)/real(l*(l+1),mcp)**2*twopi
        end if
        j = j+1
    end do

    do l=0,like%lensing_lmax
        !skip C_0 and C_1
        if (l >= 2) then
            clik_cl_and_pars(j) = acl(l,mapped_index(1))
        end if
        j = j+1
    end do

    do i=1,like%clik_nnuis
        clik_cl_and_pars(j) = DataParams(i)
        j = j+1
    end do

    !Get - ln like needed by CosmoMC
    clik_lensing_lnlike = -1.d0*clik_lensing_compute(like%clikid,clik_cl_and_pars)

    if (Feedback>1) Print*,trim(like%name)//' lnlike = ',clik_lensing_lnlike

    end function clik_lensing_lnlike

    end module cliklike