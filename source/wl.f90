! Module for galaxy weak lensing

module wl
  use settings
  use CosmologyTypes
  use CosmoTheory
  use Calculator_Cosmology
  use Likelihood_Cosmology
  implicit none
  private

  type, extends(TCosmoCalcLikelihood) :: WLLikelihood
     real(mcp), allocatable, dimension(:,:) :: wl_invcov,wl_cov(:)
   contains
     procedure :: LogLike => WL_LnLike
     procedure :: ReadIni => WL_ReadIni
     procedure, private :: WL_CFHTLENS_loglike
  end type WLLikelihood

  logical :: use_wl_lss  = .false.

  public WLLikelihood, WLLikelihood_Add, use_wl_lss
  contains

  subroutine WLLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WLLikelihood), pointer :: this
    integer numwlsets, i

    if (Ini%Read_Logical('use_WL',.false.)) then
       use_wl_lss = .true.
       numwlsets = Ini%Read_Int('wl_numdatasets',0)
       do i= 1, numwlsets
          allocate(this)
          this%needs_nonlinear_pk = .true.
          this%kmax=100.0
          call this%ReadDatasetFile(Ini%ReadFileName(numcat('wl_dataset',i)))
          this%LikelihoodType = 'WL'
          this%needs_powerspectra = .true.
          this%num_z = Ini%Read_Int('nz_wl',100)
          this%max_z = Ini%Read_Double('max_z',10.0d0)
       end do
       call LikeList%Add(this)
    end if

  end subroutine WLLikelihood_Add

  subroutine WL_ReadIni(this, Ini)
    class(WLLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: measurements_file, cov_file
    Type(TTextFile) :: F

    if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

  end subroutine WL_ReadIni

  function WL_LnLike(this, CMB, Theory, DataParams)
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WL_LnLike

    WL_LnLike=0

    if (this%name=='CFHTLENS_1bin' .or. this%name=='CFHTLENS_6bin') then
       WL_LnLike = this%WL_CFHTLENS_loglike(CMB,Theory)
    end if
    
  end function WL_LnLike

  ! CFHTLENS likelihood
  function WL_CFHTLENS_loglike(this,CMB,Theory)
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) WL_CFHTLENS_loglike
    WL_CFHTLENS_loglike = 0

    write(*,*) WL_CFHTLENS_loglike

  end function WL_CFHTLENS_loglike

end module wl
  
  
