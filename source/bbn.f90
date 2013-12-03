    !Gets predicted BBN mass fraction YHe by interpolating in pre-computed grid

    !Module originally from Jan Hamann, 4/2010
    !Rewrittten with updated table from AlterBBN by AL Dec 2013
    module bbn
    use settings
    use Interpolation
    implicit none

    type, extends(InterpGrid2D) :: TBBNPredictions
    contains
    procedure Error => BBNPredictions_error
    procedure FirstUse =>  BBN_Init
    end type TBBNPredictions

    type(TBBNPredictions) BBNPredictions

    contains

    subroutine BBN_Init(W)
    class(TBBNPredictions):: W

    if (feedback .ge. 1) print*,'Initialising BBN Helium data...'

    call W%InitFromFile(trim(DataDir)//'BBN_full_alterBBN_880.1.dat', xcol=1,ycol=3,zcol=5)
    
    if (feedback .ge. 1) print*,'Done. Interpolation table is ', W%nx,' by ',W%ny

    end subroutine BBN_Init

    subroutine BBNPredictions_error(W,S)
    class(TBBNPredictions):: W
    character(LEN=*), intent(in) :: S

    call MpiStop('BBN Error: '//trim(S))

    end subroutine BBNPredictions_error



    end module bbn

