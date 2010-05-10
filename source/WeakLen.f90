
MODULE WeakLen
 use cmbtypes
  !Dummy module

  logical :: Use_WeakLen = .false.
    
CONTAINS
  
  FUNCTION WeakLenLnLike(CMB, Theory)
    TYPE (CMBParams) CMB
    TYPE (CosmoTheory) Theory
    REAL :: WeakLenLnLike

       stop 'Weak lensing not implemented'

   WeakLenLnLike = 0
    
  END FUNCTION WeakLenLnLike

end module

