! A parser for SDSS Supernovae FITRES 
! Wessel Valkenburg, LAPTH, 2009
! E-mail: wessel.valkenburg@lapp.in2p3.fr
!
!
! Usage:
!  "use sn_sdss_parser"
!  
!  declare a structure of type sdssdata (already publicly defined in this module)
! e.g.  
!  "type(sdssdata) :: thisdata"
!
! then simply call GetSDSSSN
!  "call GetSDSSSN('a', thisdata)"
! 
! where a can be any single character a,b,c,d,e or f, specifying the subset of SNe:
!   a: SDSS only
!   b: SDSS+ESSENCE+SNLS
!   c: LOWZ+SDSS
!   d: LOWZ+SDSS+ESSENCE+SNLS
!   e: LOWZ+SDSS+ESSENCE+SNLS+HST
!   f: LOWZ+ESSENCE+SNLS
!
! or
!  "call GetSDSSSN('1 3 4 50 100', thisdata)"
! with any space-seperated list of arbitrary selectrion of the following id's:
!
! 1   = SDSS
! 3   = ESSE
! 4   = SNLS
! 50  = Nearby
! 100 = HST
!
! eg. 'e' is equivalent to '1 3 4 50 100' 
!
! the structure thisdata will be allocated and returned containing
! thisdata%z(:)      - redshift
! thisdata%zerr(:)   - error in redshift
! thisdata%mu(:)     - apparent brightness
! thisdata%muerr(:)  - error on brightness
! thisdata%sigint(:) - intrinsic error
!
! all these arrays will have exactly the size that corresponds to the number of SNe
! in the chosen subset.
!
!!
!AL sept 09: changed not to use non-standard allocatable arrays in derived types
!            .dataset read using standard INI framework
module sn_sdss_parser
  use AmlUtils
  use IniFile
  implicit none
  private

  logical :: sdss_sn_ini = .false.
  integer, parameter :: maxit = 10000, sn_len=128, sn_head=32, sn_dlen=Ini_max_string_len
  integer, parameter :: dp=KIND(1.d0)
!  character(len=32), parameter :: sntag(3)=/'NVAR:','VARNAMES:','SN:'/

  type ReallocArray2
    real(dp), dimension(:,:), pointer :: Arr
  end type ReallocArray2

  type sdssvars
     character(len=sn_dlen) :: filename(2)
     character(len=sn_head) :: tagnvar, tagvarnames, tagsn, subset
     character(len=sn_head), pointer :: colnames(:)
     character(len=1) :: presets(6)
     integer :: ntags, telids(10), set
     real(dp), pointer :: sdss_sigint(:)
  end type sdssvars
  type(sdssvars) :: sv

  type sdssdata
     real(dp), pointer :: mu(:), muerr(:), z(:), sigint(:), sigz(:)
  end type sdssdata

  public GetSDSSSn, sdssdata

contains

  subroutine GetSDSSSn(thisfilename, thedata, feedback)
    character(len=*) :: thisfilename
!    character(len=sn_head) :: thissubset, thissubsetb
!    integer :: thisset, thissetb
    integer :: feedback
    type(sdssdata) :: thedata
    character(len=128) :: mystring

    call ReadSDSSdataset(thisfilename)

!    thissubsetb = thissubset
!    thissetb = thisset
    call InitializeSN(sv%set, sv%subset)
    

    call ReadSDSSSN(sv%set, sv%subset, thedata)

    ! Feedback which telescopes are uesd.
    if(feedback>0)then
       write(*,*)'Reading: '//trim(adjustl(sv%filename(sv%set)))
       write(mystring,*)size(thedata%z)
       mystring = 'Using '//trim(adjustl(mystring))//' SNe'
       mystring = trim(adjustl(mystring))//' from subset '
       if(any(sv%telids.eq.1))mystring = trim(adjustl(mystring))//' SDSS +'
       if(any(sv%telids.eq.3))mystring = trim(adjustl(mystring))//' ESSENCE +'
       if(any(sv%telids.eq.4))mystring = trim(adjustl(mystring))//' SNLS +'
       if(any(sv%telids.eq.50))mystring = trim(adjustl(mystring))//' LOWZ +'
       if(any(sv%telids.eq.100))mystring = trim(adjustl(mystring))//' HST +'
       mystring=mystring(1:len(trim(mystring))-2)//'  '
       write(*,*)trim(adjustl(mystring))//'.'
    end if

  end subroutine GetSDSSSn


  subroutine InitializeSN(thisset, thissubset)
    implicit none
    character :: mychar
    character(len=*) :: thissubset
    integer :: ntel, thisset, maxtel


    if(sdss_sn_ini)stop 'Calling InitializeSN more than once.'

    ! All in upper case
    sv%tagnvar = 'NVAR:'
    sv%tagvarnames = 'VARNAMES:'
    sv%tagsn = 'SN:'

    ! These are the names of the columns, where MU gets the extra character that is the input in
    ! GetSDSSSn('x',...). That is, the actual column becomes sv%colnames(2)//'x'. If this is 
    ! not present in your file, then input ' ' in stead of 'x'.

    sv%ntags=5

    allocate(sv%colnames(sv%ntags))

    sv%colnames(1) = 'Z'
    sv%colnames(2) = 'MU'
    sv%colnames(3) = 'MUERR'
    sv%colnames(4) = 'ZERR'
    sv%colnames(5) = 'IDTEL'
    
    ! If SALT-II, then:
    ! attach 'e' because we want the column that contains just all the data.
    ! We then select the proper telescopes for our subset ourselves, as the 
    ! a - f subsets are hardcoded in the next lines.
    if (thisset .eq.2) then 
       sv%colnames(2) = trim(adjustl(sv%colnames(2))) //'E'    
       sv%colnames(3) = trim(adjustl(sv%colnames(3))) //'E'
    end if

    ! interpret thisset:
    thissubset = adjustl(thissubset)
    sv%telids = 0
    sv%presets(:)=(/'a','b','c','d','e','f'/)

    mychar = thissubset(1:1)
    if(ichar(thissubset(1:1)).lt.58)then ! then it's a number
       ntel = CountItemsIn(trim(adjustl(thissubset)))
       read(thissubset,*)sv%telids(1:ntel)
    else if (any(thissubset(1:1).eq.sv%presets))then
       if(thissubset(1:1).eq.'a')sv%telids(1)=1
       if(thissubset(1:1).eq.'b')sv%telids(1:3)=(/1,3,4/)
       if(thissubset(1:1).eq.'c')sv%telids(1:2)=(/1,50/)
       if(thissubset(1:1).eq.'d')sv%telids(1:4)=(/1,3,4,50/)
       if(thissubset(1:1).eq.'e')sv%telids(1:5)=(/1,3,4,50,100/)
       if(thissubset(1:1).eq.'f')sv%telids(1:3)=(/3,4,50/)
    else
       stop 'No subset defined for sn_sdss!'
    end if


    ! set intrinsic uncertainty for this dataset
    maxtel = 100 ! highest telesope id
!    if(allocated(sv%sdss_sigint))deallocate(sv%sdss_sigint)
    allocate(sv%sdss_sigint(1:maxtel))
    
    ! As specified in arXiv:0908.4274v1 
    ! table 11 for SALT and paragraph 8.1 for MLCS

    if(thisset.eq.1)then ! MLCS
       sv%sdss_sigint = 0.16 ! use same intrinsic error for all telescopes
    else if (thisset .eq.2) then !SALT-II
       sv%sdss_sigint = 0.d0
       sv%sdss_sigint(1) = 0.09 ! SDSS
       sv%sdss_sigint(3) = 0.13 ! ESSENCE
       sv%sdss_sigint(4) = 0.16 ! SNLS
       sv%sdss_sigint(50) = 0.16 ! LOWZ
       sv%sdss_sigint(100) = 0.15 ! HST
    else
       stop 'Need to define intrinsic errors for this set.'
    end if
    sdss_sn_ini = .true.
    
    
  end subroutine InitializeSN


  ! Open the actual file and process it, returns thedata containing
  ! only the SNe from subset of nset.
  subroutine ReadSDSSSN(nset, thisset, thedata)
    character(len=512) :: thisline
    character(len=sn_head) :: thistag, thisfmt, thisname
    character(len = 1) :: thisset
    integer :: fnum,nset,myiostat
    logical :: snexist

    integer :: nvar, i
    Type(ReallocArray2) ::sdss_file
    character(len=sn_head), allocatable :: header(:)
    
    type(sdssdata) :: thedata

    fnum=NewIONum()

    inquire(file=sv%filename(nset),exist=snexist)
    if(.not.snexist)then
       write(*,*)'File '//trim(sv%filename(nset))//' not found. ABORTING.'
       stop 'File not found.'
    end if 
   open(fnum, file=sv%filename(nset), status='OLD', iostat=myiostat)


    ! First get nvar
    myiostat=0
    Do while (myiostat.eq.0)

       read(fnum,'(A)',iostat=myiostat)thisline

       if(myiostat.ne.0)then
          write(*,*)'Tag '//trim(sv%tagnvar)//' not found in file '//trim(sv%filename(nset))//' not found. ABORTING.'
          stop 'Tag not found.'
       end if

       thisline=trim(adjustl(thisline))
       thistag = thisline(1:len(trim(sv%tagnvar)))
       call s_cap(thistag)
       if(trim(thistag).eq.sv%tagnvar)then
          call RemoveExp(thisline,trim(sv%tagnvar))
          read(thisline,*)nvar
          exit
       end if
    end Do


    ! Get header
    allocate(header(nvar))
    myiostat=0
    Do while (myiostat.eq.0)
       read(fnum,'(A)',iostat=myiostat)thisline
       if(myiostat.ne.0)then
          write(*,*)'Tag '//trim(sv%tagvarnames)//' not found in file '//trim(sv%filename(nset))//' not found. ABORTING.'
          stop 'Tag not found.'
       end if
       thisline=trim(adjustl(thisline))
       thistag =thisline(1:len(trim(sv%tagvarnames)))
       call s_cap(thistag)
       if(trim(thistag).eq.trim(sv%tagvarnames))then
          call RemoveExp(thisline,trim(sv%tagvarnames))
          read(thisline,*)header(:)
          exit
       end if
    end Do

    ! make header data uppercase:
    do i=1,nvar
       call s_cap(header(i))
    end do

    ! We could make the code more flexible
    ! in order to deal with CID in any other column,
    ! but why bother...
    if(trim(adjustl(header(1))).ne.'CID')then
       stop 'Wrong file layout: first column must be SN:, second column CID.'
    end if
    ! anticipate the ignoring of the CID, the real array sdss_file starts 
    ! at second column in stead of first.
    header=header(2:)
    header(size(header))=''


    ! Get data
    call AllocDble(sdss_file,nvar-1)
    myiostat=0
    i=0
    Do while (myiostat .eq.0)
       read(fnum,'(A)',iostat=myiostat)thisline
       if(myiostat.ne.0)exit
       thisline=trim(adjustl(thisline))


       ! ignore lines that do not start with 'SN:'
       if(thisline(1:len(trim(sv%tagsn))).eq.sv%tagsn)then
          
          i=i+1
          ! check array size, if necessary increase
          if(i.gt.size(sdss_file%arr(1,:)))call IncreaseDble(sdss_file)
          
          call RemoveExp(thisline,trim(sv%tagsn))
          read(thisline,*)thisname,sdss_file%arr(:,i)
       end if

    end Do

    call TrimDble(sdss_file,i)

    call TranslateSDSS(thisset, header, sdss_file, thedata)


    deallocate(header)
    deallocate(sdss_file%arr)

  end subroutine ReadSDSSSN


  ! Go from full fitres file to only interesting columns
  ! from subset.
  subroutine TranslateSDSS(thissetd, header, sdss_file, thedata)
    implicit none
    character(len=1) :: thissetd, thissetabc
    character(len=sn_head) :: header(:), tags(sv%ntags)
    Type(ReallocArray2) :: sdss_file
    type(sdssdata) :: thedata
    integer :: i,j, totsn, setsn, tagn(sv%ntags)

    thissetabc=thissetd
    call s_cap(thissetabc)

    tags(1) = trim(adjustl(sv%colnames(1)))
    tags(2) = trim(adjustl(sv%colnames(2)))
    tags(3) = trim(adjustl(sv%colnames(3)))
    ! we are taking set 'thisset', that is MUA, or MUB, or MUC etc..
    tags(4) = trim(adjustl(sv%colnames(4)))
    tags(5) = trim(adjustl(sv%colnames(5)))


    ! get column number of tags: 
    tagn=0
    do i=1,size(header)
       do j=1,sv%ntags
          if(trim(adjustl(header(i))).eq.trim(adjustl(tags(j))))tagn(j)=i
       end do
    end do
    if(any(tagn.eq.0))stop 'Tags not found in translatesdss (in sn_sdss_parser.f90).'

    ! count SN in set:
    totsn =size(sdss_file%arr(1,:)) 
    setsn=0

    do i=1,totsn
       if(any(int(sdss_file%arr(tagn(5),i)+0.001).eq.sv%telids))setsn=setsn+1
    end do


  !  if(allocated(thedata%z))deallocate(thedata%z)
  !  if(allocated(thedata%mu))deallocate(thedata%mu)
  !  if(allocated(thedata%muerr))deallocate(thedata%muerr)
  !  if(allocated(thedata%sigint))deallocate(thedata%sigint)
  !  if(allocated(thedata%sigz))deallocate(thedata%sigz)
    allocate(thedata%z(setsn))
    allocate(thedata%mu(setsn))
    allocate(thedata%muerr(setsn))
    allocate(thedata%sigint(setsn))
    allocate(thedata%sigz(setsn))

    j=0
    do i=1,totsn

!       if(sdss_file(tagn(2),i).gt.-8.9999)then
       if(any(int(sdss_file%arr(tagn(5),i)+0.001).eq.sv%telids))then
          j=j+1
          thedata%z(j)=sdss_file%arr(tagn(1),i)
          thedata%mu(j)=sdss_file%arr(tagn(2),i)
          thedata%muerr(j)=sdss_file%arr(tagn(3),i)
          thedata%sigz(j)=sdss_file%arr(tagn(4),i)
          thedata%sigint(j)=sv%sdss_sigint(int(sdss_file%arr(tagn(5),i)+0.001))
       end if

    end do

    if(j.ne.setsn) stop 'There is a bug in the code, TranslateSDSS in sn_sdss_parser.f90.'

  end subroutine TranslateSDSS

  ! Allocate a double precision array with two indices
  subroutine AllocDble(A,somesize1)
    Type(ReallocArray2) :: A
    integer :: somesize1

    !if(associated(A%arr))deallocate(A%arr)
    allocate(A%arr(somesize1,sn_len))

  end subroutine AllocDble

  ! Increase the size in second index of a double precision array with two indices
  subroutine IncreaseDble(A)
    Type(ReallocArray2) :: A
    real(dp), allocatable :: arr2(:,:)
    integer :: arrsize(2)

    arrsize(1) = size(A%arr(:,1))
    arrsize(2) = size(A%arr(1,:))
    allocate(arr2(arrsize(1),arrsize(2)))
    arr2=A%arr
    deallocate(A%arr)
    allocate(A%arr(arrsize(1),arrsize(2)+sn_len))
    A%arr=0
    A%arr(:,1:arrsize(2))=arr2(:,:)

    deallocate(arr2)
  end subroutine IncreaseDble

  ! Trim in second index of a double precision array with two indices to size 'finsize'
  subroutine TrimDble(A,finsize)
    Type(ReallocArray2) :: A
    real(dp), allocatable :: arr2(:,:)
    integer, optional :: finsize
    integer :: arrsize(2)

    if(.not.present(finsize))call fatalerror('TrimSing called without final size. Faulty code.')

    arrsize(1) = size(A%arr(:,1))
    arrsize(2) = size(A%arr(1,:))
    allocate(arr2(arrsize(1),arrsize(2)))
    arr2=A%arr
    deallocate(A%arr)
    allocate(A%arr(arrsize(1),finsize))
    A%arr=0
    A%arr=arr2(:,1:finsize)

    deallocate(arr2)
  end subroutine TrimDble

  ! return a free io-number for file-io
  function NewIONum()
    integer :: NewIONum, i
    logical :: notfree

    notfree = .true.
    i=20
    
    do while (notfree)
       i=i+1
       inquire(i,opened=notfree)
       if(i.gt.maxit)call FatalError('NewIONum could not find free filenumber.')
    end do
    NewIONum = i
    
  endfunction NewIONum
  
  ! open existing file for append, or create if nonexistent.
  subroutine OpenAppend(fnum,fname)
    character(len=*) :: fname
    integer :: fnum
    logical fpresent
    
    inquire(file=trim(adjustl(fname)),exist=fpresent)
    
    if(fpresent)then
       open(fnum,file=fname,status='OLD',access='APPEND')
    else
       open(fnum,file=fname,status='NEW')
    end if
    
  end subroutine OpenAppend
  
  
  ! split string into two pieces, first half goes into word, second half
  ! into string. Split point is set by 'seperator', which itself is thrown away
  subroutine split(string,word,seperator)
    character(len=*)::string,word,seperator
    integer :: i,sep
    
    string = adjustl(string)
    sep = len(seperator)-1
    do i = 1,len(trim(string))!+2
       if(string(i:i+sep).eq.seperator)then
          word=trim(string(1:i-1))
          string=trim(string(i+sep+1:))
          exit
       end if
    enddo
  end subroutine split
  
  ! Remove an expression from string.
  subroutine RemoveExp(string,seperator)
    character(len=*)::string,seperator
    integer :: i,sep
    
    !string = adjustl(string)
    sep = len(seperator)
    do i = 1,len(trim(string))-sep+1
       if(string(i:i+sep-1).eq.seperator)then
          if(i.gt.1)then
             string=string(:i-1)//trim(string(i+sep:))
          else
             string=trim(string(i+sep:))
          end if
       end if
    enddo
  end subroutine RemoveExp
  
  ! Replace expression 'seperator' in 'string' by 'word'.
  subroutine ReplaceExp(string,seperator,word)
    character(len=*)::string,seperator,word
    integer :: i,sep,j
    
    !if(seperator.eq.'%fnumber')then
    !debug=.false.
    !else
    !debug=.false.
    !end if
    !string = adjustl(string)
    sep = len(seperator)
    i=0
    j=0
    do
       i=i+1
       j=j+1
       !   if(debug)write(*,*)string(i:i+sep-1),',',seperator
       if(string(i:i+sep-1).eq.seperator)then
          !   if(debug)then
          !   write(*,*)'Whoopie!'
          !write(*,*)'----------------------------------------------------------'
          !    write(*,*)string(:i-1)
          !write(*,*)'----------------------------------------------------------'
          !    write(*,*)adjustl(word)
          !write(*,*)'----------------------------------------------------------'
          !    write(*,*)string(i+sep:)
          !write(*,*)'----------------------------------------------------------'
          !stop
          !    end if
          if(i.gt.1)then
             string=string(:i-1)//adjustl(word)//string(i+sep:)
          else
             string=adjustl(word)//string(i+sep:)
          end if
          if (sep.ne.len(word)) i=i+len(word)-sep
       end if
       if (i .ge. len(trim(string))-sep+1) exit
       if (j.gt. 1000*len(string))call fatalerror('ReplaceExp is hanging.')
    enddo
  end subroutine ReplaceExp
  
  ! Grab expression surrounded by initag and endtag in string. 
  ! Result is word.
  ! Stops is tags are not found.
  subroutine GetFromTags(string,word,initag,endtag)
    character(len=*), intent(in) :: string
    character(len=len(string)) :: mystring1, mystring2
    character(len=*), intent(out) :: word
    character(len=*), intent(in) :: initag
    character(len=*), intent(in) :: endtag
  
    mystring1 = string

    mystring2=''
    call split(mystring1,mystring2,initag)
    if(len(trim(mystring2)).eq.0)then
       write(*,*)'Tag '//trim(adjustl(initag))//' not found.'
       stop 'Tag not found (GetFromTags in sn_sdss_parser).'
    end if

    mystring2=''
    call split(mystring1,mystring2,endtag)
    if(len(trim(mystring2)).eq.0)then
       write(*,*)'Tag '//trim(adjustl(initag))//' not found.'
       stop 'Tag not found (GetFromTags in sn_sdss_parser).'
    end if

    word = trim(adjustl(mystring2))
    
  end subroutine GetFromTags


  ! Count the number of space-separated items in mystr
  function CountItemsIn(mystr)
    implicit none
    integer :: CountItemsIn, i, j
    character(len=*) :: mystr

    j=0
    do i=1,len(mystr)
       if(mystr(i:i).eq.' ')then
          if(i.eq.1)then
             j=j+1
          else if((i.gt.1).and.(mystr(i-1:i-1).ne.' '))then
             j=j+1
          end if
       end if
    end do
    CountItemsIn = j + 1

  end function CountItemsIn

! subroutine ch_cap taken from chrpak at
! http://orion.math.iastate.edu/burkardt/f_src/chrpak/chrpak.html
! at 08/09/2009
  subroutine ch_cap ( c )
    !
    !*******************************************************************************
    !
    !! CH_CAP capitalizes a single character.
    !
    !
    !  Modified:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character C, the character to capitalize.
    !
    implicit none
    !
    character c
    integer itemp
    !
    itemp = ichar ( c )
    
    if ( 97 <= itemp .and. itemp <= 122 ) then
       c = char ( itemp - 32 )
    end if
    
    return
  end subroutine ch_cap
  
  
! subroutine ch_cap taken from chrpak at
! http://orion.math.iastate.edu/burkardt/f_src/chrpak/chrpak.html
! at 08/09/2009
  subroutine s_cap ( s )
    !
    !*******************************************************************************
    !
    !! S_CAP replaces any lowercase letters by uppercase ones in a string.
    !
    !
    !  Modified:
    !
    !    28 June 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !
    implicit none
    !
    character c
    integer i
    integer nchar
    character ( len = * ) s
    !
    nchar = len_trim ( s )

    do i = 1, nchar

       c = s(i:i)
       call ch_cap ( c )
!       if(len(s).eq.1)then
!          s=''
!          s=c
!       else
       s(i:i) = c
!       end if

    end do
    return
  end subroutine s_cap
  
  subroutine FatalError(msg)
    character(len=*) :: msg
  
       write(*,*) msg
       call MpiStop('sn_sdsss_parser: ABORTING DUE TO FATAL ERROR.')
    
  end subroutine FatalError

  
  ! Open and parse the .dataset-file which contains
  ! the parameters needed in order to 
  ! parse the FITRES files
  ! (paths, set and subset)
  subroutine ReadSDSSdataset(thisfilename)
    use Settings
    character(len=*) :: thisfilename 
    character(len=sn_head) :: setstring
    integer :: fnum
    type(TIniFile) :: Ini
    logical bad

    
    fnum = NewIONum()
    call Ini_Open_File(Ini, thisfilename, fnum, bad, .false.)
    
    if(bad)then
       write(*,*)'File '//trim(thisfilename)//' not found. ABORTING.'
       call MpiStop('File not found.')
    end if
   
    sv%filename(1) = ReadIniFileName(Ini,'MLCSfile')
    sv%filename(2) = ReadIniFileName(Ini,'SALTfile')
    setstring = Ini_Read_String_File(Ini,'SDSS_set')
    sv%subset =  Ini_Read_String_File(Ini,'SDSS_subset')
  
    read(setstring,*)sv%set
    call Ini_Close_File(Ini)
    
  end subroutine ReadSDSSdataset

end module sn_sdss_parser
