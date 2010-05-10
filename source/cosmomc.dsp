# Microsoft Developer Studio Project File - Name="cosmomc" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=cosmomc - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cosmomc.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cosmomc.mak" CFG="cosmomc - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cosmomc - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "cosmomc - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "cosmomc - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "cosmomc - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /define:"NOWMAP" /define:"MATRIX_SINGLE" /define:"DECONLY" /fpp /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /out:"../cosmomc.exe" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "cosmomc - Win32 Release"
# Name "cosmomc - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Group "CAMB"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\camb\bessels.f90
DEP_F90_BESSE=\
	".\Debug\lvalues.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	".\Debug\Ranges.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\camb.f90
DEP_F90_CAMB_=\
	".\Debug\CAMBmain.mod"\
	".\Debug\GaugeInterface.mod"\
	".\Debug\InitialPower.mod"\
	".\Debug\lensing.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	".\Debug\SpherBessels.mod"\
	".\Debug\ThermoData.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\cmbmain.f90
DEP_F90_CMBMA=\
	".\Debug\GaugeInterface.mod"\
	".\Debug\InitialPower.mod"\
	".\Debug\lvalues.mod"\
	".\Debug\MassiveNu.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\NonLinear.mod"\
	".\Debug\Precision.mod"\
	".\Debug\SpherBessels.mod"\
	".\Debug\ThermoData.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\equations.f90
DEP_F90_EQUAT=\
	".\Debug\lvalues.mod"\
	".\Debug\MassiveNu.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	".\Debug\ThermoData.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\halofit.f90
DEP_F90_HALOF=\
	".\Debug\ModelParams.mod"\
	".\Debug\Transfer.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\inifile.f90
# End Source File
# Begin Source File

SOURCE=..\camb\lensing.f90
DEP_F90_LENSI=\
	".\Debug\InitialPower.mod"\
	".\Debug\lvalues.mod"\
	".\Debug\ModelData.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\modules.f90
DEP_F90_MODUL=\
	".\Debug\AMLutils.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\InitialPower.mod"\
	".\Debug\Precision.mod"\
	".\Debug\Ranges.mod"\
	".\Debug\RECFAST.MOD"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\power_tilt.f90
DEP_F90_POWER=\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\recfast.f90
DEP_F90_RECFA=\
	".\Debug\Precision.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\camb\subroutines.f90
DEP_F90_SUBRO=\
	".\Debug\AMLutils.mod"\
	
# End Source File
# End Group
# Begin Source File

SOURCE=.\calclike.f90
DEP_F90_CALCL=\
	".\Debug\CMB_Cls.mod"\
	".\Debug\cmbdata.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\lya.mod"\
	".\Debug\mpk.mod"\
	".\Debug\ParamDef.mod"\
	".\Debug\Random.mod"\
	".\Debug\settings.mod"\
	".\Debug\snovae.mod"\
	".\Debug\WeakLen.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\CMB_Cls_simple.f90
DEP_F90_CMB_C=\
	".\Debug\CAMB.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\LambdaGeneral.mod"\
	".\Debug\lensing.mod"\
	".\Debug\lya.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\settings.mod"\
	".\Debug\snovae.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\cmbdata.f90
DEP_F90_CMBDA=\
	".\Debug\CMBLikes.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\MatrixUtils.mod"\
	".\Debug\settings.mod"\
	
NODEP_F90_CMBDA=\
	".\Debug\WMAP_OPTIONS.mod"\
	".\Debug\WMAP_PASS2.mod"\
	".\Debug\WMAP_UTIL.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\cmbtypes.f90
DEP_F90_CMBTY=\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\conjgrad_wrapper.f90
DEP_F90_CONJG=\
	".\Debug\CalcLike.mod"\
	".\Debug\ParamDef.mod"\
	".\Debug\Random.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\driver.F90
DEP_F90_DRIVE=\
	".\Debug\CalcLike.mod"\
	".\Debug\cmbdata.mod"\
	".\Debug\ConjGradModule.mod"\
	".\Debug\EstCovmatModule.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\MatrixUtils.mod"\
	".\Debug\MonteCarlo.mod"\
	".\Debug\mpk.mod"\
	".\Debug\ParamDef.mod"\
	".\Debug\posthoc.mod"\
	".\Debug\settings.mod"\
	".\Debug\WeakLen.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\EstCovmat.f90
DEP_F90_ESTCO=\
	".\Debug\CalcLike.mod"\
	".\Debug\MatrixUtils.mod"\
	".\Debug\ParamDef.mod"\
	".\Debug\Random.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\lya.f90
DEP_F90_LYA_F=\
	".\Debug\cmbtypes.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Matrix_utils.F90
DEP_F90_MATRI=\
	".\Debug\AMLutils.mod"\
	
NODEP_F90_MATRI=\
	".\Debug\IFPORT.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MCMC.f90
DEP_F90_MCMC_=\
	".\Debug\CalcLike.mod"\
	".\Debug\ParamDef.mod"\
	".\Debug\propose.mod"\
	".\Debug\Random.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mpk.f90
DEP_F90_MPK_F=\
	".\Debug\cmbtypes.mod"\
	".\Debug\MatrixUtils.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\paramdef.F90
DEP_F90_PARAM=\
	".\Debug\CMB_Cls.mod"\
	".\Debug\cmbdata.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\MatrixUtils.mod"\
	".\Debug\Random.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\params_CMB.f90
DEP_F90_PARAMS=\
	".\Debug\CAMB.mod"\
	".\Debug\CMB_Cls.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\ModelParams.mod"\
	".\Debug\ParamDef.mod"\
	".\Debug\Precision.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Planck_like.f90
DEP_F90_PLANC=\
	".\Debug\AMLutils.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\MatrixUtils.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\postprocess.f90
DEP_F90_POSTP=\
	".\Debug\CalcLike.mod"\
	".\Debug\CMB_Cls.mod"\
	".\Debug\cmbdata.mod"\
	".\Debug\cmbtypes.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\propose.f90
DEP_F90_PROPO=\
	".\Debug\ParamDef.mod"\
	".\Debug\Random.mod"\
	".\Debug\settings.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\settings.f90
DEP_F90_SETTI=\
	".\Debug\AMLutils.mod"\
	".\Debug\IniFile.mod"\
	".\Debug\Random.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\supernovae.f90
DEP_F90_SUPER=\
	".\Debug\CAMB.mod"\
	".\Debug\cmbtypes.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\utils.F90
DEP_F90_UTILS=\
	{$(INCLUDE)}"cxml_dll_use.mod"\
	{$(INCLUDE)}"CXML_INCLUDE.F90"\
	{$(INCLUDE)}"cxml_static_use.mod"\
	
NODEP_F90_UTILS=\
	".\Debug\F90_UNIX.mod"\
	".\Debug\IFPORT.mod"\
	".\Debug\mpif.h"\
	
# End Source File
# Begin Source File

SOURCE=.\WeakLen.f90
DEP_F90_WEAKL=\
	".\Debug\cmbtypes.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
