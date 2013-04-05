# Microsoft Developer Studio Project File - Name="StereoMatch" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=StereoMatch - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "StereoMatch.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "StereoMatch.mak" CFG="StereoMatch - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "StereoMatch - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "StereoMatch - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""$/Users/Szeliski/Projects/Misc/StereoMatch", UIHCAAAA"
# PROP Scc_LocalPath "."
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "StereoMatch - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "StereoMatch - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "StereoMatch - Win32 Release"
# Name "StereoMatch - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\BoxFilter.cpp
# End Source File
# Begin Source File

SOURCE=.\Convert.cpp
# End Source File
# Begin Source File

SOURCE=.\Convolve.cpp
# End Source File
# Begin Source File

SOURCE=.\Histogram1D.cpp
# End Source File
# Begin Source File

SOURCE=.\Image.cpp
# End Source File
# Begin Source File

SOURCE=.\ImageIO.cpp
# End Source File
# Begin Source File

SOURCE=.\main.cpp
# End Source File
# Begin Source File

SOURCE=.\maxflow\maxflow.c
# End Source File
# Begin Source File

SOURCE=.\MinFilter.cpp
# End Source File
# Begin Source File

SOURCE=.\ParameterIO.cpp
# End Source File
# Begin Source File

SOURCE=.\RefCntMem.cpp
# End Source File
# Begin Source File

SOURCE=.\StcAggregate.cpp
# End Source File
# Begin Source File

SOURCE=.\StcDiffusion.cpp
# End Source File
# Begin Source File

SOURCE=.\StcEvaluate.cpp
# End Source File
# Begin Source File

SOURCE=.\StcGraphCut.cpp
# End Source File
# Begin Source File

SOURCE=.\StcOptDP.cpp
# End Source File
# Begin Source File

SOURCE=.\StcOptimize.cpp
# End Source File
# Begin Source File

SOURCE=.\StcOptSO.cpp
# End Source File
# Begin Source File

SOURCE=.\StcPreProcess.cpp
# End Source File
# Begin Source File

SOURCE=.\StcRawCosts.cpp
# End Source File
# Begin Source File

SOURCE=.\StcRefine.cpp
# End Source File
# Begin Source File

SOURCE=.\StcSimulAnn.cpp
# End Source File
# Begin Source File

SOURCE=.\StereoIO.cpp
# End Source File
# Begin Source File

SOURCE=.\StereoMatcher.cpp
# End Source File
# Begin Source File

SOURCE=.\StereoParameters.cpp
# End Source File
# Begin Source File

SOURCE=.\Warp1D.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\BoxFilter.h
# End Source File
# Begin Source File

SOURCE=.\Convert.h
# End Source File
# Begin Source File

SOURCE=.\Convolve.h
# End Source File
# Begin Source File

SOURCE=.\Copyright.h
# End Source File
# Begin Source File

SOURCE=.\Error.h
# End Source File
# Begin Source File

SOURCE=.\Histogram1D.h
# End Source File
# Begin Source File

SOURCE=.\Image.h
# End Source File
# Begin Source File

SOURCE=.\ImageIO.h
# End Source File
# Begin Source File

SOURCE=.\maxflow\maxflow.h
# End Source File
# Begin Source File

SOURCE=.\MinFilter.h
# End Source File
# Begin Source File

SOURCE=.\ParameterIO.h
# End Source File
# Begin Source File

SOURCE=.\RefCntMem.h
# End Source File
# Begin Source File

SOURCE=.\StereoIO.h
# End Source File
# Begin Source File

SOURCE=.\StereoMatcher.h
# End Source File
# Begin Source File

SOURCE=.\StereoParameters.h
# End Source File
# Begin Source File

SOURCE=.\Verbose.h
# End Source File
# Begin Source File

SOURCE=.\Warp1D.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=".\README-StereoMatch.txt"
# End Source File
# End Target
# End Project
