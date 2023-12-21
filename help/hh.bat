@echo off

rem hh.bat is part of the NONMEM DOS help system.
rem Usage:  hh topic
rem Files ~htmp5 and ~htmp7 are created by hh.bat for use by dd.bat,
rem which displays the detailed description files.

rem Written by JM Gries and AJ Boeckmann 4/96

rem Must be executed in the same directory into which the diskette was read

rem # if no arguments, print listhelp document

if '%1.x' == '.x' goto usage
goto cont1
:usage
type listhelp | more
goto exit

:cont1
rem # if help help, give general advice
if '%1' == 'help' goto helpp
if '%1' == 'nmhelp' goto helpp
if '%1' == 'HELP' goto helpp
if '%1' == 'NMHELP' goto helpp
goto cont2

:helpp
type helphelp | more
goto exit

:cont2

rem # A leading blank line seems to be helpful
echo.  > ~htmp1

rem # test argument
if '%1' == '-a' goto all
if '%1' == '-w' goto word
     
rem # this is default grep
set name=%1
grep -i " %name%" < index >> ~htmp1 
goto drop

rem # 
:all
set name=%2
grep -i  "%name%" < index >> ~htmp1
goto drop

:word
set name=%2
grep -iw  "%name%" < index >> ~htmp1

rem # drop the keywords
:drop
gawk '{FS="~";print $1 ,$2}' < ~htmp1 | sort > ~htmp2

rem # drop duplicate lines ('uniq')
gawk '{if(last!=$0)print $0;last=$0}' < ~htmp2 > ~htmp3

rem # drop the leading blank
gawk '{if(NR!=1)print $0}' < ~htmp3 > ~htmp4

rem # make a file with numbers and filenames
gawk '{print NR,$1}' < ~htmp4 > ~htmp5

rem # make a file with numbers and titles, for display
gawk '{$1="";print NR,$0}' < ~htmp4 > ~htmp6

rem # make a two-column format
gawk -f hh.awk < ~htmp6 > ~htmp7

rem # clean up

erase ~htmp1 
erase ~htmp2 
erase ~htmp3 
erase ~htmp4 

echo Help for %name% is available in the following documents:
more < ~htmp7

echo.
echo Enter "dd n" to display the nth. document, e.g., dd 1
echo to display the first document.
:end
:exit
@echo on
