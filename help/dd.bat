@echo off

rem dd.bat is part of the NONMEM DOS help system.
rem Files ~htmp5 and ~htmp7 are created by hh.bat.
rem dd displays the detailed descriptions.
rem Usage of dd:  dd n 
rem               where n refers to the nth. document listed in ~htmp7

rem Written by JM Gries and AJ Boeckmann 4/96

rem Must be executed in the same directory into which the diskette was read

@echo off

rem # if no arguments, print the list from hh
if .%1 == . goto list

rem # if argument "n" is present, create a batch file to display nth.
gawk '{if (NR==%1)print "more < " $2}' < ~htmp5 > ~htmp.bat

rem # call the batch file
call ~htmp.bat

rem #  post-display 
choice /c:yn "Display the list again? "
if errorlevel==2 goto exit

rem # answer was Y
:list
more <~htmp7
echo.
echo Enter "dd n" to display the nth. document. 

:exit
@echo on
