C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101   10:57 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::FileOC
      Subroutine FileOC(IFlag)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE     not needed					sld01
cc    FileOC is called from routines:  LHS,RDPAR2                       sld01
cc    FileOC does not call any other external routines                  sld01
c
c     -- This subroutine opens or closes all of the files LHS needs to run
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c     -- If IFlag > 0, open the files, otherwise, close them.
c
      If (IFlag <= 0) Then
c
         Close (1)
         Close (2)
         Close (3)
         CLOSE (5)
         Close (7)
         Close (8)
         Close (9)
c        -- Since this is a normal program termination, close the LHS
c        -- Error file with Status=Delete so that other processors can
c        -- know that the program terminated successfully.
c 11-96 change -----
c Note:  Closing of the error file has been moved to the main routine (i.e.
c        either LHSDLL.FOR or LHSDRV.FOR to accomodate conversion to DLLs)
c         Close (99, Status='DELETE')
c
      Else
c
c 11-96 change ----
c Note:  Opening of error file 99 has been moved to main routine (i.e.
c         either LHSDLL.FOR or LHSDRV.FOR to accomodate conversion to DLLs)
c        -- Processing to open the LHS Error File
c
c        -- We want to be sure that the error file is created anew for this
c        -- run in order to assure that it has the proper date and time stamp
c        -- information.  Thus, we open it with Status=Unknown and immediately
c        -- close it with Status=Delete.  Then we re-open the file and write
c        -- an error message to it.  We then close it with Status=Keep in
c        -- order to assure that the program buffers are flushed and the file
c        -- will really exist on the disk.  Finally, we re-open the file and
c        -- leave it open until the program executes a normal termination.
c        -- If the program crashes, the error file will remain in existance.
c        -- However, on normal termination, the error file will be deleted.
c        -- Thus, an external program can check for normal termination by
c        -- checking for the existance of this file.
c         Open (99, File='LHS.ERR', Status='UNKNOWN', Form='FORMATTED')
c         Write (99,*) 'One line into the file just to be sure...'
c         Close  (99, Status='DELETE')
c         Open (99, File='LHS.ERR', Status='NEW', Form='FORMATTED')
c         Write (99,*) 'An error occurred during LHS processing.'
c         Write (99,*) 'Consult the message file for additional ',
c     1      'information.'
c         Close (99, Status='KEEP')
c         Open (99, File='LHS.ERR', Status='OLD', Form='FORMATTED')
c
c        -- Open the LHS Scratch Files
c
c        -- Note that these are only the scratch files.  The input and
c        -- output files are opened in RDPAR and RDPAR2.
c
C$$$  Modified by slbrow, 12-23-2004
C$$$     Open (2,  Form='UNFORMATTED', Status='SCRATCH') 
C$$$     Open (3,  Form='UNFORMATTED', Status='SCRATCH') 
C$$$     Open (7,  Form='UNFORMATTED', Status='SCRATCH') 
C$$$     Open (8,  Form='UNFORMATTED', Status='SCRATCH') 
C$$$     Open (9,  Form='UNFORMATTED', Status='SCRATCH') 
C    Modified by dmdunla, 09-21-2006
C    Added file names to avoid conflicts with other Fortran temp files
         Open (2,File='LHS_2.out',Form='UNFORMATTED',Status='UNKNOWN') 
         Open (3,File='LHS_3.out',Form='UNFORMATTED',Status='UNKNOWN') 
         Open (7,File='LHS_7.out',Form='UNFORMATTED',Status='UNKNOWN') 
         Open (8,File='LHS_8.out',Form='UNFORMATTED',Status='UNKNOWN') 
         Open (9,File='LHS_9.out',Form='UNFORMATTED',Status='UNKNOWN') 
c        Modified by SFW on 3/21/2002 for Solaris 
c        Should be ifdef'd around in the future     
c        Open (2, File='S2', Form='UNFORMATTED',Status='SCRATCH')
c        Open (3, File='S3', Form='UNFORMATTED',Status='SCRATCH')
c        Open (7, File='S7', Form='FORMATTED'), Status='SCRATCH')
c        Open (8, File='S8', Form='UNFORMATTED',Status='SCRATCH')
c        Open (9, File='S9', Form='UNFORMATTED',Status='SCRATCH')
c
         Rewind 2
         Rewind 3
         Rewind 7
         Rewind 8
         Rewind 9
c
      End If
c
      Return
c
      End
