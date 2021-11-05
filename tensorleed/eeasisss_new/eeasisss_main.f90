!                                 Supplemental Material
!
!           Low-Energy Electron Diffraction With Energy Invariant Carrier Wave
!                Wavenumber Modulated by Exchange-Correlation Interaction
!
!                John Rundgren,1, Wolfgang Moritz,2, and Bo E. Sernelius,3
!
!     1,Department of Physics, KTH Royal Institute of Technology, 10691 Stockholm, Sweden
!      2,Department of Earth and Environmental Sciences, Ludwig-Maximilians-University,
!                      Theresienstrasse 41, 80333 Munich, Germany
!3,Department of Physics, Chemistry and Biology, Linköping University, 58183 Linköping, Sweden
!
  !=====================================================================
  !This program is free software under the terms of the GNU General Public
  !License as published by the Free Software Foundation.
  !Author: John O. Rundgren, jru@KTH.se ,
  !        KTH Royal Institute of Technology, Stockholm, Sweden.
  !Version: 28 March 2021.
  !-----------------------------------------------------------------------
  ! Adapted for ViPErLEED by Alexander M. Imre, October 2021
  !-----------------------------------------------------------------------
  program eeasisss_main
  !This program requires compiler adhering to Fortran 2003 standard.
  implicit none
  integer :: argc, iarg
  character(len=127) :: argv, input_file, log_file, output_dir, atom_dir
  character(len=127) :: prefix
  interface
    subroutine eeasisss(inpfile, logfile, outdir, atdir)
    character(len=127) :: inpfile, logfile, outdir, atdir
    end subroutine
  end interface

  !Set default arguments:
  input_file = 'inputX'
  log_file   = 'logX'
  output_dir = './out/'
  atom_dir   = './atlib/'  !windows 10 with user name olof.
  !Accessing file whose name begins with '~/' :
  if(atom_dir(1:2)=='~/')then
    call getenv('HOME',prefix)
    atom_dir=trim(prefix)//atom_dir(2:len_trim(atom_dir))
  endif
 
  !Check if any arguments are found.
  argc = command_argument_count()
  iarg=1
  !loop across options.
  do while(iarg <= argc)
    call get_command_argument(iarg, argv) 
    selectcase(adjustl(argv))
      case("--help", "-h")
write(*,*)"Object files are ~/bin/eeasisss and ~/bin/eeas, "
write(*,*)"working directory is ./ ."
write(*,*)"===================================================================="
write(*,*)"--help or -h shows working space organization:"
write(*,*)"FILES"
write(*,*)"--input, -i <name of input file>, default is inputX ."
write(*,*)"--log, -l <name of log file>, default is logX ."
write(*,*)"DIRECTORIES"
write(*,*)"--output_dir, -o <path to log file and phase shifts>, "
write(*,*)"  default is ./result/ ."
write(*,*)"--atom_dir, -a <path to charge density files>, default is ~/atlib/ ."
write(*,*)"--------------------------------------------------------------------"
write(*,*)"Please contact John Rundgren <jru@KTH.se> for queries and comments."
        stop
      case("--input", "-i")
        if(iarg+1 .le. argc)then
          call get_command_argument(iarg+1, argv)
          input_file = adjustl(argv)
          iarg = iarg + 1
        endif
      case("--log", "-l")
        if(iarg+1 .le. argc)then
          call get_command_argument(iarg+1, argv)
          log_file = adjustl(argv)
          iarg = iarg + 1
        endif
      case("--output_dir", "-o")
        if(iarg+1 .le. argc)then
          call get_command_argument(iarg+1, argv)
          output_dir = adjustl(argv)
          iarg = iarg + 1
        endif
      case("--atom_dir", "-a")
        if(iarg+1 .le. argc)then
          call get_command_argument(iarg+1, argv)
          atom_dir = adjustl(argv)
          iarg = iarg + 1
        endif
    endselect
    iarg = iarg + 1
  end do
  call eeasisss(input_file,log_file,output_dir,atom_dir)
  stop
  end program eeasisss_main
!------------------------------------------------------------------------------
