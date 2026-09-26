! MarfaSimple.f90 - Simplified single-file version of MARFA
! 
! Changes from original MARFA:
! - Single file instead of modular structure
! - Reads HITRAN files directly into memory (no binary intermediate files)
! - Takes all parameters from command line (no atmospheric profiles)
! - Outputs to stdout instead of binary files
! - Simplified grid calculation (keeping cascade structure but simplified implementation)
! - Single atmospheric layer calculation
! - Function-based approach to minimize side effects
!
! Usage: marfa_simple <parfile> <wv_min> <wv_max> <cutoff> <temp> <p_self> <p_foreign> <TIPS_ref> <TIPS_target>
! Example: marfa_simple CO2.par 2000 2100 125 300 0.2 1.4 286.1 291.0 > abscoef.out

module constants

    implicit none

    ! Double precision kind parameter
    integer, parameter :: DP = selected_real_kind(15, 307)
    
    ! Physical constants (CGS units consistent with original MARFA)
    real(kind=DP), parameter :: PI = 3.1415926
    real(kind=DP), parameter :: sqln2 = sqrt(log(2.))
    real(kind=DP), parameter :: PLANCK = 6.626070e-27 ! [erg*s]
    real(kind=DP), parameter :: SPL = 2.99792458e10 ! [cm/s] -- speed of light
    real(kind=DP), parameter :: BOLsgs = 1.3806503e-16 ! [erg/K] -- Boltzmann constant in CGS
    real(kind=DP), parameter :: C2 = 1.438777 ! [cm * K] -- second radiation constant
    real(kind=DP), parameter :: refTemperature = 296. ! [K] -- reference temperature for HITRAN data
    real(kind=DP), parameter :: AVOGADRO = 6.02214076e23 ! [1/mol]
    real(kind=DP), parameter :: dopplerCONST = sqrt(2*AVOGADRO*BOLsgs*log(2.)) / SPL

end module constants


module parameters

use constants

implicit none

    ! Settings
    logical, parameter :: DEBUG = .TRUE.
    logical, parameter :: OUTPUT_SPECTRA = .TRUE.
    integer, parameter :: SPECTRA_FCHAN = 1000
    integer, parameter :: DEBUG_FCHAN = 1010
    character(len=256), parameter :: DEBUG_FNAME = "marfasimple_timings.out"
    
    ! Grid parameters (same as original MARFA)
    integer, parameter :: NT0 = 10
    integer, parameter :: NT1 = NT0 * 2
    integer, parameter :: NT2 = NT1 * 2
    integer, parameter :: NT3 = NT2 * 2
    integer, parameter :: NT4 = NT3 * 2
    integer, parameter :: NT5 = NT4 * 2
    integer, parameter :: NT6 = NT5 * 2
    integer, parameter :: NT7 = NT6 * 2
    integer, parameter :: NT8 = NT7 * 2
    integer, parameter :: NT9 = NT8 * 2
    integer, parameter :: NT = NT9 * 4 + 1
    
    real(kind=DP), parameter :: deltaWV = 10.0 ! [cm-1] resolution (subinterval width)
    real(kind=DP), parameter :: STEP = 1.0
    
    ! Grid spacing parameters
    real(kind=DP), parameter :: H0 = STEP
    real(kind=DP), parameter :: H1 = H0 / 2.0
    real(kind=DP), parameter :: H2 = H1 / 2.0
    real(kind=DP), parameter :: H3 = H2 / 2.0
    real(kind=DP), parameter :: H4 = H3 / 2.0
    real(kind=DP), parameter :: H5 = H4 / 2.0
    real(kind=DP), parameter :: H6 = H5 / 2.0
    real(kind=DP), parameter :: H7 = H6 / 2.0
    real(kind=DP), parameter :: H8 = H7 / 2.0
    real(kind=DP), parameter :: H9 = H8 / 2.0
    real(kind=DP), parameter :: H = H9 / 4.0

end module parameters


program MarfaSimple
    
    use constants
    use parameters
    
    implicit none
       
    ! Spectral line data structure
    type :: SpectralLine
        real(kind=DP) :: wv          ! wavenumber [cm-1]
        real(kind=DP) :: intensity             ! line intensity at 296K [cm-1/(molec*cm-2)]
        real(kind=DP) :: gammaForeign         ! air-broadened HWHM [cm-1/atm]
        real(kind=DP) :: gammaSelf            ! self-broadened HWHM [cm-1/atm]
        real(kind=DP) :: lowerState           ! lower state energy [cm-1]
        real(kind=DP) :: tempCoeff            ! temperature dependence coefficient
        real(kind=DP) :: deltaForeign         ! pressure shift [cm-1/atm]
        integer :: molIso            ! molecule/isotope code
    end type SpectralLine

    abstract interface
        pure function chifactor(X, moleculeIntCode)
            import :: DP
            implicit none
            real :: chiFactor
            real(kind=DP), intent(in) :: X
            integer, intent(in) :: moleculeIntCode
        end function chifactor
    end interface

    abstract interface
        real function shape(X, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            import :: DP
            implicit none
            real(kind=DP), intent(in) :: X, dopHWHM, lorHWHM, lineIntensity
            procedure(chifactor), pointer :: chiFactorFuncPtr
            integer, intent(in) :: moleculeIntCode
        end function shape
    end interface
    
    ! Command line arguments
    character(len=256) :: parFile, outFile
    real(kind=DP) :: startWV, endWV, cutOff
    real(kind=DP) :: temperature, pSelf, pForeign
    real(kind=DP) :: TIPS_ref, TIPS_target
    integer :: moleculeIntCode
        
    ! Line data array
    type(SpectralLine), allocatable :: lines(:)
    integer :: nLines

    ! pointer to the specific chi factor function which will be used in grid calculations 
    procedure(chifactor), pointer :: chiFactorFuncPtr

    ! pointer to the specific line shape function which will be used in grid calculations 
    procedure(shape), pointer :: shapeFuncPtr
    
    ! Number of subintervas
    integer :: nSubintervals, totalGridPoints
    
    ! Grid arrays for absorption calculation
    real, allocatable :: RK(:)
    real, allocatable :: RK0(:), RK0L(:), RK0P(:)
    real, allocatable :: RK1(:), RK1L(:), RK1P(:)
    real, allocatable :: RK2(:), RK2L(:), RK2P(:)
    real, allocatable :: RK3(:), RK3L(:), RK3P(:)
    real, allocatable :: RK4(:), RK4L(:), RK4P(:)
    real, allocatable :: RK5(:), RK5L(:), RK5P(:)
    real, allocatable :: RK6(:), RK6L(:), RK6P(:)
    real, allocatable :: RK7(:), RK7L(:), RK7P(:)
    real, allocatable :: RK8(:), RK8L(:), RK8P(:)
    real, allocatable :: RK9(:), RK9L(:), RK9P(:)
    
    ! Other variables
    real(kind=DP) :: EPS
    integer :: i

    ! Timing variables
    real :: start_time, end_time, elapsed_time
    integer :: count_start, count_end, count_rate, count_max

    ! Set default line shape and chi factor functions
    shapeFuncPtr => voigt
    chiFactorFuncPtr => noneChi
    
    ! Calculate number of subintervals and allocate full grid
    nSubintervals = ceiling((endWV - startWV) / deltaWV)

    ! Parse command line arguments
    call parseCommandLine(parFile,nLines,moleculeIntCode,startWV,endWV, &
        cutOff,temperature,pSelf,pForeign,TIPS_ref,TIPS_target,outFile)
    
    ! Read HITRAN file into memory with timing
    call system_clock(count_start, count_rate, count_max)
    call cpu_time(start_time)
    
    call readHITRANFile(parFile, lines, nLines, startWV, endWV, cutOff)
    
    call system_clock(count_end, count_rate, count_max)
    call cpu_time(end_time)
    
    ! Open output file
    open(SPECTRA_FCHAN, file=outFile, status='unknown')
    
    ! Open debug file
    if (DEBUG) then
        open(DEBUG_FCHAN, file=DEBUG_FNAME, status='unknown')
        write(DEBUG_FCHAN,"(a)") "isubint;ctreset;wtreset;ctfind;wtfind;" // &
            "ctprocsub;wtprocsub;ctcascint;wtcascint;ctoutsub;wtoutsub;" // &
            "startDeltaWV;endDeltaWV;deltaWV;startLineIdx;lineIdx"
    end if
    
    elapsed_time = end_time - start_time
    print *, '================================'
    print *, '=== TIMING: readHITRANFile ==='
    print *, 'CPU time: ', elapsed_time, ' seconds'
    if (count_rate > 0) then
        elapsed_time = real(count_end - count_start) / real(count_rate)
        print *, 'Wall clock time: ', elapsed_time, ' seconds'
    end if
    print *, 'Lines read: ', nLines
    print *, '================================'
    
    ! Calculate number of subintervals and allocate full grid
    nSubintervals = ceiling((endWV - startWV) / deltaWV)

    ! Allocate grid arrays
    allocate(RK(NT))
    allocate(RK0(NT0), RK0L(NT0), RK0P(NT0))
    allocate(RK1(NT1), RK1L(NT1), RK1P(NT1))
    allocate(RK2(NT2), RK2L(NT2), RK2P(NT2))
    allocate(RK3(NT3), RK3L(NT3), RK3P(NT3))
    allocate(RK4(NT4), RK4L(NT4), RK4P(NT4))
    allocate(RK5(NT5), RK5L(NT5), RK5P(NT5))
    allocate(RK6(NT6), RK6L(NT6), RK6P(NT6))
    allocate(RK7(NT7), RK7L(NT7), RK7P(NT7))
    allocate(RK8(NT8), RK8L(NT8), RK8P(NT8))
    allocate(RK9(NT9), RK9L(NT9), RK9P(NT9))
    
    ! Calculate absorption cross-sections with timing
    call system_clock(count_start, count_rate, count_max)
    call cpu_time(start_time)
    
    call calculateAbsorption(lines, nLines, shapeFuncPtr, chiFactorFuncPtr, moleculeIntCode, &
        startWV, endWV, deltaWV, cutOff, temperature, pSelf, pForeign, TIPS_ref, TIPS_target, EPS, &
        NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
        H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, nSubintervals)
    
    call system_clock(count_end, count_rate, count_max)
    call cpu_time(end_time)
    
    elapsed_time = end_time - start_time
    print *, '=== TIMING: calculateAbsorption ==='
    print *, 'CPU time: ', elapsed_time, ' seconds'
    if (count_rate > 0) then
        elapsed_time = real(count_end - count_start) / real(count_rate)
        print *, 'Wall clock time: ', elapsed_time, ' seconds'
    end if
    print *, '===================================='
    
    ! Clean up
    deallocate(lines)
        
    ! Deallocate grid arrays
    deallocate(RK)
    deallocate(RK0, RK0L, RK0P)
    deallocate(RK1, RK1L, RK1P)
    deallocate(RK2, RK2L, RK2P)
    deallocate(RK3, RK3L, RK3P)
    deallocate(RK4, RK4L, RK4P)
    deallocate(RK5, RK5L, RK5P)
    deallocate(RK6, RK6L, RK6P)
    deallocate(RK7, RK7L, RK7P)
    deallocate(RK8, RK8L, RK8P)
    deallocate(RK9, RK9L, RK9P)
    
contains

    subroutine parseCommandLine(parFile,nLines,moleculeIntCode,startWV,endWV, &
        cutOff,temperature,pSelf,pForeign,TIPS_ref,TIPS_target,outFile)
        
        use constants
        
        ! Parse command line arguments
        ! Usage: marfa_simple CO2.par 2000 2100 125 300 0.2 1.4 286.1 291.0
        implicit none
        
        character(len=256),intent(out) :: parFile, outFile
        real(kind=DP),intent(out) :: startWV, endWV, cutOff
        real(kind=DP),intent(out) :: temperature, pSelf, pForeign
        real(kind=DP),intent(out) :: TIPS_ref, TIPS_target
        integer,intent(out) :: moleculeIntCode, nLines
        
        character(len=256) :: arg
        integer :: argc
        
        argc = command_argument_count()
        if (argc /= 12) then
            print *, 'Usage: marfa_simple <parfile> <maxlines> <molec_id> <wv_min> <wv_max> <cutoff> <temp> ' // &
             '<p_self> <p_foreign> <TIPS_ref> <TIPS_target> <outfile>'
            stop 1
        end if
                
        call get_command_argument(1, parFile)
        
        call get_command_argument(2, arg)
        read(arg, *) nLines
        
        call get_command_argument(3, arg)
        read(arg, *) moleculeIntCode

        call get_command_argument(4, arg)
        read(arg, *) startWV
        
        call get_command_argument(5, arg)
        read(arg, *) endWV
        
        call get_command_argument(6, arg)
        read(arg, *) cutOff
        
        call get_command_argument(7, arg)
        read(arg, *) temperature
        
        call get_command_argument(8, arg)
        read(arg, *) pSelf
        
        call get_command_argument(9, arg)
        read(arg, *) pForeign
        
        call get_command_argument(10, arg)
        read(arg, *) TIPS_ref
        
        call get_command_argument(11, arg)
        read(arg, *) TIPS_target
        
        call get_command_argument(12, outFile)
        
    end subroutine parseCommandLine
    
    
    subroutine readHITRANFile(filename, lineData, lineCount, startWV, endWV, cutOff)
        ! Read HITRAN-formatted file directly into memory with binary caching
        ! Only reads lines within the extended spectral range [startWV-cutOff, endWV+cutOff]
        implicit none
        
        character(len=*), intent(in) :: filename
        type(SpectralLine), allocatable, intent(out) :: lineData(:)
        integer, intent(inout) :: lineCount
        real(kind=DP), intent(in) :: startWV, endWV, cutOff

        ! HITRAN format parameters (based on processParFile.f90)
        character(len=81) :: parFormat = '(I2, A1, F12.6, 2E10.3, 2F5.4, F10.4, F4.2, F8.6, 4A15, 6I1, 6I2, A1, F7.1, F7.2)'
        
        ! Variables for reading HITRAN data
        integer :: MO, IS, ios, unit_num
        character(len=1) :: ISO
        real(kind=DP) :: lineWV, refLineIntensityDP, lineLowerStateDP
        real :: A, gammaForeign, gammaSelf, foreignTempCoeff, deltaForeign
        character(len=15) :: V1, V2, Q1, Q2
        integer :: IERR(6), IREF(6)
        character(len=1) :: FLAG
        real(kind=DP) :: Gu, Gl
        
        ! Temporary storage
        type(SpectralLine), allocatable :: tempLines(:)
        integer :: currentLine
        real(kind=DP) :: extStartWV, extEndWV
        
        ! Binary cache variables
        character(len=256) :: binFilename
        logical :: binFileExists
        integer :: binUnit, nLinesInBin
        real(kind=DP) :: cachedStartWV, cachedEndWV, cachedCutOff
        
        ! Define extended range
        extStartWV = startWV - cutOff
        extEndWV = endWV + cutOff
        
        ! Create binary filename
        binFilename = trim(filename) // '.bin'
        
        ! Check if binary cache exists
        inquire(file=binFilename, exist=binFileExists)
        
        if (binFileExists) then
            ! Try to read from binary cache
            binUnit = 20
            open(binUnit, file=binFilename, form='unformatted', status='old', iostat=ios)
            if (ios == 0) then
                ! Read cache header
                read(binUnit, iostat=ios) nLinesInBin, cachedStartWV, cachedEndWV, cachedCutOff
                
                ! Check if cached data matches current request
                if (ios == 0 .and. nLinesInBin > 0 .and. &
                    abs(cachedStartWV - extStartWV) < 1e-6 .and. &
                    abs(cachedEndWV - extEndWV) < 1e-6 .and. &
                    abs(cachedCutOff - cutOff) < 1e-6) then
                    
                    ! Cache matches - read data
                    allocate(lineData(nLinesInBin))
                    read(binUnit, iostat=ios) lineData
                    
                    if (ios == 0) then
                        lineCount = nLinesInBin
                        close(binUnit)
                        print *, 'Read', lineCount, 'spectral lines from binary cache file'
                        return
                    else
                        print *, 'WARNING: Error reading binary cache, falling back to text file'
                        deallocate(lineData)
                        close(binUnit)
                    end if
                else
                    print *, 'WARNING: Binary cache parameters do not match, falling back to text file'
                    close(binUnit)
                end if
            else
                print *, 'WARNING: Cannot open binary cache file, falling back to text file'
            end if
        end if
        
        ! Fall back to reading from text file
        print *, 'Reading from HITRAN text file and creating binary cache...'
        
        ! Initial allocation
        allocate(tempLines(lineCount))
        
        unit_num = 10
        open(unit_num, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR: Cannot open HITRAN file: ', trim(filename)
            stop 2
        end if
        
        currentLine = 0
        ios = 0
                
        ! Read lines and filter by wavenumber range
        do while (.not. is_iostat_end(ios))
            read(unit_num, parFormat, iostat=ios) MO, ISO, lineWV, refLineIntensityDP, A, &
                gammaForeign, gammaSelf, lineLowerStateDP, &
                foreignTempCoeff, deltaForeign, V1, V2, Q1, Q2, IERR, IREF, FLAG, Gu, Gl
            
            if (ios > 0) then
                print *, 'ERROR: Error reading HITRAN file at line', currentLine + 1
                stop 3
            end if
            
            if (ios < 0) exit  ! End of file
            
            if (currentLine>lineCount) exit
            
            ! Filter by wavenumber range
            if (lineWV >= extStartWV .and. lineWV <= extEndWV) then
                currentLine = currentLine + 1
                
                ! Store line data (simplified - assuming single isotope)
                tempLines(currentLine)%wv = lineWV
                tempLines(currentLine)%intensity = real(refLineIntensityDP)
                tempLines(currentLine)%gammaForeign = gammaForeign
                tempLines(currentLine)%gammaSelf = gammaSelf
                tempLines(currentLine)%lowerState = real(lineLowerStateDP)
                tempLines(currentLine)%tempCoeff = foreignTempCoeff
                tempLines(currentLine)%deltaForeign = deltaForeign
                tempLines(currentLine)%molIso = MO * 100 + convert_ISO(ISO)  ! Simplified encoding
            end if
        end do
        
        close(unit_num)
        
        ! Copy to final array
        lineCount = currentLine
        allocate(lineData(lineCount))
        lineData(1:lineCount) = tempLines(1:lineCount)
        
        deallocate(tempLines)
        
        ! Save binary cache for future use
        binUnit = 20
        open(binUnit, file=binFilename, form='unformatted', status='replace', iostat=ios)
        if (ios == 0) then
            ! Write cache header
            write(binUnit) lineCount, extStartWV, extEndWV, cutOff
            ! Write line data
            write(binUnit) lineData
            close(binUnit)
            print *, 'Created binary cache file: ', trim(binFilename)
        else
            print *, 'WARNING: Could not create binary cache file'
        end if
        
        print *, 'Read', lineCount, 'spectral lines from HITRAN file'
        
    end subroutine readHITRANFile
    
    
    function convert_ISO(ISOchar) result(ISOnum)
        implicit none

        ! Convert ISO character to number (from processParFile.f90)
        character(len=1), intent(in) :: ISOchar
        integer :: ISOnum
        
        select case (ISOchar)
            case ('1'); ISOnum = 1
            case ('2'); ISOnum = 2
            case ('3'); ISOnum = 3
            case ('4'); ISOnum = 4
            case ('5'); ISOnum = 5
            case ('6'); ISOnum = 6
            case ('7'); ISOnum = 7
            case ('8'); ISOnum = 8
            case ('9'); ISOnum = 9
            case ('0'); ISOnum = 10
            case ('A'); ISOnum = 11
            case default; ISOnum = 1  ! Default to most abundant
        end select
    end function convert_ISO
    
    
    subroutine calculateAbsorption(lines, nLines, shapeFuncPtr, chiFactorFuncPtr, moleculeIntCode, &
        startWV, endWV, deltaWV, cutOff, temperature, pSelf, pForeign, TIPS_ref, TIPS_target, EPS, &
        NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
        H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, nSubintervals)
        ! Main calculation routine
        ! Processes spectral interval in subintervals of width deltaWV
        implicit none

        type(SpectralLine), intent(in) :: lines(:) ! array of spectral lines
        integer, intent(in) :: nLines ! number of lines in the array
        
        real(kind=DP), intent(in) :: startWV, endWV, deltaWV, cutOff
        real(kind=DP), intent(in) :: temperature, pSelf, pForeign, TIPS_ref, TIPS_target

        real(kind=DP), intent(inout) :: EPS
        
        procedure(shape), pointer :: shapeFuncPtr ! pointer to the specific line shape function which will be used in grid calculations 
        procedure(chifactor), pointer :: chiFactorFuncPtr ! pointer to the specific chi factor function which will be used in grid calculations 
        
        integer, intent(in) :: moleculeIntCode ! integer code of the molecule
        
        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)
        real(kind=DP), intent(in) :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H
        
        integer, intent(in) :: nSubintervals
                
        real(kind=DP) :: startDeltaWV, endDeltaWV
        real(kind=DP) :: extStartDeltaWV, extEndDeltaWV
        integer :: lineIdx
        
        ! Timing variables (common)
        real(kind=DP) :: start_time, end_time, ctime, wtime ! clock time & wall time
        integer :: count_start, count_end, count_rate, count_max
        
        ! Timing variables (subprogram-specific)
        real elapsed_time1, wall_clock_time1 
        real elapsed_time2, wall_clock_time2 
        real elapsed_time3, wall_clock_time3
        real elapsed_time4, wall_clock_time4
        real elapsed_time5, wall_clock_time5
        
        integer :: iSubinterval
        
        integer :: startLineIdx
                
        ! Output header
        write(SPECTRA_FCHAN,'(A)') '# Wavenumber [cm-1]    Absorption cross-section [cm2/molec]'
        
        ! Process each subinterval
        startDeltaWV = startWV
        endDeltaWV = startDeltaWV + deltaWV
        lineIdx = 1
        
        elapsed_time1 = 0.0
        wall_clock_time1 = 0.0

        elapsed_time2 = 0.0
        wall_clock_time2 = 0.0

        elapsed_time3 = 0.0
        wall_clock_time3 = 0.0

        elapsed_time4 = 0.0
        wall_clock_time4 = 0.0

        elapsed_time5 = 0.0
        wall_clock_time5 = 0.0
                
        !do while (startDeltaWV < endWV)
        do iSubinterval = 1, nSubintervals
            ! Define extended subinterval boundaries
            extStartDeltaWV = startDeltaWV - cutOff
            extEndDeltaWV = endDeltaWV + cutOff
            
            if (DEBUG) write(DEBUG_FCHAN,"(i6,';')",advance="no") iSubinterval
            
            ! Reset grid arrays with timing            
            call system_clock(count_start, count_rate, count_max)
            call cpu_time(start_time)
            
            call resetGrids(EPS, NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
                RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
                RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P)

            call system_clock(count_end, count_rate, count_max)
            call cpu_time(end_time)
            
            ctime = end_time - start_time
            wtime = dble(count_end - count_start) / dble(count_rate)
                                                            
            elapsed_time1 = elapsed_time1 + ctime
            if (count_rate > 0) then
                wall_clock_time1 = wall_clock_time1 + wtime
            end if

            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") ctime
            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") wtime

            ! Find starting line for this subinterval with timing
            call system_clock(count_start, count_rate, count_max)
            call cpu_time(start_time)
            
            !print *,'lineIdx begin>>>',lineIdx
            !! Find starting line for this subinterval
            ! Ensure lineIdx is within valid bounds before the loop
            if (lineIdx < 1) lineIdx = 1
            if (lineIdx > nLines) lineIdx = nLines
            
            !print *, 'DEBUG: iSubinterval =', iSubinterval, 'lineIdx =', lineIdx, 'nLines =', nLines
            
            ! Safe loop: never access lines(0)
            do while (lineIdx > 1)
                if (lines(lineIdx-1)%wv >= extStartDeltaWV) then
                    lineIdx = lineIdx - 1
                else
                    exit
                end if
            end do

            call system_clock(count_end, count_rate, count_max)
            call cpu_time(end_time)
            
            ctime = end_time - start_time
            wtime = dble(count_end - count_start) / dble(count_rate)
                        
            elapsed_time2 = elapsed_time2 + ctime
            if (count_rate > 0) then
                wall_clock_time2 = wall_clock_time2 + wtime
            end if

            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") ctime
            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") wtime

            ! Process lines in this subinterval with timing
            call system_clock(count_start, count_rate, count_max)
            call cpu_time(start_time)
                        
            startLineIdx = lineIdx
            call processSubinterval(startDeltaWV, endDeltaWV, lineIdx, shapeFuncPtr, chiFactorFuncPtr, &
                moleculeIntCode, lines, nLines, cutOff, temperature, pSelf, pForeign, TIPS_ref, TIPS_target, EPS, &
                NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
                RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, RK5, &
                RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
                H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, deltaWV)
            if (lineIdx < 1) lineIdx = 1
            if (lineIdx > nLines) lineIdx = nLines
                
            call system_clock(count_end, count_rate, count_max)
            call cpu_time(end_time)
            
            ctime = end_time - start_time
            wtime = dble(count_end - count_start) / dble(count_rate)
                        
            elapsed_time3 = elapsed_time3 + ctime
            if (count_rate > 0) then
                wall_clock_time3 = wall_clock_time3 + wtime
            end if

            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") ctime
            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") wtime

            !Cascade interpolation with timing
            call system_clock(count_start, count_rate, count_max)
            call cpu_time(start_time)
                                    
            call cascadeInterpolation(NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
                RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
                RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P)
                
            call system_clock(count_end, count_rate, count_max)
            call cpu_time(end_time)
            
            ctime = end_time - start_time
            wtime = dble(count_end - count_start) / dble(count_rate)
                        
            elapsed_time4 = elapsed_time4 + ctime
            if (count_rate > 0) then
                wall_clock_time4 = wall_clock_time4 + wtime
            end if

            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") ctime
            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") wtime
            
            ! Output results for this subinterval with timing
            call system_clock(count_start, count_rate, count_max)
            call cpu_time(start_time)

            if (OUTPUT_SPECTRA) call outputSubintervalResults(startDeltaWV, endDeltaWV, NT, RK, H)

            call system_clock(count_end, count_rate, count_max)
            call cpu_time(end_time)

            ctime = end_time - start_time
            wtime = dble(count_end - count_start) / dble(count_rate)
            
            elapsed_time5 = elapsed_time5 + ctime
            if (count_rate > 0) then
                wall_clock_time5 = wall_clock_time5 + wtime
            end if
            
            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") ctime
            if (DEBUG) write(DEBUG_FCHAN,"(E12.4,';')",advance="no") wtime
            
            if (DEBUG) then
                write(DEBUG_FCHAN,"(E12.4,';')",advance="no") startDeltaWV
                write(DEBUG_FCHAN,"(E12.4,';')",advance="no") endDeltaWV
                write(DEBUG_FCHAN,"(E12.4,';')",advance="no") deltaWV
                write(DEBUG_FCHAN,"(I10,';')",advance="no") startLineIdx
                write(DEBUG_FCHAN,"(I10)") lineIdx
            end if 
            
            ! Move to next subinterval
            startDeltaWV = startDeltaWV + deltaWV
            endDeltaWV = endDeltaWV + deltaWV
        end do
        
        print *, '=== TIMING: resetGrids (total) ==='
        print *, 'CPU time: ', elapsed_time1, ' seconds'
        print *, 'Wall clock time: ', wall_clock_time1, ' seconds'        

        print *, '=== TIMING: finding starting line (total) ==='
        print *, 'CPU time: ', elapsed_time2, ' seconds'
        print *, 'Wall clock time: ', wall_clock_time2, ' seconds'        

        print *, '=== TIMING: processSubinterval (total) ==='
        print *, 'CPU time: ', elapsed_time3, ' seconds'
        print *, 'Wall clock time: ', wall_clock_time3, ' seconds'        

        print *, '=== TIMING: cascadeInterpolation (total) ==='
        print *, 'CPU time: ', elapsed_time4, ' seconds'
        print *, 'Wall clock time: ', wall_clock_time4, ' seconds'        

        print *, '=== TIMING: outputSubintervalResults (total) ==='
        print *, 'CPU time: ', elapsed_time5, ' seconds'
        print *, 'Wall clock time: ', wall_clock_time5, ' seconds'        
        
    end subroutine calculateAbsorption
    
    
    subroutine processSubinterval(startDelta, endDelta, lineStartIdx, shapeFuncPtr, chiFactorFuncPtr, &
        moleculeIntCode, lines, nLines, cutOff, temperature, pSelf, pForeign, TIPS_ref, TIPS_target, EPS, &
        NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
        H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, deltaWV)
        ! Process all lines contributing to current subinterval
        implicit none
        
        real(kind=DP), intent(in) :: startDelta, endDelta, cutOff
        real(kind=DP), intent(in) :: temperature, pSelf, pForeign, TIPS_ref, TIPS_target, EPS
        
        integer, intent(inout) :: lineStartIdx

        type(SpectralLine), intent(in) :: lines(:) ! array of spectral lines
        integer, intent(in) :: nLines ! number of lines in the array
        
        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)
        real(kind=DP), intent(in) :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H

        real(kind=DP), intent(in) :: deltaWV

        integer, intent(in) :: moleculeIntCode

        procedure(shape), pointer :: shapeFuncPtr
        procedure(chiFactor), pointer :: chiFactorFuncPtr

        real(kind=DP) :: extStartDelta, extEndDelta
        real(kind=DP) :: shiftedWV
        real(kind=DP) :: dopHWHM, lorHWHM, lineIntensity
        integer :: idx
        
        extStartDelta = startDelta - cutOff
        extEndDelta = endDelta + cutOff
        
        ! Process each line
        idx = lineStartIdx
        do while (idx <= nLines)
            if (lines(idx)%wv > extEndDelta) exit
            
            ! Calculate line parameters at target conditions
            shiftedWV = lines(idx)%wv + lines(idx)%deltaForeign * (pSelf + pForeign)
            
            ! Skip if outside extended range
            if (shiftedWV < extStartDelta .or. shiftedWV > extEndDelta) then
                idx = idx + 1
                cycle
            end if
            
            ! Calculate line shape parameters
            dopHWHM = calculateDopplerHWHM(shiftedWV, temperature, lines(idx)%molIso)
            lorHWHM = calculateLorentzHWHM(lines(idx), temperature, pSelf, pForeign)
            lineIntensity = calculateIntensity(lines(idx), temperature, TIPS_ref, TIPS_target)
            
            ! Calculate contribution on grid
            if (shiftedWV < startDelta) then
                call calculateLeftWing(startDelta, shiftedWV, shapeFuncPtr, EPS, dopHWHM, lorHWHM, &
                    lineIntensity, chiFactorFuncPtr, moleculeIntCode, &
                    NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
                    RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
                    RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
                    H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, cutOff)
            else if (shiftedWV > endDelta) then
                call calculateRightWing(startDelta, shiftedWV, shapeFuncPtr, EPS, dopHWHM, lorHWHM, &
                    lineIntensity, chiFactorFuncPtr, moleculeIntCode, &
                    NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
                    RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
                    RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
                    H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, cutOff, deltaWV)
            else
                call calculateCenter(startDelta, shiftedWV, shapeFuncPtr, EPS, dopHWHM, lorHWHM, lineIntensity, &
                    chiFactorFuncPtr, moleculeIntCode, &
                    NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
                    RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, RK5, &
                    RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
                    H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, deltaWV)
            end if
            
            idx = idx + 1
        end do
        
        ! Update line index for next subinterval
        lineStartIdx = idx
        
        ! Ensure lineStartIdx is within valid bounds
        if (lineStartIdx < 1) lineStartIdx = 1
        if (lineStartIdx > nLines) lineStartIdx = nLines
        
    end subroutine processSubinterval
    
    
    function calculateDopplerHWHM(wv, temp, molIso) result(dopHWHM)
        ! Calculate Doppler half-width at half maximum
        implicit none

        real(kind=DP), intent(in) :: wv
        real(kind=DP), intent(in) :: temp
        integer, intent(in) :: molIso
        real :: dopHWHM
        
        real :: molarMass
        
        ! Get molar mass (simplified - using approximate values)
        molarMass = getMolarMass(molIso / 100)
        
        dopHWHM = dopplerCONST * wv * sqrt(temp / molarMass)
        
    end function calculateDopplerHWHM
    
    
    function calculateLorentzHWHM(line, temp, pS, pF) result(lorHWHM)
        ! Calculate pressure-broadened Lorentz half-width
        implicit none
        
        type(SpectralLine), intent(in) :: line
        real(kind=DP), intent(in) :: temp, pS, pF
        real :: lorHWHM
        
        lorHWHM = ((refTemperature / temp)**line%tempCoeff) * &
                  (line%gammaForeign * pF + line%gammaSelf * pS)
        
    end function calculateLorentzHWHM
    
    
    function calculateIntensity(line, temp, TIPS_ref, TIPS_target) result(intensity)
        ! Calculate line intensity at target temperature
        implicit none
        
        type(SpectralLine), intent(in) :: line
        real(kind=DP), intent(in) :: temp, TIPS_ref, TIPS_target
        real(kind=DP) :: intensity
        
        real(kind=DP) :: TIPSFactor, boltzmannFactor, stimulatedEmissionFactor
        
        ! TIPS factor (using provided partition sums)
        TIPSFactor = TIPS_ref / TIPS_target
        
        ! Boltzmann factor
        boltzmannFactor = exp(-C2 * line%lowerState / temp) / &
                         exp(-C2 * line%lowerState / refTemperature)
        
        ! Stimulated emission correction
        stimulatedEmissionFactor = (1.0 - exp(-C2 * line%wv / temp)) / &
                                  (1.0 - exp(-C2 * line%wv / refTemperature))
        
        intensity = line%intensity * TIPSFactor * boltzmannFactor * stimulatedEmissionFactor
        
    end function calculateIntensity
    
    
    function getMolarMass(molNum) result(mass)
        ! Get molar mass for molecule (simplified)
        implicit none

        integer, intent(in) :: molNum
        real :: mass
        
        select case (molNum)
            case (1); mass = 18.015  ! H2O
            case (2); mass = 44.01   ! CO2
            case (3); mass = 47.998  ! O3
            case (4); mass = 44.013  ! N2O
            case (5); mass = 28.01   ! CO
            case (6); mass = 16.04   ! CH4
            case (7); mass = 31.999  ! O2
            case (8); mass = 30.01   ! NO
            case (9); mass = 64.066  ! SO2
            case default; mass = 44.01  ! Default to CO2
        end select
        
    end function getMolarMass
    
    
    subroutine resetGrids(EPS, NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, RK5, RK5L, &
        RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P)
        ! Reset all grid arrays to zero
        implicit none

        real(kind=DP), intent(inout) :: EPS

        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)

        RK = 0.0
        RK0 = 0.0; RK0L = 0.0; RK0P = 0.0
        RK1 = 0.0; RK1L = 0.0; RK1P = 0.0
        RK2 = 0.0; RK2L = 0.0; RK2P = 0.0
        RK3 = 0.0; RK3L = 0.0; RK3P = 0.0
        RK4 = 0.0; RK4L = 0.0; RK4P = 0.0
        RK5 = 0.0; RK5L = 0.0; RK5P = 0.0
        RK6 = 0.0; RK6L = 0.0; RK6P = 0.0
        RK7 = 0.0; RK7L = 0.0; RK7P = 0.0
        RK8 = 0.0; RK8L = 0.0; RK8P = 0.0
        RK9 = 0.0; RK9L = 0.0; RK9P = 0.0
        EPS = 0.0D0
        
    end subroutine resetGrids    


    subroutine calculateLeftWing(FREQ, UL, FSHAPE, EPS, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode, &
        NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
        H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, cutOff)
        ! calculation contributions to the subinterval on each grid
        ! FROM the LEFT PART of the extended subinterval: [startDeltaWV-cutOff; startDeltaWV]

        ! TODO: when refactoring this file, deal with EPS. In legacy it is under `save` 

        !IMPLICIT INTEGER*4 (I-N)
        implicit none
        
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL, UU
        real(kind=DP) :: EPS
        real(kind=DP) :: XXX
        real :: FF
        integer :: I
        
        procedure(shape), pointer :: FSHAPE
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode
        real(kind=DP) :: lineIntensity
        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)
        real(kind=DP), intent(in) :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H
        real(kind=DP), intent(in) :: dopHWHM, lorHWHM, cutOff
        
        UU = UL - FREQ
        
        ! Early return conditions
        if (UU >= 0.0_DP) return
        if (-UU > cutOff) return
        
        FF = FSHAPE(UU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        if (FF < EPS) return
        
        ! Add contribution to main grid
        RK(1) = FF + RK(1)
        
        ! Branch based on distance from grid boundary
        if (-UU < H0) then
            ! Case: Very close to grid boundary (within H0)
            ! Continue to main processing loop
        else if (-UU < H1) then
            ! Case: Within H1 distance from boundary
            RK1P(1) = RK1P(1) + FF
            FF = FSHAPE(UU - H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1(1) = RK1(1) + FF
            FF = FSHAPE(UU - H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1L(1) = RK1L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H2) then
            ! Case: Within H2 distance from boundary
            RK2P(1) = RK2P(1) + FF
            FF = FSHAPE(UU - H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2(1) = RK2(1) + FF
            FF = FSHAPE(UU - H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2L(1) = RK2L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H3) then
            ! Case: Within H3 distance from boundary
            RK3P(1) = RK3P(1) + FF
            FF = FSHAPE(UU - H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3(1) = RK3(1) + FF
            FF = FSHAPE(UU - H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3L(1) = RK3L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H4) then
            ! Case: Within H4 distance from boundary
            RK4P(1) = RK4P(1) + FF
            FF = FSHAPE(UU - H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4(1) = RK4(1) + FF
            FF = FSHAPE(UU - H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4L(1) = RK4L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H5) then
            ! Case: Within H5 distance from boundary
            RK5P(1) = RK5P(1) + FF
            FF = FSHAPE(UU - H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5(1) = RK5(1) + FF
            FF = FSHAPE(UU - H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5L(1) = RK5L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H6) then
            ! Case: Within H6 distance from boundary
            RK6P(1) = RK6P(1) + FF
            FF = FSHAPE(UU - H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6(1) = RK6(1) + FF
            FF = FSHAPE(UU - H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6L(1) = RK6L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H7) then
            ! Case: Within H7 distance from boundary
            RK7P(1) = RK7P(1) + FF
            FF = FSHAPE(UU - H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7(1) = RK7(1) + FF
            FF = FSHAPE(UU - H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7L(1) = RK7L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H8) then
            ! Case: Within H8 distance from boundary
            RK8P(1) = RK8P(1) + FF
            FF = FSHAPE(UU - H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8(1) = RK8(1) + FF
            FF = FSHAPE(UU - H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8L(1) = RK8L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else if (-UU < H9) then
            ! Case: Within H9 distance from boundary
            RK9P(1) = RK9P(1) + FF
            FF = FSHAPE(UU - H - H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9(1) = RK9(1) + FF
            FF = FSHAPE(UU - H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9L(1) = RK9L(1) + FF
            if (FF < EPS) return
            ! Continue to cascade processing
        else
            ! Case: Beyond H9 distance - add to main grid points
            RK(2) = RK(2) + FSHAPE(UU - H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK(3) = RK(3) + FSHAPE(UU - H - H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK(4) = RK(4) + FSHAPE(UU + H - H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            FF = FSHAPE(UU - H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK(5) = RK(5) + FF
        end if
        
        ! Cascade processing - propagate contributions to higher grid levels
        if (-UU >= H9) then
            ! Cascade level 9
            RK9P(2) = RK9P(2) + FF
            FF = FSHAPE(UU - H9 - H - H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9(2) = RK9(2) + FF
            FF = FSHAPE(UU - H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9L(2) = RK9L(2) + FF
        end if
        
        if (-UU >= H8) then
            ! Cascade level 8
            RK8P(2) = RK8P(2) + FF
            FF = FSHAPE(UU - H8 - H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8(2) = RK8(2) + FF
            FF = FSHAPE(UU - H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8L(2) = RK8L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H7) then
            ! Cascade level 7
            RK7P(2) = RK7P(2) + FF
            FF = FSHAPE(UU - H7 - H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7(2) = RK7(2) + FF
            FF = FSHAPE(UU - H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7L(2) = RK7L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H6) then
            ! Cascade level 6
            RK6P(2) = RK6P(2) + FF
            FF = FSHAPE(UU - H6 - H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6(2) = RK6(2) + FF
            FF = FSHAPE(UU - H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6L(2) = RK6L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H5) then
            ! Cascade level 5
            RK5P(2) = RK5P(2) + FF
            FF = FSHAPE(UU - H5 - H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5(2) = RK5(2) + FF
            FF = FSHAPE(UU - H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5L(2) = RK5L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H4) then
            ! Cascade level 4
            RK4P(2) = RK4P(2) + FF
            FF = FSHAPE(UU - H4 - H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4(2) = RK4(2) + FF
            FF = FSHAPE(UU - H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4L(2) = RK4L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H3) then
            ! Cascade level 3
            RK3P(2) = RK3P(2) + FF
            FF = FSHAPE(UU - H3 - H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3(2) = RK3(2) + FF
            FF = FSHAPE(UU - H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3L(2) = RK3L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H2) then
            ! Cascade level 2
            RK2P(2) = RK2P(2) + FF
            FF = FSHAPE(UU - H2 - H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2(2) = RK2(2) + FF
            FF = FSHAPE(UU - H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2L(2) = RK2L(2) + FF
            if (FF < EPS) return
        end if
        
        if (-UU >= H1) then
            ! Cascade level 1
            RK1P(2) = RK1P(2) + FF
            FF = FSHAPE(UU - H1 - H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1(2) = RK1(2) + FF
            FF = FSHAPE(UU - H0, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1L(2) = RK1L(2) + FF
        end if
        
        ! Main processing loop for grid level 0
        if (-UU >= H0) then
            RK0P(1) = RK0P(1) + FF
            FF = FSHAPE(UU - H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0(1) = RK0(1) + FF
            FF = FSHAPE(UU - H0, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0L(1) = RK0L(1) + FF
        end if
        
        ! Continue with remaining grid points
        XXX = H0
        do I = 2, NT0
            RK0P(I) = RK0P(I) + FF
            FF = FSHAPE(UU - XXX - H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0(I) = RK0(I) + FF
            XXX = XXX + H0
            FF = FSHAPE(UU - XXX, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0L(I) = RK0L(I) + FF
            if (FF < EPS) exit
        end do
        
    end subroutine calculateLeftWing    
    
    
    subroutine calculateRightWing(FREQ, UL, FSHAPE, EPS, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, &
        moleculeIntCode, NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
        H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, cutOff, deltaWV)
        ! calculation contributions to the subinterval on each grid
        ! FROM the RIGHT PART of the extended subinterval: [endDeltaWV; endDeltaWV+cutOff]

        !IMPLICIT INTEGER*4 (I-N)
        implicit none
        
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL, UU
        real(kind=DP) :: EPS
        real(kind=DP) :: XXX
        real :: FF
        integer :: I, N, N2
        
        procedure(shape), pointer :: FSHAPE
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode
        real(kind=DP) :: lineIntensity
        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)
        real(kind=DP), intent(in) :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H
        real(kind=DP), intent(in) :: cutOff
        real(kind=DP), intent(in) :: dopHWHM, lorHWHM, deltaWV

        UU = UL - FREQ - deltaWV ! CHANGE DIAP
        
        ! Early return conditions
        if (UU >= cutOff) return
        FF = FSHAPE(UU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        if (FF < EPS) return
        
        ! Branch based on distance from grid boundary
        if (UU >= H0) then
            ! Case: Beyond H0 distance from boundary
            RK0L(NT0) = RK0L(NT0) + FF
            FF = FSHAPE(UU + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0(NT0) = RK0(NT0) + FF
            FF = FSHAPE(UU + H0, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0P(NT0) = RK0P(NT0) + FF
        else if (UU >= H1) then
            ! Case: Within H1 distance from boundary
            RK1L(NT1) = RK1L(NT1) + FF
            FF = FSHAPE(UU + H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1(NT1) = RK1(NT1) + FF
            FF = FSHAPE(UU + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1P(NT1) = RK1P(NT1) + FF
            if (FF < EPS) return
        else if (UU >= H2) then
            ! Case: Within H2 distance from boundary
            RK2L(NT2) = RK2L(NT2) + FF
            FF = FSHAPE(UU + H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2(NT2) = RK2(NT2) + FF
            FF = FSHAPE(UU + H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2P(NT2) = RK2P(NT2) + FF
            if (FF < EPS) return
        else if (UU >= H3) then
            ! Case: Within H3 distance from boundary
            RK3L(NT3) = RK3L(NT3) + FF
            FF = FSHAPE(UU + H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3(NT3) = RK3(NT3) + FF
            FF = FSHAPE(UU + H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3P(NT3) = RK3P(NT3) + FF
            if (FF < EPS) return
        else if (UU >= H4) then
            ! Case: Within H4 distance from boundary
            RK4L(NT4) = RK4L(NT4) + FF
            FF = FSHAPE(UU + H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4(NT4) = RK4(NT4) + FF
            FF = FSHAPE(UU + H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4P(NT4) = RK4P(NT4) + FF
            if (FF < EPS) return
        else if (UU >= H5) then
            ! Case: Within H5 distance from boundary
            RK5L(NT5) = RK5L(NT5) + FF
            FF = FSHAPE(UU + H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5(NT5) = RK5(NT5) + FF
            FF = FSHAPE(UU + H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5P(NT5) = RK5P(NT5) + FF
            if (FF < EPS) return
        else if (UU >= H6) then
            ! Case: Within H6 distance from boundary
            RK6L(NT6) = RK6L(NT6) + FF
            FF = FSHAPE(UU + H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6(NT6) = RK6(NT6) + FF
            FF = FSHAPE(UU + H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6P(NT6) = RK6P(NT6) + FF
            if (FF < EPS) return
        else if (UU >= H7) then
            ! Case: Within H7 distance from boundary
            RK7L(NT7) = RK7L(NT7) + FF
            FF = FSHAPE(UU + H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7(NT7) = RK7(NT7) + FF
            FF = FSHAPE(UU + H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7P(NT7) = RK7P(NT7) + FF
            if (FF < EPS) return
        else if (UU >= H8) then
            ! Case: Within H8 distance from boundary
            RK8L(NT8) = RK8L(NT8) + FF
            FF = FSHAPE(UU + H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8(NT8) = RK8(NT8) + FF
            FF = FSHAPE(UU + H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8P(NT8) = RK8P(NT8) + FF
            if (FF < EPS) return
        else if (UU >= H9) then
            ! Case: Within H9 distance from boundary
            RK9L(NT9) = RK9L(NT9) + FF
            FF = FSHAPE(UU + H + H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9(NT9) = RK9(NT9) + FF
            FF = FSHAPE(UU + H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9P(NT9) = RK9P(NT9) + FF
            if (FF < EPS) return
        else
            ! Case: Beyond H9 distance - add to main grid points
            RK(NT) = RK(NT) + FSHAPE(UU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK(NT - 1) = RK(NT - 1) + FSHAPE(UU + H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK(NT - 2) = RK(NT - 2) + FSHAPE(UU + H + H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK(NT - 3) = RK(NT - 3) + FSHAPE(UU + H9 - H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        end if
        
        ! Cascade processing - propagate contributions to higher grid levels
        if (UU >= H9) then
            ! Cascade level 9
            N2 = NT9 - 1
            FF = FSHAPE(UU + H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9L(N2) = RK9L(N2) + FF
            FF = FSHAPE(UU + H9 + H + H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9(N2) = RK9(N2) + FF
            FF = FSHAPE(UU + H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK9P(N2) = RK9P(N2) + FF
        end if
        
        if (UU >= H8) then
            ! Cascade level 8
            N = NT8 - 1
            RK8L(N) = RK8L(N) + FF
            FF = FSHAPE(UU + H8 + H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8(N) = RK8(N) + FF
            FF = FSHAPE(UU + H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK8P(N) = RK8P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H7) then
            ! Cascade level 7
            N = NT7 - 1
            RK7L(N) = RK7L(N) + FF
            FF = FSHAPE(UU + H7 + H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7(N) = RK7(N) + FF
            FF = FSHAPE(UU + H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK7P(N) = RK7P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H6) then
            ! Cascade level 6
            N = NT6 - 1
            RK6L(N) = RK6L(N) + FF
            FF = FSHAPE(UU + H6 + H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6(N) = RK6(N) + FF
            FF = FSHAPE(UU + H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK6P(N) = RK6P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H5) then
            ! Cascade level 5
            N = NT5 - 1
            RK5L(N) = RK5L(N) + FF
            FF = FSHAPE(UU + H5 + H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5(N) = RK5(N) + FF
            FF = FSHAPE(UU + H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK5P(N) = RK5P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H4) then
            ! Cascade level 4
            N = NT4 - 1
            RK4L(N) = RK4L(N) + FF
            FF = FSHAPE(UU + H4 + H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4(N) = RK4(N) + FF
            FF = FSHAPE(UU + H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK4P(N) = RK4P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H3) then
            ! Cascade level 3
            N = NT3 - 1
            RK3L(N) = RK3L(N) + FF
            FF = FSHAPE(UU + H3 + H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3(N) = RK3(N) + FF
            FF = FSHAPE(UU + H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK3P(N) = RK3P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H2) then
            ! Cascade level 2
            N = NT2 - 1
            RK2L(N) = RK2L(N) + FF
            FF = FSHAPE(UU + H2 + H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2(N) = RK2(N) + FF
            FF = FSHAPE(UU + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK2P(N) = RK2P(N) + FF
            if (FF < EPS) return
        end if
        
        if (UU >= H1) then
            ! Cascade level 1
            N = NT1 - 1
            RK1L(N) = RK1L(N) + FF
            FF = FSHAPE(UU + H1 + H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1(N) = RK1(N) + FF
            FF = FSHAPE(UU + H0, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK1P(N) = RK1P(N) + FF
        end if
        
        ! Main processing loop for grid level 0
        if (UU >= H0) then
            RK0L(NT0) = RK0L(NT0) + FF
            FF = FSHAPE(UU + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0(NT0) = RK0(NT0) + FF
            FF = FSHAPE(UU + H0, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0P(NT0) = RK0P(NT0) + FF
        end if
        
        ! Continue with remaining grid points
        XXX = H0
        do I = NT0 - 1, 1, -1
            RK0L(I) = RK0L(I) + FF
            FF = FSHAPE(UU + XXX + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0(I) = RK0(I) + FF
            XXX = XXX + H0
            FF = FSHAPE(UU + XXX, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            RK0P(I) = RK0P(I) + FF
            if (FF < EPS) exit
        end do
        
        RK(1) = RK(1) + FF
        
    end subroutine calculateRightWing
    

    subroutine calculateCenter(FREQ, UL, FSHAPE, EPS, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, &
        moleculeIntCode, NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P, &
        H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H, deltaWV)
        ! calculation contributions to the subinterval on each grid
        ! FROM the CENTRAL PART of the extended subinterval: [startDeltaWV; endDeltaWV]
        
        !IMPLICIT INTEGER*4 (I-N)
        implicit none
        
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL, UU, CONSER, UUU
        real(kind=DP) :: EPS, EPS4
        real :: FF, FA
        integer :: I, IB, ICON, II, III, NPOINT
        
        procedure(shape), pointer :: FSHAPE
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode
        real(kind=DP) :: lineIntensity
        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)
        real(kind=DP), intent(in) :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H
        
        real(kind=DP), intent(in) :: dopHWHM, lorHWHM, deltaWV
        
        UU = UL - FREQ
        
        ! Early return if outside central region
        if (UU >= deltaWV) return
        
        FF = FSHAPE(0.0_DP, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        if (FF < EPS) return
        
        ! Initialize processing
        NPOINT = 1
        CONSER = UU - H
        FA = FSHAPE(UU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        EPS4 = EPS * 0.25_DP
        
        if (FA > EPS4) RK(1) = RK(1) + FA
        
        ! Process left-right side
        if (UU >= H) then
            I = 0
            UUU = UU
            
            ! Process grid level 0
            if (UUU >= H0 + H0) then
                do I = 1, NT0
                    UUU = UUU - H0
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK0P(I) = RK0P(I) + FA
                    RK0(I) = RK0(I) + FSHAPE(UUU + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK0L(I) = RK0L(I) + FF
                    FA = FF
                    if (UUU - H0 < H0) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 1
            if (UUU >= H0) then
                IB = I + 1
                do I = IB, NT1
                    UUU = UUU - H1
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK1P(I) = RK1P(I) + FA
                    RK1(I) = RK1(I) + FSHAPE(UUU + H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK1L(I) = RK1L(I) + FF
                    FA = FF
                    if (UUU - H1 < H1) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 2
            if (UUU >= H1) then
                IB = I + 1
                do I = IB, NT2
                    UUU = UUU - H2
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK2P(I) = RK2P(I) + FA
                    RK2(I) = RK2(I) + FSHAPE(UUU + H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK2L(I) = RK2L(I) + FF
                    FA = FF
                    if (UUU - H2 < H2) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 3
            if (UUU >= H2) then
                IB = I + 1
                do I = IB, NT3
                    UUU = UUU - H3
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK3P(I) = RK3P(I) + FA
                    RK3(I) = RK3(I) + FSHAPE(UUU + H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK3L(I) = RK3L(I) + FF
                    FA = FF
                    if (UUU - H3 < H3) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 4
            if (UUU >= H3) then
                IB = I + 1
                do I = IB, NT4
                    UUU = UUU - H4
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK4P(I) = RK4P(I) + FA
                    RK4(I) = RK4(I) + FSHAPE(UUU + H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK4L(I) = RK4L(I) + FF
                    FA = FF
                    if (UUU - H4 < H4) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 5
            if (UUU >= H4) then
                IB = I + 1
                do I = IB, NT5
                    UUU = UUU - H5
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK5P(I) = RK5P(I) + FA
                    RK5(I) = RK5(I) + FSHAPE(UUU + H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK5L(I) = RK5L(I) + FF
                    FA = FF
                    if (UUU - H5 < H5) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 6
            if (UUU >= H5) then
                IB = I + 1
                do I = IB, NT6
                    UUU = UUU - H6
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK6P(I) = RK6P(I) + FA
                    RK6(I) = RK6(I) + FSHAPE(UUU + H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK6L(I) = RK6L(I) + FF
                    FA = FF
                    if (UUU - H6 < H6) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 7
            if (UUU >= H6) then
                IB = I + 1
                do I = IB, NT7
                    UUU = UUU - H7
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK7P(I) = RK7P(I) + FA
                    RK7(I) = RK7(I) + FSHAPE(UUU + H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK7L(I) = RK7L(I) + FF
                    FA = FF
                    if (UUU - H7 < H7) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 8
            if (UUU >= H7) then
                IB = I + 1
                do I = IB, NT8
                    UUU = UUU - H8
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK8P(I) = RK8P(I) + FA
                    RK8(I) = RK8(I) + FSHAPE(UUU + H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK8L(I) = RK8L(I) + FF
                    FA = FF
                    if (UUU - H8 < H8) exit
                end do
            end if
            
            I = I * 2
            
            ! Process grid level 9
            if (UUU >= H8) then
                IB = I + 1
                do I = IB, NT9
                    UUU = UUU - H9
                    FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    if (FF < EPS) exit
                    RK9P(I) = RK9P(I) + FA
                    RK9(I) = RK9(I) + FSHAPE(UUU + H + H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                    RK9L(I) = RK9L(I) + FF
                    FA = FF
                    if (UUU - H9 < H9) exit
                end do
            end if
            
            I = I * 4
            IB = I + 2
            CONSER = UU - (IB - 1) * H
            
            do ICON = IB, NT
                RK(ICON) = RK(ICON) + FSHAPE(CONSER, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                CONSER = CONSER - H
                if (CONSER < 0.0_DP) then
                    NPOINT = ICON
                    exit
                end if
            end do
        end if
        
        ! Process right-left side
        NPOINT = NPOINT + 1
        UUU = deltaWV - UU
        FA = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        III = 0
        
        ! Process grid level 0 (reverse)
        if (UUU >= H0 + H0) then
            do I = NT0, 1, -1
                III = III + 1
                UUU = UUU - H0
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK0L(I) = RK0L(I) + FA
                RK0(I) = RK0(I) + FSHAPE(UUU + H1, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK0P(I) = RK0P(I) + FF
                FA = FF
                if (UUU - H0 < H0) exit
            end do
        end if
        
        if (UUU >= H0) then
            III = III * 2
            IB = NT1 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H1
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK1L(I) = RK1L(I) + FA
                RK1(I) = RK1(I) + FSHAPE(UUU + H2, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK1P(I) = RK1P(I) + FF
                FA = FF
                if (UUU - H1 < H1) exit
            end do
        end if
        
        ! Process grid level 2 (reverse)
        if (UUU >= H1) then
            III = III * 2
            IB = NT2 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H2
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK2L(I) = RK2L(I) + FA
                RK2(I) = RK2(I) + FSHAPE(UUU + H3, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK2P(I) = RK2P(I) + FF
                FA = FF
                if (UUU - H2 < H2) exit
            end do
        end if
        
        ! Process grid level 3 (reverse)
        if (UUU >= H2) then
            III = III * 2
            IB = NT3 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H3
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK3L(I) = RK3L(I) + FA
                RK3(I) = RK3(I) + FSHAPE(UUU + H4, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK3P(I) = RK3P(I) + FF
                FA = FF
                if (UUU - H3 < H3) exit
            end do
        end if
        
        ! Process grid level 4 (reverse)
        if (UUU >= H3) then
            III = III * 2
            IB = NT4 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H4
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK4L(I) = RK4L(I) + FA
                RK4(I) = RK4(I) + FSHAPE(UUU + H5, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK4P(I) = RK4P(I) + FF
                FA = FF
                if (UUU - H4 < H4) exit
            end do
        end if
        
        ! Process grid level 5 (reverse)
        if (UUU >= H4) then
            III = III * 2
            IB = NT5 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H5
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK5L(I) = RK5L(I) + FA
                RK5(I) = RK5(I) + FSHAPE(UUU + H6, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK5P(I) = RK5P(I) + FF
                FA = FF
                if (UUU - H5 < H5) exit
            end do
        end if
        
        ! Process grid level 6 (reverse)
        if (UUU >= H5) then
            III = III * 2
            IB = NT6 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H6
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK6L(I) = RK6L(I) + FA
                RK6(I) = RK6(I) + FSHAPE(UUU + H7, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK6P(I) = RK6P(I) + FF
                FA = FF
                if (UUU - H6 < H6) exit
            end do
        end if
        
        ! Process grid level 7 (reverse)
        if (UUU >= H6) then
            III = III * 2
            IB = NT7 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H7
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK7L(I) = RK7L(I) + FA
                RK7(I) = RK7(I) + FSHAPE(UUU + H8, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK7P(I) = RK7P(I) + FF
                FA = FF
                if (UUU - H7 < H7) exit
            end do
        end if
        
        ! Process grid level 8 (reverse)
        if (UUU >= H7) then
            III = III * 2
            IB = NT8 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H8
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK8L(I) = RK8L(I) + FA
                RK8(I) = RK8(I) + FSHAPE(UUU + H9, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK8P(I) = RK8P(I) + FF
                FA = FF
                if (UUU - H8 < H8) exit
            end do
        end if
        
        ! Process grid level 9 (reverse)
        if (UUU >= H8) then
            III = III * 2
            IB = NT9 - III
            do I = IB, 1, -1
                III = III + 1
                UUU = UUU - H9
                FF = FSHAPE(UUU, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                if (FF < EPS) exit
                RK9L(I) = RK9L(I) + FA
                RK9(I) = RK9(I) + FSHAPE(UUU + H + H, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
                RK9P(I) = RK9P(I) + FF
                FA = FF
                if (UUU - H9 < H9) exit
            end do
        end if
        
        ! Final processing for main grid
        III = III * 4
        I = NT - III
        do II = NPOINT, I
            RK(II) = RK(II) + FSHAPE(CONSER, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            CONSER = CONSER - H
        end do
        
    end subroutine calculateCenter

! Content of Spectroscopy.f90

    pure function lorentz(X, lorHWHM, lineIntensity) result(lorentzShape)
        implicit none

        real(kind=DP) :: lorentzShape
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: lorHWHM ! doppler HWHM
        real(kind=DP), intent(in) :: lineIntensity

        lorentzShape = lorHWHM / (pi*(X**2 + lorHWHM**2)) * lineIntensity
    end function lorentz

    
    pure function doppler(X, dopHWHM, lineIntensity) result(dopplerShape)
        implicit none

        real(kind=DP) :: dopplerShape
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: dopHWHM ! doppler HWHM 
        real(kind=DP), intent(in) :: lineIntensity

        dopplerShape = sqln2 / (sqrt(pi) * dopHWHM) * exp(-(X/dopHWHM)**2 * log(2.)) * lineIntensity
    end function doppler

    
    pure function chiCorrectedLorentz(X, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode) result(chiCorrectedLorentzShape)
        ! Lorentz line shape with a χ-corrected wing
        implicit none

        real(kind=DP) :: chiCorrectedLorentzShape
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: lorHWHM ! Lorentz HWHM
        real(kind=DP), intent(in) :: lineIntensity
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode

        chiCorrectedLorentzShape = lorentz(X, lorHWHM, lineIntensity) * chiFactorFuncPtr(X, moleculeIntCode)

    end function chiCorrectedLorentz

    
    pure function voigtAsymptotic1(X, lorHWHM, VX, lineIntensity, chiFactorFuncPtr, moleculeIntCode) result(voigtAsymptotic1Value)
        ! Lorentz leading, Doppler-influenced asymptotic correction for the Voigt function (y<<1, x>>1)
        implicit none

        real(kind=DP) :: voigtAsymptotic1Value
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: lorHWHM
        real(kind=DP), intent(in) :: VX ! Voigt function K(x,y): x parameter
        real(kind=DP), intent(in) :: lineIntensity
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode

        voigtAsymptotic1Value = chiCorrectedLorentz(X, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode) * (1 + 1.5/VX**2)
    end function voigtAsymptotic1


    pure function voigtAsymptotic2(X, lorHWHM, VX, dopHWHM, lineIntensity) result(voigtAsymptotic2Value)
        ! Voigt ≈ Doppler + Lorentz + correction (y<<1, x>1)
        implicit none 

        real(kind=DP) :: voigtAsymptotic2Value
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: lorHWHM, dopHWHM
        real(kind=DP), intent(in) :: VX ! Voigt function K(x,y): x parameter
        real(kind=DP), intent(in) :: lineIntensity
        real(kind=DP), parameter :: U(9) = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
        real(kind=DP), parameter :: W(9) = [-0.688, 0.2667, 0.6338, 0.4405, 0.2529, 0.1601, 0.1131, 0.0853, 0.068]
        real(kind=DP) :: F
        integer :: I

        I = VX/0.5 - 1.00001
        F = 2. * (W(I)*(U(I+1)-VX) + W(I+1)*(VX-U(I)))
        voigtAsymptotic2Value = doppler(X, dopHWHM, lineIntensity) + (lorHWHM/(pi*X**2) * (1.+F)) * lineIntensity
    end function voigtAsymptotic2

! end of Content of Spectroscopy.f90

! content of ChiFactors.f90

    pure function noneChi(X, moleculeIntCode) result(noFactor)
        
        implicit none
        real :: noFactor
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !
        noFactor = 1.
    end function noneChi

    pure function tonkov(X, moleculeIntCode) result(tonkovFactor)
        ! MV Tonkov et al. "Measurements and empirical modeling of pure CO2 absorption in the 2.3-μm
        ! region at room temperature: far wings, allowed and collision-induced bands". In: Applied optics 35.24
        ! (1996), pp. 4863–4870.

        implicit none
        real :: tonkovFactor
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !

        tonkovFactor = 1. ! default value
        if (moleculeIntCode == 2) then ! CO2
            if (abs(X) > 3.) then
                if (abs(X) <= 150.) then
                    tonkovFactor = 1.084 * exp(-0.027*abs(X))
                else if (abs(x) <= 300.) then
                    tonkovFactor = 0.208 * exp(-0.016*abs(X))
                else
                    tonkovFactor = 0.025 * exp(-0.009*abs(X))
                end if
            end if
        end if
    end function tonkov

    pure function pollack(X, moleculeIntCode) result(pollackFactor)
        ! James B. Pollack et al. "Near-Infrared Light from Venus' Nightside: A Spectroscopic Analysis". In:
        ! Icarus 103 (1 1993), pp. 1–42. ISSN: 10902643. DOI: 10.1006/icar.1993.1055
        
        ! Pollack(1993):  recommended cutoff is 125 cm-1 or 160 cm-1 for the 1.2 um window complex (see p. 5.4 in the paper)
        
        implicit none
        real :: pollackFactor
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !

        pollackFactor = 1. ! default value of the chi-factor
        
        if (moleculeIntCode == 2) then
            if (abs(X) >= 3.) then
                if (abs(X) < 10.) then
                    pollackFactor = 1.35 * exp(-abs(X)/10.)
                else
                    pollackFactor = 0.614 * exp(-abs(X)/47.)
                    ! pollackFactor = 0.614 * exp(abs(X)/90.) ! for the 1.18 um window
                end if
            end if
        end if
    end function pollack

    pure function perrin(X, moleculeIntCode) result(perrinFactor)
        ! MY Perrin and JM Hartmann. "Temperature-dependent measurements and modeling of absorption
        ! by CO2-N2 mixtures in the far line-wings of the 4.3 μm CO2 band". In: Journal of Quantitative
        ! Spectroscopy and Radiative Transfer 42.4 (1989), pp. 311–317.
        
        implicit none
        real :: perrinFactor
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !
        real :: B1, B2, B3, S1, S2, S3

        perrinFactor = 1. ! default value for the chi-factor

        ! temperature dependence
        ! B(T) = alpha + beta * exp(-epsilon*T)
        B1 = 0.0888 - 0.16*exp(-0.00410*temperature)
        B2 = 0.0526 * exp(-0.00152*temperature)
        B3 = 0.0232

        ! boundaries
        S1 = 3.
        S2 = 30.
        S3 = 120.

        if (moleculeIntCode == 2) then 
            ! CO2
            if (abs(X) > S1) then
                if (abs(X) < S2) then
                    perrinFactor = exp(-B1 * (abs(X)-S1))
                else 
                    if (abs(X) < S3) then
                        perrinFactor = exp(-B1*(S2-S1) - B2*(abs(X)-S2))
                    else
                        perrinFactor = exp(-B1*(S2-S1) - B2*(S3-S2) - B3*(abs(X)-S3))
                    end if
                end if 
            end if
        end if
    end function perrin

! end of Content of ChiFactors.f90


    real function voigt(X, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode) ! new thread-safe version
        ! Implementation is based on the Humlicek method, improved by Kuntz:
        ! M. Kuntz. "A new implementation of the Humlicek algorithm for the calculation of the Voigt profile
        ! function". In: Journal of Quantitative Spectroscopy and Radiative Transfer 57.6 (1997), pp. 819–824.

        ! Additional features and deviations:
        ! - scheme is not recursive
        ! - add region 0 where pure Lorentz is used
        ! - the boundaries between regions 3 and 4 are deviated from Kuntz's initial study.
        ! - in region 4 asymptotical analytical approximations are used for increasing speed
        
        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: dopHWHM
        real(kind=DP), intent(in) :: lorHWHM
        real(kind=DP), intent(in) :: lineIntensity
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !

        real(kind=DP) :: VX, VY ! x and y parameters in the K(x,y) function
        
        real(kind=DP) :: VXsquared ! x**2
        real(kind=DP) :: Y1=0, Y2=0, Y3=0
        real(kind=DP) :: Y_2
        real(kind=DP) :: A1, B1, A2, B2, A3, B3, C3, D3, A4, B4, C4, D4, A5, B5, C5, D5, E5, &
                            A6, B6, C6, D6, E6
        
        VX = abs(sqln2 * X / dopHWHM)
        VY = sqln2 * lorHWHM / dopHWHM

        ! REGION 0: Lorentz domination: pure Lorentz with χ-corrected wing (if set)
        if (VX >= 15.) then
            voigt = chiCorrectedLorentz(X, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            return
        end if

        VXsquared = VX ** 2

        ! REGION 1: Voigt rational approximation 1
        if (VX + VY >= 15.0) then 
            
            Y1 = VY
            Y_2 = Y1 ** 2
            A1 = (0.2820948 + 0.5641896*Y_2) * Y1
            B1 = 0.5641896 * Y1
            A2 = 0.25 + Y_2 + Y_2**2
            B2 = Y_2 + Y_2 - 1.
            
            ! rational approximation for Voigt using A1, A2, B1, B2 coefficients
            voigt = (A1+B1*VXsquared) / (A2+B2*VXsquared+VXsquared**2) / sqrt(pi) / (dopHWHM/sqln2) * lineIntensity
        
        ! REGION 2: Voigt rational approximation 2       
        else if (VX + VY >= 5.5) then

            Y2 = VY
            Y_2 = Y2**2
            A3 = Y2 * (((0.56419*Y_2+3.10304)*Y_2+4.65456)*Y_2+1.05786)
            B3 = Y2 * ((1.69257*Y_2+0.56419)*Y_2+2.962)
            C3 = Y2 * (1.69257*Y_2-2.53885)
            D3 = Y2*0.56419
            A4 = (((Y_2+6.0)*Y_2+10.5)*Y_2+4.5)*Y_2+0.5625
            B4 = ((4.0*Y_2+6.0)*Y_2+9.0)*Y_2-4.5
            C4 = 10.5+6.0*(Y_2-1.0)*Y_2
            D4 = 4.0*Y_2-6.0

            voigt = (((D3*VXsquared+C3)*VXsquared+B3)*VXsquared+A3) / ((((VXsquared+D4)*VXsquared+C4)*VXsquared+B4)*VXsquared+A4) &
                            / sqrt(pi) / (dopHWHM/sqln2) * lineIntensity
        
        ! REGION 3: Voigt rational approximation 3
        else if (VX <= 1.0 .OR. VY >= 0.02) then 

            Y3 = VY
            A5 = ((((((((0.564224*Y3+7.55895)*Y3+49.5213)*Y3+204.510)*Y3+	&
                581.746)*Y3+1174.8)*Y3+1678.33)*Y3+1629.76)*Y3+973.778)*Y3+272.102
            B5 = ((((((2.25689*Y3+22.6778)*Y3+100.705)*Y3+247.198)*Y3+336.364)*	&
                Y3+220.843)*Y3-2.34403)*Y3-60.5644
            C5 = ((((3.38534*Y3+22.6798)*Y3+52.8454)*Y3+42.5683)*Y3+18.546)*Y3+	&
                4.58029
            D5 = ((2.25689*Y3+7.56186)*Y3+1.66203)*Y3-0.128922
            E5 = 0.971457E-3+0.564224*Y3
            A6 = (((((((((Y3+13.3988)*Y3+88.2674)*Y3+369.199)*Y3+1074.41)*Y3+	&
                2256.98)*Y3+3447.63)*Y3+3764.97)*Y3+2802.87)*Y3+1280.83)*Y3+	&
                272.102
            B6 = (((((((5.*Y3+53.5952)*Y3+266.299)*Y3+793.427)*Y3+1549.68)*Y3+	&
                2037.31)*Y3+1758.34)*Y3+902.306)*Y3+211.678
            C6 = (((((10.*Y3+80.3928)*Y3+269.292)*Y3+479.258)*Y3+497.302)*Y3+	&
                308.186)*Y3+78.866
            D6 = (((10.*Y3+53.5952)*Y3+92.7586)*Y3+55.0293)*Y3+22.0353
            E6 = (5.0*Y3+13.3988)*Y3+1.49645
            
            voigt = ((((E5*VXsquared+D5)*VXsquared+C5)*VXsquared+B5)*VXsquared+A5)/	&
                    (((((VXsquared+E6)*VXsquared+D6)*VXsquared+C6)*VXsquared+B6)*VXsquared+A6) / sqrt(pi) &
                            / (dopHWHM/sqln2) * lineIntensity
    
        else
            ! REGION 4: ASYMPTOTIC REGION: (1. < VX < 5.5 and y < 0.02)
            ! where instead of approximation of Voigt function, 
            ! can be used analytical asymptotical expressions based on series expansions 
            ! of Lorentz and Doppler shape functions	
            if (VX > 5.) then
                voigt = voigtAsymptotic1(X, lorHWHM, VX, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            else if (VX > sqrt(1.4)) then
                voigt = voigtAsymptotic2(X, lorHWHM, VX, dopHWHM, lineIntensity)
            else 
                voigt = doppler(X, dopHWHM, lineIntensity)
            end if

        end if
    end function voigt


    real function voigt_(X, dopHWHM, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
        ! Implementation is based on the Humlicek method, improved by Kuntz:
        ! M. Kuntz. "A new implementation of the Humlicek algorithm for the calculation of the Voigt profile
        ! function". In: Journal of Quantitative Spectroscopy and Radiative Transfer 57.6 (1997), pp. 819–824.

        ! Additional features and deviations:
        ! - scheme is not recursive
        ! - add region 0 where pure Lorentz is used
        ! - the boundaries between regions 3 and 4 are deviated from Kuntz's initial study.
        ! - in region 4 asymptotical analytical approximations are used for increasing speed
        
        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real(kind=DP), intent(in) :: dopHWHM
        real(kind=DP), intent(in) :: lorHWHM
        real(kind=DP), intent(in) :: lineIntensity
        procedure(chifactor), pointer :: chiFactorFuncPtr
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !

        real(kind=DP) :: VX, VY ! x and y parameters in the K(x,y) function
        
        real(kind=DP) :: VXsquared ! x**2
        real(kind=DP) :: Y1=0, Y2=0, Y3=0
        real(kind=DP) :: Y_2
        real(kind=DP) :: A1, B1, A2, B2, A3, B3, C3, D3, A4, B4, C4, D4, A5, B5, C5, D5, E5, &
                            A6, B6, C6, D6, E6
        
        ! TODO: figure out if this save can be removed (do when applying atmospheric level parallelization)
        save A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6, C3, C4, C5, C6, D3, D4, D5, D6, E5, E6
        
        VX = abs(sqln2 * X / dopHWHM)
        VY = sqln2 * lorHWHM / dopHWHM

        ! REGION 0: Lorentz domination: pure Lorentz with χ-corrected wing (if set)
        if (VX >= 15.) then
            voigt_ = chiCorrectedLorentz(X, lorHWHM, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            return
        end if

        VXsquared = VX ** 2

        ! REGION 1: Voigt rational approximation 1
        if (VX + VY >= 15.0) then 

            if (VY /= Y1) then
                Y1 = VY
                Y_2 = Y1 ** 2
                A1 = (0.2820948 + 0.5641896*Y_2) * Y1
                B1 = 0.5641896 * Y1
                A2 = 0.25 + Y_2 + Y_2**2
                B2 = Y_2 + Y_2 - 1.
            end if
            ! rational approximation for Voigt using A1, A2, B1, B2 coefficients
            voigt_ = (A1+B1*VXsquared) / (A2+B2*VXsquared+VXsquared**2) / sqrt(pi) / (dopHWHM/sqln2) * lineIntensity
        
        ! REGION 2: Voigt rational approximation 2       
        else if (VX + VY >= 5.5) then
            if (VY /= Y2) then
                Y2 = VY
                Y_2 = Y2**2
                A3 = Y2 * (((0.56419*Y_2+3.10304)*Y_2+4.65456)*Y_2+1.05786)
                B3 = Y2 * ((1.69257*Y_2+0.56419)*Y_2+2.962)
                C3 = Y2 * (1.69257*Y_2-2.53885)
                D3 = Y2*0.56419
                A4 = (((Y_2+6.0)*Y_2+10.5)*Y_2+4.5)*Y_2+0.5625
                B4 = ((4.0*Y_2+6.0)*Y_2+9.0)*Y_2-4.5
                C4 = 10.5+6.0*(Y_2-1.0)*Y_2
                D4 = 4.0*Y_2-6.0
            end if 

            voigt_ = (((D3*VXsquared+C3)*VXsquared+B3)*VXsquared+A3) / ((((VXsquared+D4)*VXsquared+C4)*VXsquared+B4)*VXsquared+A4) &
                            / sqrt(pi) / (dopHWHM/sqln2) * lineIntensity
        
        ! REGION 3: Voigt rational approximation 3
        else if (VX <= 1.0 .OR. VY >= 0.02) then 
            if (VY /= Y3) then
                Y3 = VY
                A5 = ((((((((0.564224*Y3+7.55895)*Y3+49.5213)*Y3+204.510)*Y3+	&
                    581.746)*Y3+1174.8)*Y3+1678.33)*Y3+1629.76)*Y3+973.778)*Y3+272.102
                B5 = ((((((2.25689*Y3+22.6778)*Y3+100.705)*Y3+247.198)*Y3+336.364)*	&
                    Y3+220.843)*Y3-2.34403)*Y3-60.5644
                C5 = ((((3.38534*Y3+22.6798)*Y3+52.8454)*Y3+42.5683)*Y3+18.546)*Y3+	&
                    4.58029
                D5 = ((2.25689*Y3+7.56186)*Y3+1.66203)*Y3-0.128922
                E5 = 0.971457E-3+0.564224*Y3
                A6 = (((((((((Y3+13.3988)*Y3+88.2674)*Y3+369.199)*Y3+1074.41)*Y3+	&
                    2256.98)*Y3+3447.63)*Y3+3764.97)*Y3+2802.87)*Y3+1280.83)*Y3+	&
                    272.102
                B6 = (((((((5.*Y3+53.5952)*Y3+266.299)*Y3+793.427)*Y3+1549.68)*Y3+	&
                    2037.31)*Y3+1758.34)*Y3+902.306)*Y3+211.678
                C6 = (((((10.*Y3+80.3928)*Y3+269.292)*Y3+479.258)*Y3+497.302)*Y3+	&
                    308.186)*Y3+78.866
                D6 = (((10.*Y3+53.5952)*Y3+92.7586)*Y3+55.0293)*Y3+22.0353
                E6 = (5.0*Y3+13.3988)*Y3+1.49645
            end if
            
            voigt_ = ((((E5*VXsquared+D5)*VXsquared+C5)*VXsquared+B5)*VXsquared+A5)/	&
                    (((((VXsquared+E6)*VXsquared+D6)*VXsquared+C6)*VXsquared+B6)*VXsquared+A6) / sqrt(pi) &
                            / (dopHWHM/sqln2) * lineIntensity
    
        else
            ! REGION 4: ASYMPTOTIC REGION: (1. < VX < 5.5 and y < 0.02)
            ! where instead of approximation of Voigt function, 
            ! can be used analytical asymptotical expressions based on series expansions 
            ! of Lorentz and Doppler shape functions	
            if (VX > 5.) then
                voigt_ = voigtAsymptotic1(X, lorHWHM, VX, lineIntensity, chiFactorFuncPtr, moleculeIntCode)
            else if (VX > sqrt(1.4)) then
                voigt_ = voigtAsymptotic2(X, lorHWHM, VX, dopHWHM, lineIntensity)
            else 
                voigt_ = doppler(X, dopHWHM, lineIntensity)
            end if

        end if
    end function voigt_


    subroutine interpolation_step(coarserP, coarser, coarserL, finerP, finer, finerL)
        real, intent(in)    :: coarserP(:), coarser(:), coarserL(:)
        real, intent(inout) :: finerP(:), finer(:), finerL(:)
        integer :: J, I, M, N
        real(kind=DP) :: A

        N = size(coarserP)
        do J = 1, N
            I = 2 * J - 1
            finerP(I) = finerP(I) + coarserP(J)
            A = 0.375 * coarserP(J) + 0.75 * coarser(J) - 0.125 * coarserL(J)
            if (A < coarserP(J) .or. A < coarser(J)) A = 0.5 * (coarserP(J) + coarser(J))
            finer(I) = finer(I) + A
            finerL(I) = finerL(I) + coarser(J)
            
            M = I + 1
            finerP(M) = finerP(M) + coarser(J)
            A = 0.375 * coarserL(J) + 0.75 * coarser(J) - 0.125 * coarserP(J)
            if (A < coarserL(J) .or. A < coarser(J)) A = 0.5 * (coarserL(J) + coarser(J))
            finer(M) = finer(M) + A
            finerL(M) = finerL(M) + coarserL(J)
        end do
    end subroutine interpolation_step

    
    subroutine cascadeInterpolation(NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9, &
        RK, RK0, RK0L, RK0P, RK1, RK1L, RK1P, RK2, RK2L, RK2P, RK3, RK3L, RK3P, RK4, RK4L, RK4P, &
        RK5, RK5L, RK5P, RK6, RK6L, RK6P, RK7, RK7L, RK7P, RK8, RK8L, RK8P, RK9, RK9L, RK9P)

        implicit none

        integer, intent(in) :: NT, NT0, NT1, NT2, NT3, NT4, NT5, NT6, NT7, NT8, NT9
        real, intent(inout) :: RK(:), RK0(:), RK1(:), RK2(:), RK3(:), RK4(:), RK5(:), &
                               RK6(:), RK7(:), RK8(:), RK9(:)
        real, intent(inout) :: RK0L(:), RK1L(:), RK2L(:), RK3L(:), RK4L(:), RK5L(:), &
                               RK6L(:), RK7L(:), RK8L(:), RK9L(:)
        real, intent(inout) :: RK0P(:), RK1P(:), RK2P(:), RK3P(:), RK4P(:), RK5P(:), &
                               RK6P(:), RK7P(:), RK8P(:), RK9P(:)

        integer :: I, J
        real :: A

        ! Interpolate through each grid level
        call interpolation_step(RK0P, RK0, RK0L, RK1P, RK1, RK1L)
        call interpolation_step(RK1P, RK1, RK1L, RK2P, RK2, RK2L)
        call interpolation_step(RK2P, RK2, RK2L, RK3P, RK3, RK3L)
        call interpolation_step(RK3P, RK3, RK3L, RK4P, RK4, RK4L)
        call interpolation_step(RK4P, RK4, RK4L, RK5P, RK5, RK5L)
        call interpolation_step(RK5P, RK5, RK5L, RK6P, RK6, RK6L)
        call interpolation_step(RK6P, RK6, RK6L, RK7P, RK7, RK7L)
        call interpolation_step(RK7P, RK7, RK7L, RK8P, RK8, RK8L)
        call interpolation_step(RK8P, RK8, RK8L, RK9P, RK9, RK9L)

        ! Final interpolation step to the finest grid (RK)
        I = 1
        do J = 1, NT9
            I = I + 1
            A = RK9P(J)*0.375 + RK9(J)*0.75 - RK9L(J)*0.125
            if (A < RK9P(J) .or. A < RK9(J)) A = 0.5 * (RK9P(J) + RK9(J))
            RK(I) = RK(I) + A

            I = I + 1
            RK(I) = RK(I) + RK9(J)

            I = I + 1
            A = RK9L(J)*0.375 + RK9(J)*0.75 - RK9P(J)*0.125
            if (A < RK9L(J) .or. A < RK9(J)) A = 0.5 * (RK9L(J) + RK9(J))
            RK(I) = RK(I) + A

            I = I + 1
            RK(I) = RK(I) + RK9L(J)
        end do
    end subroutine cascadeInterpolation
    
    
    subroutine outputSubintervalResults(subStartWV, subEndWV, NT, RK, H)
        ! Output results for current subinterval
        implicit none
        
        real(kind=DP), intent(in) :: subStartWV, subEndWV
        integer, intent(in) :: NT
        real, intent(in) :: RK(:)
        real(kind=DP), intent(in) :: H
        real(kind=DP) :: wv
        integer :: i
        
        ! Output data points for this subinterval
        do i = 1, NT
            wv = subStartWV + (i-1) * H
            if (wv > subEndWV) exit
            
            write(SPECTRA_FCHAN,'(F12.6, E15.6)') wv, RK(i)
        end do
        
    end subroutine outputSubintervalResults
    
    subroutine outputResults()
        ! This subroutine is no longer used
    end subroutine outputResults
    
    
    function is_iostat_end(ios) result(is_end)
        ! Check if iostat indicates end of file
        integer, intent(in) :: ios
        logical :: is_end
        
        is_end = (ios < 0)
    end function is_iostat_end

end program MarfaSimple 
