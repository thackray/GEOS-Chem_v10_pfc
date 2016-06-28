!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_PFCA_mod.F90
!
! !DESCRIPTION: Defines the HEMCO extension for the GEOS-Chem PFCA
!    specialty simulation.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_PFCA_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State   ! Derived type for HEMCO state
  USE HCOX_State_Mod, ONLY : Ext_State   ! Derived type for External state

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HcoX_PFCA_Run
  PUBLIC  :: HcoX_PFCA_Init
  PUBLIC  :: HcoX_PFCA_Final
!
! !REMARKS:
!
!  References:
!  ============================================================================
!  (1 ) Zhang, Y., and S. Tao, Global atmospheric emission inventory of
!        polycyclic aromatic hydrocarbons (PAHs) for 2004. Atm Env, 43, 812-819,
!        2009.
!  (2 ) Friedman, C.L, and N.E. Selin, Long-Range Atmospheric Transport of
!        Polycyclic Aromatic Hydrocarbons: A Global 3-D Model Analysis
!        Including Evaluation of Arctic Sources, Environ. Sci. Technol., 46(17),
!        9501-9510, 2012.
!  (3 ) Friedman, C.L., Y. Zhang, and N.E. Selin, Climate change and
!        emissions impacts on atmospheric PAH transport to the Arctic, Environ.
!        Sci. Technol., 48, 429-437, 2014.
!  (4 ) Friedman, C.L., J.R. Pierce, and N.E. Selin, Assessing the influence of
!        secondary organic versus primary carbonaceous aerosols on long-range
!        atmospheric polycyclic aromatic hydrocarbon transport, Environ. Sci.
!        Technol., 48(6), 3293-3302, 2014.
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin    - Initial Version
!  04 Jan 2011 - C.L. Friedman - Expansion on initial version
!  19 Aug 2014 - M. Sulprizio  - Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  ! Fields required by module
  INTEGER                       :: ExtNr     = 115   ! HEMCO Extension number
  INTEGER                       :: IDTFTI = 80      ! Index # for FTI tracer
  INTEGER                       :: IDTFTOH   = 79   ! Index # for FTOH tracer

  REAL*8,  PARAMETER            :: SMALLNUM = 1D-20

  ! Pointers to emission arrays read from disk
  REAL(sp), POINTER             :: FTOH_TOT_EM(:,:) => NULL()
  REAL(sp), POINTER             :: FTI_TOT_EM(:,:) => NULL()

  ! Calculated emissions of FTOH and FTI
  REAL(hp), ALLOCATABLE, TARGET :: EFTOH(:,:,:)
  REAL(hp), ALLOCATABLE, TARGET :: EFTI(:,:,:)

  ! For diagnostics
  REAL(hp), ALLOCATABLE, TARGET :: SUM_G_EM  (:,:)
  REAL(hp), ALLOCATABLE, TARGET :: SUM_OF_ALL(:,:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_PFCA_run 
!
! !DESCRIPTION: Subroutine HcoX\_PFCA\_Run computes emissions of FTOH, FTI
!  for the GEOS-Chem PFCA specialty simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_PFCA_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    ! HEMCO modules
    USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
    USE HCO_FluxArr_Mod,   ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
    TYPE(Ext_State),  POINTER       :: ExtState    ! Options for PFCA sim
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version based on EMISSMERCURY
!  29 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  14 Apr 2014 - R. Yantosca - Prevent div-by-zero error w/ SUM_OF_ALL
!  19 Aug 2014 - M. Sulprizio- Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
    REAL*8, PARAMETER :: R          = 8.31d0  


!
! !LOCAL VARIABLES:
!
    INTEGER           :: I, J, L
    INTEGER           :: PBL_MAX
    INTEGER           :: MONTH,            YEAR
    REAL*8            :: F_OF_PBL,         TK
    REAL*8            :: T_FTOH, T_FTI
    REAL*8            :: AIR_VOL
    LOGICAL, SAVE     :: FIRST = .TRUE.
    LOGICAL           :: aIR

    ! Pointers for diagnostics
    REAL(hp), POINTER :: Arr3D(:,:,:) => NULL()

    !=======================================================================
    ! HCOX_PFCA_RUN begins here!
    !=======================================================================

    ! Return if extension not turned on
    IF ( .NOT. ExtState%PFCA ) RETURN

    ! Enter
    CALL HCO_ENTER( 'HCOX_PFCA_Run (hcox_PFCA_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! am I root? 
    aIR = am_I_Root

    !=======================================================================
    ! Get pointers to gridded data imported through config. file
    !=======================================================================
    IF ( FIRST ) THEN

       CALL HCO_GetPtr( aIR, 'PFCA_FTOH', FTOH_TOT_EM, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( aIR, 'PFCA_FTI', FTI_TOT_EM, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FIRST = .FALSE.

    ENDIF

    ! Maximum extent of the PBL [model level]
    IF ( .NOT. ASSOCIATED(ExtState%PBL_MAX) ) THEN
       CALL HCO_ERROR ( 'PBL_MAX not defined in ExtState!', RC )
       RETURN
    ELSE
       PBL_MAX = DBLE( ExtState%PBL_MAX )
    ENDIF

    ! Loop over grid boxes
    DO J = 1, HcoState%Ny
    DO I = 1, HcoState%Nx

       F_OF_PBL = 0d0 
       T_FTOH = 0d0       
       T_FTI = 0d0       

       ! Here, save the total from the emissions array
       ! into the T_FTX variable [kg/m2/s]
       T_FTOH = FTOH_TOT_EM(I,J)
       T_FTI = FTI_TOT_EM(I,J)

       !====================================================================
       !          
       ! Then partition emis throughout PBL; store into STT [kg]
       ! Now make sure STT does not underflow (cdh, bmy, 4/6/06; eck 9/20/10)
       !====================================================================

       ! Loop up to max PBL level
       DO L = 1, PBL_MAX

          ! Get temp [K]
          TK = ExtState%TK%Arr%Val(I,J,L)


          ! Get air volume (m^3)
          AIR_VOL   = ExtState%AIRVOL%Arr%Val(I,J,L) 
            
          ! Fraction of PBL that box (I,J,L) makes up [unitless]
          F_OF_PBL = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) 

          ! FTOH
          EFTOH(I,J,L) = F_OF_PBL * T_FTOH

          ! FTI
          EFTI(I,J,L)  = F_OF_PBL * T_FTI

       ENDDO

       !==================================================================
       ! 
       ! through bottom layer to top of PBL for storage in ND53 diagnostic
       !==================================================================
       SUM_G_EM(I,J)  = SUM(EFTOH(I,J,1:PBL_MAX)) + SUM(EFTI(I,J,1:PBL_MAX))  
       
       SUM_OF_ALL(I,J) = SUM_G_EM(I,J)

    ENDDO
    ENDDO

    !=======================================================================
    ! Add FTX emissions to HEMCO data structure & diagnostics
    !=======================================================================

  
    !----------------------
    ! GASEOUS EMISSIONS
    !----------------------
    IF ( IDTFTOH > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => EFTOH(:,:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr3D, IDTFTOH, RC, ExtNr=ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: EFTOH', RC )
          RETURN 
       ENDIF

    ENDIF

    IF ( IDTFTI > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => EFTI(:,:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr3D, IDTFTI, RC, ExtNr=ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: EFTI', RC )
          RETURN 
       ENDIF

    ENDIF

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_PFCA_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_PFCA_Init 
!
! !DESCRIPTION: Subroutine HcoX\_PFCA\_Init initializes the HEMCO
!  PFCA extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_PFCA_Init( am_I_Root, HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,   ONLY : GetExtNr
    USE HCO_STATE_MOD,     ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! Hemco state 
    INTEGER,          INTENT(INOUT)  :: RC 

! !REVISION HISTORY:
!  19 Aug 2014 - M. Sulprizio- Initial version
!  01 May 2015 - R. Yantosca - Bug fix: need to zero arrays after allocating
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, nSpc
    CHARACTER(LEN=255)             :: MSG 

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=======================================================================
    ! HCOX_PFCA_INIT begins here!
    !=======================================================================

    ! Get the extension number
    ExtNr = GetExtNr( TRIM( ExtName ) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter HEMCO
    CALL HCO_ENTER( 'HcoX_PFCA_Init (hcox_PFCA_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set species IDs      
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use PFCA emissions module (extension module)'
       CALL HCO_MSG( MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(MSG)
       ENDDO
    ENDIF

    ! Set up tracer indices
    DO N = 1, nSpc
       SELECT CASE( TRIM( SpcNames(N) ) )
          CASE( 'FTOH' )
             IDTFTOH   = HcoIDs(N)
          CASE( 'FTI' )
             IDTFTI = HcoIDs(N)
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! ERROR: FTOH tracer is not found!
    IF ( IDTFTOH <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find FTOH tracer in list of species!', RC )
       RETURN
    ENDIF
    
    ! ERROR! FTI tracer is not found
    IF ( IDTFTI <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find FTI tracer in list of species!', RC )
       RETURN
    ENDIF

    ! Activate met fields required by this extension
    ExtState%AIRVOL%DoUse      = .TRUE. 
    ExtState%FRAC_OF_PBL%DoUse = .TRUE. 
    ExtState%TK%DoUse          = .TRUE. 

    ! Activate this extension
    ExtState%PFCA           = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    ALLOCATE( EFTOH ( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EFTOH', RC )
       RETURN
    ENDIF 
    EFTOH = 0.0e0_hp

    ALLOCATE( EFTI( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EFTI', RC )
       RETURN
    ENDIF 
    EFTI = 0.0e0_hp

    ALLOCATE( SUM_G_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate SUM_G_EM', RC )
       RETURN
    ENDIF     
    SUM_G_EM = 0.0e0_hp

    ALLOCATE( SUM_OF_ALL( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate SUM_OF_ALL', RC )
       RETURN
    ENDIF 
    SUM_OF_ALL = 0.0e0_hp

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HCOX_PFCA_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_PFCA_Final
!
! !DESCRIPTION: Subroutine HcoX\_PFCA\_Final finalizes the HEMCO
!  extension for the GEOS-Chem PFCA specialty simulation.  All module
!  arrays will be deallocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_PFCA_Final()
!
! !REVISION HISTORY:
!  19 Aug 2014 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_PFCA_FINAL begins here!
    !=======================================================================
    IF ( ALLOCATED( EFTOH      ) ) DEALLOCATE( EFTOH      )
    IF ( ALLOCATED( EFTI       ) ) DEALLOCATE( EFTI       )
    IF ( ALLOCATED( SUM_G_EM   ) ) DEALLOCATE( SUM_G_EM   )
    IF ( ALLOCATED( SUM_OF_ALL ) ) DEALLOCATE( SUM_OF_ALL )

  END SUBROUTINE HCOX_PFCA_Final
!EOC
END MODULE HCOX_PFCA_Mod
