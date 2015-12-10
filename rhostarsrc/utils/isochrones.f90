module isochrones
    implicit none
    integer, parameter :: niso =9
    character(80), dimension(niso) :: isofile
    data isofile /"isochrones/UBVRIJHKsKp/fehm25afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehm20afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehm15afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehm10afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehm05afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehp00afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehp02afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehp03afep0.UBVRIJHKsKp", &
                  "isochrones/UBVRIJHKsKp/fehp05afep0.UBVRIJHKsKp"/
end module isochrones
