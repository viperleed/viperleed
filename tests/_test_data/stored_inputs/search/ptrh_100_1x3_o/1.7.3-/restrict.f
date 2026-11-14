
C  This subroutine serves to restrict individual parameters to certain values
C  inside the search. These values are stored in the PARIND(IPARAM,IPOP) array.
C  Parameters are counted as listed in control file search.steu, including
C  a concentration parameter after each atomic site, and the domain weight
C  parameters for each domain at the end of the PARIND array.
C
C  perform restrictions in the following fashion: e.g.
C
C     PARIND(1,IPOP) = 5
C
C  or
C
C     PARIND(5,IPOP) = PARIND(1,IPOP)
C
C  etc.

      subroutine restrict(NPRMK,NPS,PARIND,IPOP)

      INTEGER NPRMK,NPS
      INTEGER PARIND
      DIMENSION PARIND(NPRMK,NPS)

C  begin restrictions

      PARIND(15,IPOP) = 1

C  end restrictions

      RETURN

      END
