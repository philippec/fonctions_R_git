C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, prepare the seq. of aggloms. and "horiz."    C
C  order of objects for plotting the dendrogram using S routine C
C  'plclust'.                                                   C
C                                                               C
C  Parameters:                                                  C
C                                                               C
C  IA, IB:       vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  IIA, IIB:     used to store IA and IB values differently     C
C                (in form needed for S command `plclust`        C
C  IORDER:       "horiz." order of objects for dendrogram       C
C                                                               C
C  F. Murtagh, ESA/ESO/STECF, Garching, June 1991               C
C                                                               C
C  HISTORY                                                      C
C                                                               C
C  Adapted from routine HCASS, which additionally determines    C
C   cluster assignments at all levels, at extra comput. expense C
C                                                               C
C---------------------------------------------------------------C
      SUBROUTINE HCASS2(N,IA,IB,IORDER,IIA,IIB)
c Args
      INTEGER N,IA(N),IB(N),IORDER(N),IIA(N),IIB(N)
c Var
      INTEGER I, J, K, K1, K2, LOC
C
C     Following bit is to get seq. of merges into format acceptable to plclust
C     I coded clusters as lowest seq. no. of constituents; S's `hclust' codes
C     singletons as -ve numbers, and non-singletons with their seq. nos.
C
      DO 912 I=1,N
         IIA(I)=IA(I)
         IIB(I)=IB(I)
  912 CONTINUE
      DO 915 I=1,N-2
C        In the following, smallest (+ve or -ve) seq. no. wanted
         K=MIN(IA(I),IB(I))
         DO 913 J=I+1, N-1
            IF(IA(J).EQ.K) IIA(J)=-I
            IF(IB(J).EQ.K) IIB(J)=-I
  913    CONTINUE
  915 CONTINUE
      DO 916 I=1,N-1
         IIA(I)=-IIA(I)
         IIB(I)=-IIB(I)
  916 CONTINUE
      DO 917 I=1,N-1
         IF (IIA(I).GT.0.AND.IIB(I).LT.0) THEN
            K = IIA(I)
            IIA(I) = IIB(I)
            IIB(I) = K
         ENDIF
         IF (IIA(I).GT.0.AND.IIB(I).GT.0) THEN
            K1 = MIN(IIA(I),IIB(I))
            K2 = MAX(IIA(I),IIB(I))
            IIA(I) = K1
            IIB(I) = K2
         ENDIF
  917 CONTINUE
C
C
C     NEW PART FOR `ORDER'
C
      IORDER(1) =IIA(N-1)
      IORDER(2) =IIB(N-1)
      LOC=2
      DO 175 I=N-2,1,-1
        DO 169 J=1,LOC
          IF(IORDER(J).EQ.I) THEN
C           REPLACE IORDER(J) WITH IIA(I) AND IIB(I)
            IORDER(J)=IIA(I)
            IF (J.EQ.LOC) THEN
                LOC=LOC+1
                IORDER(LOC)=IIB(I)
                GOTO 171
            ENDIF
            LOC=LOC+1
            DO 95 K=LOC,J+2,-1
               IORDER(K)=IORDER(K-1)
  95        CONTINUE
            IORDER(J+1)=IIB(I)
            GOTO 171
          ENDIF
 169    CONTINUE
C       SHOULD NEVER REACH HERE
 171    CONTINUE
 175  CONTINUE
C
C
      DO 181 I=1,N
         IORDER(I) = -IORDER(I)
 181  CONTINUE
C
C
      RETURN
      END
