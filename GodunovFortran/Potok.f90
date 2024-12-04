subroutine POTOK(AK,P1,R1,U1,P2,R2,U2,PBIG,UBIG,RBIG1,RBIG2,DL1,DL2,DP1,DP2,UDOT,Pf,Uf,Rf) bind(c, name = 'potok') 
    use iso_c_binding  
    IMPLICIT NONE   
!
! INPUT 
!
      REAL(c_double) AK,P1,R1,U1,P2,R2,U2,UDOT
      REAL(c_double) PBIG,UBIG                 ! parameters on contact surface
      REAL(c_double) RBIG1,RBIG2               ! density right and left of the surface
      REAL(c_double) DL1,DL2,DP1,DP2           ! speed of characteristics
!      
! OUTPUT 
!
      REAL(c_double) Pf,Uf,Rf
!
! LOCALS
!
      REAL(c_double) AK1,AK2,A1,A2
      REAL(c_double) CZJM,CZJP,CZT

      AK1=AK
      AK2=AK
      Uf=UBIG
      Pf=PBIG
      if (UBIG-UDOT > 0. ) then    

        if (DL1-DL2 .eq. 0.) then  
          if (DL1-UDOT >= 0.) then 
            Rf=R1
            Uf=U1
            Pf=P1
          else                     
            Rf=RBIG1
          endif
        else                       
          if (DL2-UDOT <= 0.) then
            Rf=RBIG1             
          else
            if (DL1-UDOT >= 0.) then  
              Rf=R1              
              Uf=U1
              Pf=P1
            else 
              CZJM=SQRT(AK1*P1/R1)
              CZT=((AK1-1.)*(U1-UDOT)+2.*CZJM)/(AK1+1.)
              Uf=CZT+UDOT
              Pf=P1*(CZT/CZJM)**(2.*AK1/(AK1-1.))
              Rf=AK1*Pf/(CZT*CZT)
            endif      
          endif      
        endif      

      else                        

        if (DP1-DP2 .eq. 0.) then  
          if (DP1-UDOT <= 0.) then
            Rf=R2
            Uf=U2
            Pf=P2
          else                  
            Rf=RBIG2
          endif
        else                       
          if (DP1-UDOT >=0.) then
            Rf=RBIG2             
          else
            if (DP2-UDOT <= 0.) then  
              Rf=R2             
              Uf=U2
              Pf=P2
            else                 
              CZJP=SQRT(AK2*P2/R2)
              CZT=(-(AK2-1.)*(U2-UDOT)+2.*CZJP)/(AK2+1.)
              Uf=-CZT+UDOT
              Pf=P2*(CZT/CZJP)**(2.*AK2/(AK2-1.))
              Rf=AK2*Pf/(CZT*CZT)
            endif      
          endif      
        endif      
      
      endif
      
      return
      end