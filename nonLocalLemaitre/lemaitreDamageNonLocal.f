! User material subroutine for plasticity combined with non-local Lemaitre damage, compatible with Abaqus/2020 or newer 
! The code is distributed under a BSD license     
      
      
! Andrew Whelan (andrew.whelan@ucdconnect.ie)
! University College Dublin
      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc'    

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
      
      dimension eelas(ntens),eplas(ntens),stress0(ntens),flow(ntens), olds(ntens),oldpl(ntens)
      
      parameter(toler=1.d-6,newton=20)

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      Sy=props(3) ! Yield stress multiplier
      xn=props(4) ! Strain hardening exponent
      sDam=props(5) ! Lemaitre damage denominator
      bDam=props(6) ! Lemaitre damage exponent
      charLength=props(7) ! Characterisitic length
      damage=statev(2+2*ntens) ! Local damage
      damageNL=temp+dtemp ! Non-local damage
      damageNLOld=temp+dtemp ! Old value of the non-local damage
      damageOld=damage ! Old value of the local damage
      
      gd=1.d0-damageNL ! Degredation due to damage

      do k1=1,ntens
            stress0(k1)=statev(k1+2+2*ntens) !Recover undegraded stress    
      enddo  
      
      call rotsig(statev(1),drot,eelas,2,ndi,nshr)
      call rotsig(statev(ntens+1),drot,eplas,2,ndi,nshr)
      eqplas=statev(1+2*ntens)
      olds=stress
      oldpl=eplas

!     Build elastic stiffness matrix
      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      ek=E/(3.d0*(1.d0-2.d0*xnu))
      
      do i=1,3
       do j=1,3
        ddsdde(j,i)=elam
       end do
       ddsdde(i,i)=2.d0*eg+elam
      end do
      do i=4,ntens
       ddsdde(i,i)=eg
      end do

!     Calculate predictor stress and elastic strain
      stress0=stress0+matmul(ddsdde,dstran)
      eelas=eelas+dstran

!     Calculate equivalent von Mises stress
      Smises=(stress0(1)-stress0(2))**2+(stress0(2)-stress0(3))**2
     1 +(stress0(3)-stress0(1))**2
      do i=4,ntens
       Smises=Smises+6.d0*stress0(i)**2
      end do
      Smises=sqrt(Smises/2.d0) 
      
!     Get yield stress from the specified hardening curve
      Sf=Sy*(0.0001+eqplas)**xn
      
!     Determine if active yielding
      if (Smises.gt.(1.d0+toler)*Sf) then

!     Calculate the flow direction
       Sh=(stress0(1)+stress0(2)+stress0(3))/3.d0
       flow(1:3)=(stress0(1:3)-Sh)/Smises
       flow(4:ntens)=stress0(4:ntens)/Smises
       
       deqpl=0.d0
       Et=xn*Sy*(0.0001+eqplas)**(xn-1.d0)

!     Newton loop to calculate plastic multiplier
       do kewton=1,newton
        rhs=Smises-(3.d0*eg*deqpl)/gd-Sf
        deqpl=deqpl+rhs/((3.d0*eg/gd)+Et)
        Sf=Sy*(0.0001+eqplas+deqpl)**xn
        Et=xn*Sy*(0.0001+eqplas+deqpl)**(xn-1.d0)

        if(abs(rhs).lt.toler*Sy) exit
       end do
       if (kewton.eq.newton) write(7,*)'WARNING: plasticity loop failed'

       if (deqpl.lt.0.d0) then
            deqpl=0.d0
       endif

!     Energy release rate associated with the damage
       Y=-((Sf**2)/(6*eg)+(Sh**2)/(2*ek))      
       YStar=((-Y/sDam)**bDam) 

!     Update local damage
       damage=damageOld+YStar*deqpl/gd

       if (damage.GT.1.d0) then
            damage=0.99
       endif
     

! update stresses and strains
       stress0(1:3)=flow(1:3)*Sf+Sh
       eplas(1:3)=eplas(1:3)+3.d0/2.d0*flow(1:3)*deqpl/gd
       eelas(1:3)=eelas(1:3)-3.d0/2.d0*flow(1:3)*deqpl/gd
       stress0(4:ntens)=flow(4:ntens)*Sf
       eplas(4:ntens)=eplas(4:ntens)+3.d0*flow(4:ntens)*deqpl/gd
       eelas(4:ntens)=eelas(4:ntens)-3.d0*flow(4:ntens)*deqpl/gd
       eqplas=eqplas+deqpl

       stress=stress0


!    Calculate the plastic strain energy density
       do i=1,ntens
        spd=spd+deqpl*Sf
       end do

       sse = 0.d0
       do k1 = 1, ntens
          sse=sse+stress0(k1)*eelas(k1)/2.d0  !Update specific elastic strain energy
       enddo 
      
!     Formulate the jacobian (material tangent)   
       effg=eg*Sf/Smises
       efflam=(E/(1.d0-2.d0*xnu)-2.d0*effg)/3.d0
       effhrd=3.d0*eg*Et/(3.d0*eg+Et)-3.d0*effg
       do i=1,3
        do j=1,3
         ddsdde(j,i)=efflam
        end do
        ddsdde(i,i)=2.d0*effg+efflam
       end do
       do i=4,ntens
        ddsdde(i,i)=effg
       end do

       do i=1,ntens
        do j=1,ntens
         ddsdde(j,i)=ddsdde(j,i)+effhrd*flow(j)*flow(i)
        end do
       end do
      endif

!     Degrade stress and material tangent
      stress=gd*stress0
      ddsdde=gd*ddsdde

!    Update state variables
      statev(1:ntens)=eelas
      statev((ntens+1):2*ntens)=eplas
      statev(1+2*ntens)=eqplas
      statev(2+2*ntens)=damage

      statev(3+2*ntens:2+2*ntens+ntens)=stress0
      statev(3+2*ntens+ntens)=damageNLOld

!    Variables for calculating of non-local damage
      rpl=(damage-damageNL)/charLength**2 
      drpldt=-1.d0/charLength**2

      return
      end
