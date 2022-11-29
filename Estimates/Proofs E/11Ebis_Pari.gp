
  /* BORDER
     y = 4/5: (1+x)*(1+x*(4/5))*(1-2*g*x*(4/5)*z+(x*(4/5)*z)^2)*(1-(4/5))*(1+2*g*(4/5)*z+((4/5)*z)^2)
  */

  L= (1+x)*(1+x*(4/5))*(1-2*g*x*(4/5)*z+(x*(4/5)*z)^2)*(1-(4/5))*(1+2*g*(4/5)*z+((4/5)*z)^2)
  Tx=factor(deriv(L,x))
  Tz=factor(deriv(L,z))
  Tg=factor(deriv(L,g))
  Lx=Tx[2,1]
  Lz=Tz[3,1]
  Lg=Tg[4,1]

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Lz,Lg,g))
  /* All factors in T2 are strictly positive */

  /* BORDER */
  /* x = 0: already discussed in greater generality (it is <= 32/27) */  
  ....
  /* z = 0: (1+x)*(1+x*(4/5))*(1-(4/5)) <= 2*(1+(4/5))*(1-(4/5)) <= 18/25 */
  ....
  /* g = 0: (1+x)*(1+x*(4/5))*(1+(x*(4/5)*z)^2)*(1-(4/5))*(1+((4/5)*z)^2) increasing in x and z, maximized for x = z = 1
            2*(1+(4/5))*(1+(4/5)^2)*(1-(4/5))*(1+(4/5)^2) <= 1.94  */  
  ....
  /* x = 1: 2*(1+(4/5))*(1-2*g*(4/5)*z+(1*(4/5)*z)^2)*(1-(4/5))*(1+2*g*(4/5)*z+((4/5)*z)^2) maximized for g = 0
            2*(1+(4/5))*(1+(1*(4/5)*z)^2)*(1-(4/5))*(1+((4/5)*z)^2) maximized for z = 1 
            2*(1+(4/5))*(1+(4/5)^2)*(1-(4/5))*(1+(4/5)^2) <= 1.94 */  
  ....
  /* z = 1: (1+x)*(1+x*(4/5))*(1-2*g*x*(4/5)+(x*(4/5))^2)*(1-(4/5))*(1+2*g*(4/5)+((4/5))^2) */
    L=(1+x)*(1+x*(4/5))*(1-2*g*x*(4/5)+(x*(4/5))^2)*(1-(4/5))*(1+2*g*(4/5)+((4/5))^2)
    Tx=factor(deriv(L,x))
    Tg=factor(deriv(L,g))
    Lx=Tx[2,1]
    Lg=Tg[3,1]

    T1=factor(polresultant(Lx,Lg,g))
    res1y=T1[2,1]
    
    u=polrootsreal(res1y,[0,1]) /* length 1 */
    x0=u[1]
    v=polrootsreal(subst(Lg,x,x0),[0,1]) /* length 1 */
    g0=v[1]
    subst(subst(L,x,x0),g,g0)
    \\ 0.63007726327792155461307328966221487569  <----   
    
    /* x = 0: already discussed in greater generality (it is <= 32/27) */
    /* g = 0: already discussed in greater generality (it is <= 1.94) */
    /* x = 1: already discussed in greater generality (it is <= 1.94) */
    /* g = 1: (1+x)*(1+x*(4/5))*(1-2*g*x*(4/5)+(x*(4/5))^2)*(1-(4/5))*(1+2*g*(4/5)+((4/5))^2) <= 0.7 */ /* PARI my(g=1);ploth(x=0,1,(1+x)*(1+x*(4/5))*(1-2*g*x*(4/5)+(x*(4/5))^2)*(1-(4/5))*(1+2*g*(4/5)+((4/5))^2)) */
  ....
  /* g = 1: (1+x)*(1+x*(4/5))*(1-2*x*(4/5)*z+(x*(4/5)*z)^2)*(1-(4/5))*(1+2*(4/5)*z+((4/5)*z)^2) */  
    L=(1+x)*(1+x*(4/5))*(1-2*x*(4/5)*z+(x*(4/5)*z)^2)*(1-(4/5))*(1+2*(4/5)*z+((4/5)*z)^2)
    Tx=factor(deriv(L,x))
    Tz=factor(deriv(L,z))    
    Lx=Tx[3,1]
    Lz=Tz[5,1]

    T1=factor(polresultant(Lx,Lz,z))
    res1z=T1[2,1]
    
    u=polrootsreal(res1z,[0,1]) /* length 1 */
    x0=u[1]
    v=polrootsreal(subst(Lz,x,x0),[0,1]) /* length 1 */
    z0=v[1]
    subst(subst(L,x,x0),z,z0)
    \\ 0.53110634391677907976262347514974616005  <----
    
    /* x = 0: already discussed in greater generality (it is <= 32/27) */
    /* z = 0: already discussed in greater generality (it is <= 18/25) */
    /* x = 1: already discussed in greater generality (it is <= 1.94) */
    /* z = 1: (1+x)*(1+x*(4/5))*(1-2*x*(4/5)*z+(x*(4/5)*z)^2)*(1-(4/5))*(1+2*(4/5)*z+((4/5)*z)^2) <= 0.7 */ /* PARI my(z=1);ploth(x=0,1,(1+x)*(1+x*(4/5))*(1-2*x*(4/5)*z+(x*(4/5)*z)^2)*(1-(4/5))*(1+2*(4/5)*z+((4/5)*z)^2)) */
.........................................................................

  /* BORDER
     z = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+2*g*y+y^2) 
  */

  L= (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+2*g*y+y^2)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tg=factor(deriv(L,g))
  Lx=Tx[3,1]
  Ly=Ty[2,1]
  Lg=Tg[5,1]

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[3,1]
  res2g=T2[5,1]

  T11=factor(polresultant(res1g,res2g,y))
  res1y=T11[3,1]

  R(x,y,g)=(1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+2*g*y+y^2)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2g,x,VX[ix]),[4/5,1]);
           if (#VY > 0,

              for(iy=1,#Vy,
                 VG=polrootsreal(subst(subst(Lg,x,VX[ix]),y,VY[iy]),[0,1]);

                 if(#VG>0,
                    for(ig=1,#VG,
                        tot+=1;
                        value=R(VX[ix],VY[iy],VG[ig]);
                        if(value>thresold,
                            print("[x,y,g]=",[VX[ix],VY[iy],VG[ig]]);
                            print("R(x,y,g)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(2);
  /* No stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 32/27) */  
  /* y = 4/5: already discussed in greater generality (<= 1.94) */  
  /* g = 0: already discussed in greater generality (<= 1.94) */  
  /* x = 1: already discussed in greater generality (<= 1.94) */  
  /* y = 1: NULL */  
  /* g = 1: (1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-y)*(1+2*y+y^2) */
    L=(1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-y)*(1+2*y+y^2)
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))    
    Lx=Tx[4,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (it is <= 32/27) */
    /* y = 4/5: already discussed in greater generality (it is <= 1.94) */
    /* x = 1: already discussed in greater generality (it is <= 1.94) */
    /* y = 1: NULL */
.........................................................................

  /* BORDER
     g = 1: (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-y)*(1+(y*z)^2)
  */

  L= (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-y)*(1+(y*z)^2)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Lx=Tx[3,1]
  Ly=Ty[2,1]
  Lz=Tz[6,1]

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z))
  res1z=T1[4,1]
  res2z=T2[5,1]

  T11=factor(polresultant(res1z,res2z,y))
  res1y=T11[3,1]

  R(x,y,z)=(1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-y)*(1+(y*z)^2)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2z,x,VX[ix]),[4/5,1]);
           if (#VY > 0,

              for(iy=1,#Vy,
                 VG=polrootsreal(subst(subst(Lz,x,VX[ix]),y,VY[iy]),[0,1]);

                 if(#VZ>0,
                    for(iz=1,#VZ,
                        tot+=1;
                        value=R(VX[ix],VY[iy],VZ[iz]);
                        if(value>thresold,
                            print("[x,y,z]=",[VX[ix],VY[iy],VZ[iz]]);
                            print("R(x,y,z)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(2);
  /* No stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 32/27) */  
  /* y = 4/5: already discussed in greater generality (<= 1.94) */  
  /* z = 0: already discussed in greater generality (<= 18/25) */  
  /* x = 1: already discussed in greater generality (<= 1.94) */  
  /* y = 1: NULL */  
  /* z = 1: already discussed in greater generality (<= 1.94) */
.........................................................................
