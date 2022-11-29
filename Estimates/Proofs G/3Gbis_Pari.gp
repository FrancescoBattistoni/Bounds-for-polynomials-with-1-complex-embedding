

/* x = 0, t = 1: (1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) */
  L=(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2);  /* times 2*sqrt(1-g^2) */

  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tg=factor(deriv(L,g)*(1-g^2)-L*g)
  Ly=Ty[3,1];
  Lz=Tz[1,1];
  Lg=Tg[2,1];

  T1=factor(polresultant(Ly,Lz,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[2,1];
  res2g=T2[4,1];

  T11=factor(polresultant(res1g,res2g,z))
  res1z=T11[4,1]*T11[8,1];

  R(y,z,g)=(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2)
  search(thresold)={
    tot=0;
    VY=polrootsreal(res1z,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VZ = polrootsreal(subst(res2g,y,VY[iy]),[0,1]);
           if (#VZ > 0,

              for(iz=1,#VZ,
                 VG=polrootsreal(subst(subst(Lg,y,VY[iy]),z,VZ[iz]),[0,1]);

                 if(#VG>0,
                    for(ig=1,#VG,
                        tot+=1;
                        value=R(VY[iy],VZ[iz],VG[ig]);
                        if(value>thresold,
                            print("[y,z,g]=",[VY[iy],VZ[iz],VG[ig]]);
                            print("R(y,z,g)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(3.8);
  /* 3 stationary points with value <= 3.8 (the value is small, we do not need to check the condition y*z > 1/10) */

   /* x = 0, t = 1:
     BORDER of (y,z), case 1: y in [1/10,1], z = 1:
     (1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2) */
     L=(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g);  /* times 2*sqrt(1-g^2)  */
     Ty=factor(deriv(L,y))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Ly=Ty[3,1];
     Lg=Tg[4,1];

     T1=factor(polresultant(Ly,Lg,g))
     res1g=T1[3,1];

     u=polrootsreal(res1g,[1/10,1]) /* no roots */

      /* y = 1/10: (1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2) <= 7.32  */ /* PARI my(y=1/10); ploth(g=0   ,1,(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2)) */
      /* g = 0   : (1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2) <= 7.28  */ /* PARI my(g=0)   ; ploth(y=1/10,1,(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2)) */
      /* y = 1   : NULL */
      /* g = 1   : NULL */
   ....
   /* x = 0, t = 1:
     BORDER of (y,z), case 2: y = 1, z in [1/10,1] 
     2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) */
     L=2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2);  /* times 2*sqrt(1-g^2) */
     Tz=factor(deriv(L,z))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Lz=Tz[2,1];
     Lg=Tg[3,1];

     T1=factor(polresultant(Lz,Lg,g))
     res1g=T1[5,1];

     u=polrootsreal(res1g,[1/10,1]) /* no roots */

      /* z = 1/10: 2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) <= 4.8  */ /* PARI my(z=1/10); ploth(g=0   ,1,2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2)) */
      /* g = 0   : 2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) <= 3.7  */ /* PARI my(g=0)   ; ploth(z=1/10,1,2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2)) */
      /* z = 1   : NULL */
      /* g = 1   : NULL */
   ....
   /* x = 0, t = 1:
     BORDER of (y,z), case 3: z=1/10/y and y in [1/10,1] 
     (1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) */
     L=(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2);  /* times 2*sqrt(1-g^2) */
     Ty=factor(deriv(L,y))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Ly=Ty[2,1];
     Lg=Tg[2,1];

     T1=factor(polresultant(Ly,Lg,g))
     res1g=T1[4,1]*T1[6,1]*T1[9,1];

     u=polrootsreal(res1g,[1/10,1]) /* length 3 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\ 0  (this happens because g0 = 1) 
     y0=u[2]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\  2.5331568993570300804988465402802940360  <---- 
     y0=u[3]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\  2.5973051334026966805439712237843690405  <---- 

      /* y = 1/10: (1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.32 */ /* PARI my(y=1/10); ploth(g=0   ,1,(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* g = 0   : (1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.28 */ /* PARI my(g=0)   ; ploth(y=1/10,1,(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* y = 1   : (1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 4.8  */ /* PARI my(y=1)   ; ploth(g=0   ,1,(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* g = 1   : NULL */
   ....
   /* x = 0, t = 1, g = 0: (1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 */
     L=(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2
     Ty=factor(deriv(L,y))
     Tz=factor(deriv(L,z))
     Ly=Ty[2,1];
     Lz=Tz[2,1];

     T1=factor(polresultant(Ly,Lz,z))
     res1z=T1[3,1];

     u=polrootsreal(res1z,[0,1]) /* length 1 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     z0=v[1]
     y0*z0  /* TEST: is it >= 1/10? */
     \\ 0.46666666666666666666666666666666666667
     subst(subst(L,y,y0),z,z0)
     \\ 3.1068918518518518518518518518518518519  <----

     /* BORDER case 1: y in [1/10,1], z = 1      : (1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 <= 7.28 */ /* PARI my(z=1); ploth(y=1/10,1,(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2)     */
     /* BORDER case 2: y = 1, z in [1/10,1]      : (1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 <= 3.7  */ /* PARI my(y=1); ploth(z=1/10,1,(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2)     */
     /* BORDER case 3: z=1/10/y and y in [1/10,1]: (1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 <= 7.28 */ /* PARI ploth(y=1/10,1,my(z=1/10/y);(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2) */
   ....
   /* x = 0, t = 1, g = 1: the function is zero */
.........................................................................

/* x = 1, t = 0: (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*2*sqrt(1-g^2) */
  L=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2);  /* times 2*sqrt(1-g^2) */

  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tg=factor(deriv(L,g)*(1-g^2)-L*g)
  Ly=Ty[2,1];
  Lz=Tz[1,1];
  Lg=Tg[2,1];

  T1=factor(polresultant(Ly,Lz,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[2,1];
  res2g=T2[6,1];

  T11=factor(polresultant(res1g,res2g,z))
  res1z=T11[4,1]*T11[7,1]*T11[8,1];

  R(y,z,g)=(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2)
  search(thresold)={
    tot=0;
    VY=polrootsreal(res1z,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VZ = polrootsreal(subst(res2g,y,VY[iy]),[0,1]);
           if (#VZ > 0,

              for(iz=1,#VZ,
                 VG=polrootsreal(subst(subst(Lg,y,VY[iy]),z,VZ[iz]),[0,1]);

                 if(#VG>0,
                    for(ig=1,#VG,
                        tot+=1;
                        value=R(VY[iy],VZ[iz],VG[ig]);
                        if(value>thresold,
                            print("[y,z,g]=",[VY[iy],VZ[iz],VG[ig]]);
                            print("R(y,z,g)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(3.8);
  /* 4 stationary points with value <= 4.6 (the value is small, we do not need to check the condition y*z > 1/10) */

   /* x = 1, t = 0:
     BORDER of (y,z), case 1: y in [1/10,1], z = 1:               
     (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*sqrt(1-g^2) */

     L=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g);  /* times 2*sqrt(1-g^2)  */
     Ty=factor(deriv(L,y))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Ly=Ty[2,1];
     Lg=Tg[3,1];

     T1=factor(polresultant(Ly,Lg,g))
     res1g=T1[5,1];

     u=polrootsreal(res1g,[1/10,1]) /* no roots */

      /* y = 1/10: (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*sqrt(1-g^2) <= 4.8 */ /* PARI my(y=1/10); ploth(g=0   ,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*sqrt(1-g^2)) */
      /* g = 0   : (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*sqrt(1-g^2) <= 3.7 */ /* PARI my(g=0)   ; ploth(y=1/10,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*sqrt(1-g^2)) */
      /* y = 1   : NULL */
      /* g = 1   : NULL */
   ....
   /* x = 1, t = 0:
     BORDER of (y,z), case 2: y = 1, z in [1/10,1] 
     2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*2*sqrt(1-g^2) */
     L=2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2);  /* times 2*sqrt(1-g^2) */
     Tz=factor(deriv(L,z))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Lz=Tz[3,1];
     Lg=Tg[4,1];

     T1=factor(polresultant(Lz,Lg,g))
     res1g=T1[3,1];

     u=polrootsreal(res1g,[1/10,1]) /* no roots */

      /* z = 1/10: 2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*2*sqrt(1-g^2) <= 7.32 */ /* PARI my(z=1/10); ploth(g=0   ,1,2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*2*sqrt(1-g^2)) */
      /* g = 0   : 2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*2*sqrt(1-g^2) <= 7.28 */ /* PARI my(g=0)   ; ploth(z=1/10,1,2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*2*sqrt(1-g^2)) */
      /* z = 1   : NULL */
      /* g = 1   : NULL */
   ....
   /* x = 1, t = 0:
     BORDER of (y,z), case 3: z=1/10/y and y in [1/10,1] 
     (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) */
     L=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2);  /* times 2*sqrt(1-g^2) */
     Ty=factor(deriv(L,y))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Ly=Ty[2,1];
     Lg=Tg[2,1];

     T1=factor(polresultant(Ly,Lg,g))
     res1g=T1[4,1]*T1[7,1]*T1[9,1];

     u=polrootsreal(res1g,[1/10,1]) /* length 4 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* no roots */
     y0=u[2]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\  2.6003669535187453872447605959081846947  <---- 
     y0=u[3]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\  2.5973051334026966805439712237843690405  <---- 
     y0=u[4]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\  2.5331568993570300804988465402802940360  <---- 

      /* y = 1/10: (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 4.8  */ /* PARI my(y=1/10); ploth(g=0   ,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* g = 0   : (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.28 */ /* PARI my(g=0)   ; ploth(y=1/10,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* y = 1   : (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.32 */ /* PARI my(y=1)   ; ploth(g=0   ,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* g = 1   : NULL */
   ....

   /* x = 1, t = 0, g = 0: (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2 */
     L=(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2
     Ty=factor(deriv(L,y))
     Tz=factor(deriv(L,z))
     Ly=Ty[3,1];
     Lz=Tz[1,1];

     T1=factor(polresultant(Ly,Lz,z))
     res1z=T1[3,1];

     u=polrootsreal(res1z,[0,1]) /* length 1 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     z0=v[1]
     y0*z0  /* TEST: is it >= 1/10? */
     \\ 0.46666666666666666666666666666666666667
     subst(subst(L,y,y0),z,z0)
     \\ 3.1068918518518518518518518518518518519  <----

     /* BORDER case 1: y in [1/10,1], z = 1      : (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2 <= 3.7  */ /* PARI my(z=1); ploth(y=1/10,1,(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2)     */
     /* BORDER case 2: y = 1, z in [1/10,1]      : (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2 <= 7.28 */ /* PARI my(y=1); ploth(z=1/10,1,(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2)     */
     /* BORDER case 3: z=1/10/y and y in [1/10,1]: (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2 <= 7.28 */ /* PARI ploth(y=1/10,1,my(z=1/10/y);(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*2) */
   ....
   /* x = 1, t = 0, g = 1: the function is zero */
.........................................................................

/* x = 1, t = 1: (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) */
  L=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2);  /* times 2*sqrt(1-g^2) */

  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tg=factor(deriv(L,g)*(1-g^2)-L*g)
  Ly=Ty[3,1];
  Lz=Tz[1,1];
  Lg=Tg[3,1];

  T1=factor(polresultant(Ly,Lz,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[4,1];
  res2g=T2[5,1];

  T11=factor(polresultant(res1g,res2g,z))
  /* All factors are strictly positive in the open set */

   /* x = 1, t = 1:
     BORDER of (y,z), case 1: y in [1/10,1], z = 1:               
     (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2) */
     L=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g);  /* times 2*sqrt(1-g^2)  */
     Ty=factor(deriv(L,y))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Ly=Ty[3,1];
     Lg=Tg[5,1];

     T1=factor(polresultant(Ly,Lg,g))
     /* All factors are strictly positive in  the open set */

      /* y = 1/10: (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(y=1/10); ploth(g=0   ,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2)) */
      /* g = 0   : (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(g=0)   ; ploth(y=1/10,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y)*2*(1+g)*2*(1-g)*2*sqrt(1-g^2)) */
      /* y = 1   : NULL */
      /* g = 1   : NULL */
   ....
   /* x = 1, t = 1:
     BORDER of (y,z), case 2: y = 1, z in [1/10,1] 
     2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) */                       
     L=2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2);  /* times 2*sqrt(1-g^2) */
     Tz=factor(deriv(L,z))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Lz=Tz[3,1];
     Lg=Tg[5,1];

     T1=factor(polresultant(Lz,Lg,g))
     /* All factors are strictly positive in the open set */

      /* z = 1/10: 2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(z=1/10); ploth(g=0   ,1,2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2)) */
      /* g = 0   : 2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(g=0)   ; ploth(z=1/10,1,2*(1-g)*2*(1+g)*(1-z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2)) */
      /* z = 1   : NULL */
      /* g = 1   : NULL */
   ....
   /* x = 1, t = 1:
     BORDER of (y,z), case 3: z=1/10/y and y in [1/10,1] 
     (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) */
     L=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2);  /* times 2*sqrt(1-g^2) */
     Ty=factor(deriv(L,y))
     Tg=factor(deriv(L,g)*(1-g^2)-L*g)
     Ly=Ty[2,1]*Ty[4,1];
     Lg=Tg[3,1];

     T1=factor(polresultant(Ly,Lg,g))
     res1g=T1[4,1]*T1[6,1];

     u=polrootsreal(res1g,[1/10,1]) /* length 2 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 but the root is 1: we discard it */
     y0=u[2]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     g0=v[1]
     subst(subst(L,y,y0),g,g0)*2*sqrt(1-g0^2)
     \\ 1.1905507630552942004048190130837329859  <---- 

      /* y = 1/10: (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(y=1/10); ploth(g=0   ,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* g = 0   : (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(g=0)   ; ploth(y=1/10,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* y = 1   : (1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2) <= 7.345 */ /* PARI my(y=1)   ; ploth(g=0   ,1,(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*(1/10/y))*(1+2*(1/10/y)*g+(1/10/y)^2)*(1-2*(1/10/y)*g+(1/10/y)^2)*2*sqrt(1-g^2)) */
      /* g = 1   : NULL */
   ....

   /* x = 1, t = 1, g = 0: (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 */
     L=(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2
     Ty=factor(deriv(L,y))
     Tz=factor(deriv(L,z))
     Ly=Ty[3,1];
     Lz=Tz[2,1];

     T1=factor(polresultant(Ly,Lz,z))
     res1z=T1[3,1];

     u=polrootsreal(res1z,[0,1]) /* length 1 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     z0=v[1]
     y0*z0  /* TEST: is it >= 1/10? */
     \\ 0.6000000000000000000000000000000000000
     subst(subst(L,y,y0),z,z0)
     \\ 5.2428800000000000000000000000000000000  <----

     /* BORDER case 1: y in [1/10,1], z = 1      : (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 <= 7.345 */ /* PARI my(z=1); ploth(y=1/10,1,(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2)     */
     /* BORDER case 2: y = 1, z in [1/10,1]      : (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 <= 7.345 */ /* PARI my(y=1); ploth(z=1/10,1,(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2)     */
     /* BORDER case 3: z=1/10/y and y in [1/10,1]: (1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2 <= 7.345 */ /* PARI ploth(y=1/10,1,my(z=1/10/y);(1+y^2)*(1+y^2)*(1-y*z)*(1+z^2)*(1+z^2)*2) */
   ....
   /* x = 1, t = 1, g = 1: the function is zero */
.........................................................................
