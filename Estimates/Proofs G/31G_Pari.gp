
/* INNER */

 R(x,y,z,t)=(1+x)*(1-x*y*z)*(1+x*y*z*t)*(1+y)*(1+y*z)*(1-y*z*t)*(1+z*t)*(1+t)
 
 {res1y=
 x^17 - 793/72*x^16 + 10091/162*x^15 - 19264/81*x^14 + 660523/972*x^13 - 5955065/3888*x^12 + 
     32717155/11664*x^11 - 99013513/23328*x^10 + 124910581/23328*x^9 - 131732371/23328*x^8 + 
     115847623/23328*x^7 - 84308761/23328*x^6 + 49992139/23328*x^5 - 5887741/5832*x^4 + 
     8492203/23328*x^3 - 550307/5832*x^2 + 5125/324*x - 329/288
 ;}

 {res1z=
 x^2*y^4 - x^2*y^3 + x^2*y^2 + 3*x*y^4 - 4/3*x*y^3 + 1/3*x*y^2 - 2/3*x*y + 2/3*x + y^4 - 1/3*y^3 + 
     4/3*y^2 - 2/3*y - 1/3
 ;}

 {res1t=
 x^2*y^3*z^2 - x^2*y*z + 3*x*y^3*z^2 + 5/3*x*y^2*z^2 + x*y^2*z + 1/3*x*y*z + 2/3*x + y^3*z^2 + 
     2/3*y^2*z^2 - y^2*z - 2/3*y*z - y - 1/3
 ;}

 {Lx=
 x^2*y^2*z^2*t + 2/3*x*y^2*z^2*t - 2/3*x*y*z*t + 2/3*x*y*z - 1/3*y*z*t + 1/3*y*z - 1/3
 ;}

  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res1z,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VZ = polrootsreal(subst(subst(res1t,x,VX[ix]),y,VY[iy]),[3/4,1]);

                    if(#VZ>0,
                       for(iz=1,#VZ,
                             VT=polrootsreal(subst(subst(subst(Lx,x,VX[ix]),y,VY[iy]),z,VZ[iz]),[0,1]);

                           if(#VT>0,
                              for(it=1,#VT,
                                  tot+=1;
                                  value=R(VX[ix],VY[iy],VZ[iz],VT[it]);
                                  if(value>thresold,
                                      print("[x,y,z,t]=",[VX[ix],VY[iy],VZ[iz],VT[it]]);
                                      print("R(x,y,z,t)=",value)
                                    )
                                  )  
                              );
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(2);
  /* No stationary points */
.........................................................................

/* x = 1: 2*(1-y*z)*(1+y*z*t)*(1+y)*(1+y*z)*(1-y*z*t)*(1+z*t)*(1+t) */
  L=2*(1-y*z)*(1+y*z*t)*(1+y)*(1+y*z)*(1-y*z*t)*(1+z*t)*(1+t) 

  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tt=factor(deriv(L,t))
  Ly=Ty[3,1];
  Lz=Tz[3,1];
  Lt=Tt[4,1];

  T1=factor(polresultant(Ly,Lt,t))
  T2=factor(polresultant(Lz,Lt,t))
  res1t=T1[6,1];
  res2t=T2[7,1];

  T11=factor(polresultant(res1t,res2t,z))
  res1z=T11[4,1];

  R(y,z,t)=2*(1-y*z)*(1+y*z*t)*(1+y)*(1+y*z)*(1-y*z*t)*(1+z*t)*(1+t)
  search(thresold)={
    tot=0;
    VY=polrootsreal(res1z,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VZ = polrootsreal(subst(res1t,y,VY[iy]),[0,1]);
           if (#VZ > 0,

              for(iz=1,#VZ,
                 VT=polrootsreal(subst(subst(Lt,y,VY[iy]),z,VZ[iz]),[0,1]);

                 if(#VT>0,
                    for(it=1,#VT,
                        tot+=1;
                        value=R(VY[iy],VZ[iz],VT[it]);
                        if(value>thresold,
                            print("[y,z,t]=",[VY[iy],VZ[iz],VT[it]]);
                            print("R(y,z,t)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* no stationary points */

   /* x = 1, y = 0: already discussed in greaster generality (<= 8) */
   ....
   /* x = 1, z = 0: already discussed in greaster generality (<= 8) */
   ....
   /* x = 1, t = 0: already discussed in greaster generality (<= 4) */
   ....
   /* x = 1, y = 1: 2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t) */
     L=2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t)
     Tz=factor(deriv(L,z))
     Tt=factor(deriv(L,t))
     Lz=Tz[3,1];
     Lt=Tt[4,1];

     T1=factor(polresultant(Lz,Lt,t))
     res1t=T1[3,1];
     /* All factors are strictly positive in th eopen set */

      /* z = 0: 2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t) <= 8    */ /* PARI my(z=0); ploth(t=0,1,2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t)) */
      /* t = 0: 2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t) <= 4    */ /* PARI my(t=0); ploth(z=0,1,2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t)) */
      /* z = 1: NULL */
      /* t = 1: 2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t) <= 8.85 */ /* PARI my(t=1); ploth(z=0,1,2*(1-z)*(1+z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t)) */
   ....
   /* x = 1, z = 1: 2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t) */
     L=2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t) 
     Ty=factor(deriv(L,y))
     Tt=factor(deriv(L,t))
     Ly=Ty[3,1];
     Lt=Tt[4,1];

     T1=factor(polresultant(Ly,Lt,t))
     res1t=T1[4,1];

     u=polrootsreal(res1t,[0,1]) /* no roots */

      /* y = 0: 2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t) <= 8    */ /* PARI my(y=0); ploth(t=0,1,2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t)) */
      /* t = 0: 2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t) <= 2.4  */ /* PARI my(t=0); ploth(y=0,1,2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t)) */
      /* y = 1: NULL */
      /* t = 1: 2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t) <= 8.85 */ /* PARI my(t=1); ploth(y=0,1,2*(1-y)*(1+y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t)) */
   ....
   /* x = 1, t = 1: 2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2 */
     L=2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2
     Ty=factor(deriv(L,y))
     Tz=factor(deriv(L,z))
     Ly=Ty[4,1];
     Lz=Tz[4,1];

     T1=factor(polresultant(Ly,Lz,z))
     res1z=T1[3,1];

     u=polrootsreal(res1z,[0,1]) /* length 1 */
     y0=u[1]
     v=polrootsreal(subst(Ly,y,y0),[0,1]) /* length 1 */
     z0=v[1]
     subst(subst(L,y,y0),z,z0)
     \\ 7.9443193401381668329773433967436422556  <----

      /* y = 0: 2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2 <= 8    */ /* PARI my(y=0); ploth(z=0,1,2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2) */
      /* z = 0: 2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2 <= 8    */ /* PARI my(z=0); ploth(y=0,1,2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2) */
      /* y = 1: 2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2 <= 8.85 */ /* PARI my(y=1); ploth(z=0,1,2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2) */
      /* z = 1: 2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2 <= 8.85 */ /* PARI my(z=1); ploth(y=0,1,2*(1-y*z)*(1+y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2) */
.........................................................................

  /* y = 1: (1+x)*(1-x*z)*(1+x*z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t) */
  L=(1+x)*(1-x*z)*(1+x*z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t)

  Tx=factor(deriv(L,x))
  Tz=factor(deriv(L,z))
  Tt=factor(deriv(L,t))
  Lx=Tx[5,1];
  Lz=Tz[3,1];
  Lt=Tt[4,1];

  T1=factor(polresultant(Lx,Lt,t))
  T2=factor(polresultant(Lz,Lt,t))
  res1t=T1[4,1];
  res2t=T2[7,1];

  T11=factor(polresultant(res1t,res2t,z))
  res1z=T11[4,1];

  R(x,z,t)=(1+x)*(1-x*z)*(1+x*z*t)*2*(1+z)*(1-z*t)*(1+z*t)*(1+t)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1z,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VZ = polrootsreal(subst(res1t,x,VX[ix]),[0,1]);
           if (#VZ > 0,

              for(iz=1,#VZ,
                 VT=polrootsreal(subst(subst(Lt,x,VX[ix]),z,VZ[iz]),[0,1]);

                 if(#VT>0,
                    for(it=1,#VT,
                        tot+=1;
                        value=R(VX[ix],VZ[iz],VT[it]);
                        if(value>thresold,
                            print("[x,z,t]=",[VX[ix],VZ[iz],VT[it]]);
                            print("R(x,z,t)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* 1 stationary point with value <= 5 */

   /* y = 1, x = 0: already discussed in greaster generality (<= 128/27) */
   ....
   /* y = 1, z = 0: already discussed in greaster generality (<= 8) */
   ....
   /* y = 1, t = 0: already discussed in greaster generality (<= 4) */
   ....
   /* y = 1, x = 1: already discussed in greaster generality (<= 8.85) */
   ....
   /* y = 1, z = 1: (1+x)*(1-x)*(1+x*t)*2*2*(1-t)*(1+t)*(1+t) */
     L=(1+x)*(1-x)*(1+x*t)*2*2*(1-t)*(1+t)*(1+t)
     Tx=factor(deriv(L,x))
     Tt=factor(deriv(L,t))
     Lx=Tx[3,1];
     Lt=Tt[4,1];

     T1=factor(polresultant(Lx,Lt,t))
     res1t=T1[3,1];

     u=polrootsreal(res1t,[0,1]) /* length 1 */
     x0=u[1]
     v=polrootsreal(subst(Lx,x,x0),[0,1]) /* length 1 */
     t0=v[1]
     subst(subst(L,x,x0),t,t0)
     \\ 4.8827530506765024137497709482063874798  <----

      /* x = 0: (1+x)*(1-x)*(1+x*t)*2*2*(1-t)*(1+t)*(1+t) <= 4.8 */ /* PARI my(x=0); ploth(t=0,1,(1+x)*(1-x)*(1+x*t)*2*2*(1-t)*(1+t)*(1+t)) */
      /* t = 0: (1+x)*(1-x)*(1+x*t)*2*2*(1-t)*(1+t)*(1+t) <= 4   */ /* PARI my(t=0); ploth(x=0,1,(1+x)*(1-x)*(1+x*t)*2*2*(1-t)*(1+t)*(1+t)) */
      /* x = 1: NULL */
      /* t = 1: NULL */
   ....
   /* y = 1, t = 1: (1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2 */
     L=(1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2
     Tx=factor(deriv(L,x))
     Tz=factor(deriv(L,z))
     Lx=Tx[3,1];
     Lz=Tz[3,1];

     T1=factor(polresultant(Lx,Lz,z))
     res1z=T1[3,1];

     u=polrootsreal(res1z,[0,1]) /* no roots */

      /* x = 0: (1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2 <= 4.8  */ /* PARI my(x=0); ploth(z=0,1,(1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2) */
      /* z = 0: (1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2 <= 8    */ /* PARI my(z=0); ploth(x=0,1,(1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2) */
      /* x = 1: (1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2 <= 8.85 */ /* PARI my(x=1); ploth(z=0,1,(1+x)*(1-x*z)*(1+x*z)*2*(1+z)*(1-z)*(1+z)*2) */
      /* z = 1: NULL */
.........................................................................

  /* z = 1: (1+x)*(1-x*y)*(1+x*y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t) */
  L=(1+x)*(1-x*y)*(1+x*y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t)

  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Lx=Tx[4,1];
  Ly=Ty[4,1];
  Lt=Tt[5,1];

  T1=factor(polresultant(Lx,Lt,t))
  T2=factor(polresultant(Ly,Lt,t))
  res1t=T1[4,1];
  res2t=T2[6,1];

  T11=factor(polresultant(res1t,res2t,y))
  res1y=T11[3,1];

  R(x,y,t)=(1+x)*(1-x*y)*(1+x*y*t)*(1+y)*(1+y)*(1-y*t)*(1+t)*(1+t)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res1t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VT=polrootsreal(subst(subst(Lt,x,VX[ix]),y,VY[iy]),[0,1]);

                 if(#VT>0,
                    for(it=1,#VT,
                        tot+=1;
                        value=R(VX[ix],VY[iy],VT[it]);
                        if(value>thresold,
                            print("[x,y,t]=",[VX[ix],VY[iy],VT[it]]);
                            print("R(x,y,t)=",value)
                          )
                       )  
                    );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* 1 stationary point with value <= 5 */

   /* z = 1, x = 0: already discussed in greaster generality (<= 128/27) */
   ....
   /* z = 1, y = 0: already discussed in greaster generality (<= 8) */
   ....
   /* z = 1, t = 0: already discussed in greaster generality (<= 4) */
   ....
   /* z = 1, x = 1: already discussed in greaster generality (<= 8.85) */
   ....
   /* z = 1, y = 1: already discussed in greaster generality (<= 8.85) */
   ....
   /* z = 1, t = 1: (1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2 */
     L=(1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2
     Tx=factor(deriv(L,x))
     Ty=factor(deriv(L,y))
     Lx=Tx[3,1];
     Ly=Ty[3,1];

     T1=factor(polresultant(Lx,Ly,y))
     res1y=T1[3,1];

     u=polrootsreal(res1y,[0,1]) /* no roots */

      /* x = 0: (1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2 <= 4.8  */ /* PARI my(x=0); ploth(y=0,1,(1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2) */
      /* y = 0: (1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2 <= 8    */ /* PARI my(y=0); ploth(x=0,1,(1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2) */
      /* x = 1: (1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2 <= 8.85 */ /* PARI my(x=1); ploth(y=0,1,(1+x)*(1-x*y)*(1+x*y)*(1+y)*(1+y)*(1-y)*2*2) */
      /* y = 1: NULL */
.........................................................................

  /* t = 1: (1+x)*(1-x*y*z)*(1+x*y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2 */
  L=(1+x)*(1-x*y*z)*(1+x*y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2

  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Lx=Tx[5,1];
  Ly=Ty[3,1];
  Lz=Tz[3,1];

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z))
  res1z=T1[4,1];
  res2z=T2[6,1];

  T11=factor(polresultant(res1z,res2z,y))
  res1y=T11[3,1]*T11[4,1];

  R(x,y,z)=(1+x)*(1-x*y*z)*(1+x*y*z)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*2
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res1z,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VZ=polrootsreal(subst(subst(Lz,x,VX[ix]),y,VY[iy]),[0,1]);

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

  search(4);
  /* 1 stationary point with value <= 6.4 */

   /* t = 1, x = 0: already discussed in greaster generality (<= 128/27) */
   ....
   /* t = 1, y = 0: already discussed in greaster generality (<= 8) */
   ....
   /* t = 1, z = 0: already discussed in greaster generality (<= 8) */
   ....
   /* t = 1, x = 1: already discussed in greaster generality (<= 8.85) */
   ....
   /* t = 1, y = 1: already discussed in greaster generality (<= 8.85) */
   ....
   /* t = 1, z = 1: already discussed in greaster generality (<= 8.85) */
.........................................................................
