
  /* z = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) */
  /* INNER */

  R(x,y,t,g)= (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t);

  {res1x=
  y^5 - 343/150*y^4 + 244/225*y^3 + 8/675*y^2 + 16/225*y + 896/6075
  ;}

  {res2t=
  x^2*y^3 - 5/6*x^2*y^2 + 5/6*x*y^2 - 2/27*x*y - 8/27*x + 8/27
  ;}

  {res3g=
  x*y*t + 2/3*x - 2/3
  ;}

  {Lg=
  x^2*y^2 - x*y^2 - 16/3*x*y*g - 16/9*x + 16/9
  ;}

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1x,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VX = polrootsreal(subst(res2t,y,VY[iy]),[0,1]);
           if (#VX > 0,

              for(ix=1,#VX,
                    VT = polrootsreal(subst(subst(res3g,x,VX[ix]),y,VY[iy]),[3/4,1]);

                    if(#VT>0,
                       for(it=1,#VT,
                             VG=polrootsreal(subst(subst(subst(Lg,x,VX[ix]),y,VY[iy]),t,VT[it]),[0,1]);

                           if(#VG>0,
                              for(ig=1,#VG,
                                  tot+=1;
                                  value=R(VX[ix],VY[iy],VT[it],VG[ig]);
                                  if(value>thresold,
                                      print("[x,y,t,g]=",[VX[ix],VY[iy],VT[it],VG[ig]]);
                                      print("R(x,y,t,g)=",value)
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

  /* z = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) 
  BORDER 
     t = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) */

  L= (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tg=factor(deriv(L,g))
  Lx=Tx[4,1]
  Ly=Ty[2,1]
  Lg=Tg[7,1]

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[3,1]
  res2g=T2[5,1]

  T11=factor(polresultant(res1g,res2g,y))
  res1y=T11[5,1]

  R(x,y,g)=(1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2g,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
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
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  /* y = 0: already discussed in greater generality (<= 2) */  
  /* g = 0: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) */  
    L=(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[2,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4))) */ 
    /* y = 1: NULL */
  ....
 
  /* x = 1: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) */  
    L=2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) 
    Ty=factor(deriv(L,y))
    Tg=factor(deriv(L,g))
    /* All factors in Tg are strictly positive in the open set */

    /* y = 0: already discussed in greater generality (<= 2) */
    /* g = 0: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(g=0);2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))) */ 
    /* y = 1: NULL */
    /* g = 1: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(g=1);2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))) */ 
  ....
  /* y = 1: NULL */  
  ....
  
  /* g = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) */  
    L=(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[5,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* length 1 */
    x0=u[1]
    v=polrootsreal(subst(Ly,x,x0),[0,1]) /* length 1 */
    y0=v[1]
    subst(subst(L,x,x0),y,y0)
    \\ 1.4425282255012222669584892119365115339  <----
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* z = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) 
  BORDER 
     g = 0: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*t) */

  L= (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*t) 
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Lx=Tx[4,1]
  Ly=Ty[2,1]
  Lt=Tt[6,1]

  T1=factor(polresultant(Lx,Lt,t))
  T2=factor(polresultant(Ly,Lt,t))
  res1t=T1[3,1]
  res2t=T2[4,1]

  T11=factor(polresultant(res1t,res2t,y))
  res1y=T11[3,1]

  R(x,y,t)=(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*t) 
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VT=polrootsreal(subst(subst(Lt,x,VX[ix]),y,VY[iy]),[3/4,1]);

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

  search(2);
  /* No stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  /* y = 0: */
  /* t = 3/4: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) */
    L=(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[2,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4))) */ 
    /* y = 1: NULL */
  ....
  /* x = 1: 2*(1+y)*(1+(y*(3/4))^2)*(1-y*(3/4)*t)*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*t) 
    the function decreases in t, so it is maximized for t = 3/4:
    2*(1+y)*(1+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4))  <= 2 */
  ....
  /* y = 1: NULL */
  ....
  /* t = 1: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)) */
    L=(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[2,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* z = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) 
  BORDER 
     x = 1: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) 
  decreasing in t, so it is maximized for t = 3/4.   
         2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) 
  decreasing in g, so it is maximized for g = 0.   
         2*(1+y)*(1+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) < = 2 */
.........................................................................

  /* z = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) 
  BORDER 
     t = 1: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) */

  L= (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tg=factor(deriv(L,g))
  Lx=Tx[4,1]
  Ly=Ty[2,1]
  Lg=Tg[7,1]

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[3,1]
  res2g=T2[5,1]

  T11=factor(polresultant(res1g,res2g,y))
  res1y=T11[4,1]

  R(x,y,g)=(1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2g,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
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
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  /* y = 0: already discussed in greater generality (<= 2) */  
  /* g = 0: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)) */ 
    L=(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[2,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+(y*(3/4))^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
  ....
  /* x = 1: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) */  
    L=2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))
    Ty=factor(deriv(L,y))
    Tg=factor(deriv(L,g))
    /* All factors in Tg have no zeros in the open set */
    
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* g = 0: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(g=0);2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
    /* g = 1: 2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(g=1);2*(1+y)*(1-2*g*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4))*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))) */ 
  ....
  /* y = 1: NULL */  
  /* g = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) */ 
    L=(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* z = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*g*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) 
  BORDER 
     g = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) */

  L= (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tt=factor(deriv((1-x*y*(3/4)*t)*(1+y*(3/4)*t),t))  /* These are the factors in L containing t: the simplification is needed otherwise PARI is not able to compute it: the excluded factors are not necessary since every zero of these factors is a zero for L */
  Lx=Tx[5,1]
  Ly=Ty[4,1]
  Lt=Tt[2,1]

  T1=factor(polresultant(Lx,Lt,t))
  T2=factor(polresultant(Ly,Lt,t))
  res1t=T1[3,1]
  res2t=T2[4,1]

  T11=factor(polresultant(res1t,res2t,y))
  res1y=T11[3,1]

  R(x,y,t)=(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*t)*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VT=polrootsreal(subst(subst(Lt,x,VX[ix]),y,VY[iy]),[3/4,1]);

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

  search(2);
  /* No stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  /* y = 0: already discussed in greater generality (<= 2) */
  /* t = 3/4: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) */
    L=(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[5,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* Length 1 */
    x0=u[1]
    v=polrootsreal(subst(Ly,x,x0),[0,1]) /* Length 1 */
    y0=v[1]
    subst(subst(L,x,x0),y,y0)
    \\ 1.4425282255012222669584892119365115339  <----
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4))) */ 
    /* y = 1: NULL */
  ....  
  /* x = 1: 2*(1+y)*(1-2*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*t)*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*t) decreasing in t, so it is maximized for t = 3/4 
            2*(1+y)*(1-2*y*(3/4)+(y*(3/4))^2)*(1-y*(3/4)*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)*(3/4)) <= 2 */
  ....
  /* y = 1: NULL */
  /* t = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) */
    L=(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1-2*x*y*(3/4)+(x*y*(3/4))^2)*(1-x*y*(3/4))*(1-y)*(1+2*y*(3/4)+(y*(3/4))^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* t = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) */
  /* INNER */

  R(x,y,z,g)= (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4));

  {res1x=
  y^5 - 171/125*y^4 - 1633/4500*y^3 + 7171/12000*y^2 + 1201/12000*y - 3791/144000
  ;}

  {res2z=
  x*y - 1/2*x + 1/2
  ;}

  {res3g=
  x^2*y^3*z^3 + 10/9*x^2*y^2*z^2 - 10/9*x*y^2*z^2 - 23/27*x*y*z + 2/9*x - 2/9
  ;}

  {Lg=
  x^2*y^2*z^2 - x*y^2*z^2 - 4*x*y*z*g - x + 1
  ;}

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1x,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VX = polrootsreal(subst(res2z,y,VY[iy]),[0,1]);
           if (#VX > 0,

              for(ix=1,#VX,
                    VZ = polrootsreal(subst(subst(res3g,x,VX[ix]),y,VY[iy]),[3/4,1]);

                    if(#VZ>0,
                       for(iz=1,#VZ,
                             VG=polrootsreal(subst(subst(subst(Lg,x,VX[ix]),y,VY[iy]),z,VZ[iz]),[0,1]);

                           if(#VG>0,
                              for(ig=1,#VG,
                                  tot+=1;
                                  value=R(VX[ix],VY[iy],VZ[iz],VG[ig]);
                                  if(value>thresold,
                                      print("[x,y,z,g]=",[VX[ix],VY[iy],VZ[iz],VG[ig]]);
                                      print("R(x,y,z,g)=",value)
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

  /* t = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) 
  BORDER 
     g = 0: (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) */
  L=(1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Lx=Tx[4,1]
  Ly=Ty[2,1]  /* Ly:= (72*z^6*y^7 + (-63*z^6 + 84*z^5)*y^6 + (-72*z^5 + 54*z^4)*y^5 + (-45*z^4 + 60*z^3)*y^4 - 48*z^3*y^3)*x^4 + ((63*z^6 - 84*z^5)*y^6 + (-54*z^6 + 144*z^5 - 96*z^4)*y^5 + (-60*z^5 + 125*z^4 - 60*z^3)*y^4 + (-36*z^4 + 96*z^3 - 64*z^2)*y^3 + (-36*z^3 + 48*z^2)*y^2)*x^3 + ((-72*z^5 + 54*z^4)*y^5 + (60*z^5 - 125*z^4 + 60*z^3)*y^4 + (64*z^4 - 96*z^3 + 36*z^2)*y^3 + (36*z^3 - 75*z^2 + 36*z)*y^2 + (32*z^2 - 24*z)*y)*x^2 + ((45*z^4 - 60*z^3)*y^4 + (-36*z^4 + 96*z^3 - 64*z^2)*y^3 + (-36*z^3 + 75*z^2 - 36*z)*y^2 + (-18*z^2 + 48*z - 32)*y + (-12*z + 16))*x + (-48*z^3*y^3 + (36*z^3 - 48*z^2)*y^2 + (32*z^2 - 24*z)*y + (12*z - 16)); */
  Lz=Tz[5,1]  /* Lz:= (27*z^5*y^5 + 30*z^4*y^4 + 18*z^3*y^3 + 18*z^2*y^2)*x^3 + (-30*z^4*y^4 - 32*z^3*y^3 - 18*z^2*y^2 - 16*z*y)*x^2 + (18*z^3*y^3 + 18*z^2*y^2 + 9*z*y + 6)*x + (-18*z^2*y^2 - 16*z*y - 6); */

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z)) /* computed with Magma, Pari needs too much memory: Factorization(Resultant(Ly,Lz,z)) */
  res1z=T1[4,1]
  res2z= x*y - 1/2*x + 1/2
  
  T11=factor(polresultant(res1z,res2z,y))
  res1y=T11[3,1]

  R(x,y,z)=(1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2z,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VT=polrootsreal(subst(subst(Lz,x,VX[ix]),y,VY[iy]),[3/4,1]);

                 if(#VT>0,
                    for(it=1,#VT,
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
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  /* y = 0: already discussed in greater generality (<= 2) */
  /* z = 3/4: already discussed in greater generality (<= 2) */
  ....
  /* x = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) */
    L=2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))
    Ty=factor(deriv(L,y))
    Tz=factor(deriv(L,z))
    Ly=Ty[3,1]
    Lz=Tz[6,1]
    
    T1=factor(polresultant(Ly,Lz,z))
    /* All factors in T1 are strictly positive in the open set */
    
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* z = 3/4: already discussed in greater generality (<= 2) */ 
    /* y = 1: NULL */
    /* z = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) <= 2.0594 */ /* PARI ploth(y=0,1,my(z=1);2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))) */ 
  ....
  /* y = 1: NULL */
  /* z = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4)) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[2,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4)) <= 2.0594 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* t = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) 
  BORDER 
     x = 1: 2*(1+y)*(1-2*g*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) maximized for g = 0
            2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) */
  L=2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Ly=Ty[3,1]
  Lz=Tz[6,1]
    
  T1=factor(polresultant(Ly,Lz,z))
  /* All factors in Tz have no zeros in the open set */

  /* y = 0: already discussed in greater generality (<= 2) */ 
  /* z = 3/4: already discussed in greater generality (<= 2) */ 
  /* y = 1: NULL */
  /* z = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) <= 2.0594 */ /* PARI ploth(y=0,1,my(z=1);2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))) */ 
.........................................................................

  /* t = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) 
  BORDER 
     z = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4)) */
  L=(1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4))
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tg=factor(deriv(L,g))
  Lx=Tx[4,1]
  Ly=Ty[2,1]  
  Lg=Tg[7,1]  

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Ly,Lg,g)) 
  res1g=T1[3,1]
  res2g=T2[5,1]
  
  T11=factor(polresultant(res1g,res2g,y))
  res1y=T11[4,1]

  R(x,y,g)=(1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4))
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2g,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
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
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */
  /* y = 0: already discussed in greater generality (<= 2) */
  /* g = 0: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4)) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[2,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4)) <= 2.0594 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
  ....
  /* x = 1: 2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4)) maximized for g = 0 
            2*(1+y)*(1+y^2)*(1-y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4)) <= 2.0594 */
  ....
  /* y = 1: NULL */
  ....
  /* g = 1: (1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4)) */
    L=(1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[5,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* length 1 */
    x0=u[1]
    v=polrootsreal(subst(Ly,x,x0),[0,1]) /* length 1 */
    y0=v[1]
    subst(subst(L,x,x0),y,y0)
    \\ 1.6266141799099942374112434696288952527  <----
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* t = 3/4: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) 
  BORDER 
     g = 1: (1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4)) */
  L=(1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4))
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Lx=Tx[5,1]
  Ly=Ty[4,1]  
  Lz=Tz[7,1]  

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z)) 
  /* All factors in T2 have no zeros in the open set */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */
  ....
  /* y = 0: already discussed in greater generality (<= 2) */
  ....
  /* z = 3/4: already discussed in greater generality (<= 2) */
  ....
  /* x = 1: 2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4)) */
    L=2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4))
    Ty=factor(deriv(L,y))
    Tz=factor(deriv(L,z))
    /* All factors in Tz have no zeros in the open set */
    
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* z = 3/4: 2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(z=3/4);2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4))) */ 
    /* y = 1: NULL */
    /* z = 1: 2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4))   <= 2 */ /* PARI ploth(y=0,1,my(z=1);2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4)))   */ 
  ....
  /* y = 1: NULL */
  /* z = 1: (1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4)) */
    L=(1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4))
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[5,1]
    Ly=Ty[4,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* length 1 */
    x0=u[1]
    v=polrootsreal(subst(Ly,x,x0),[0,1]) /* length 1 */
    y0=v[1]
    subst(subst(L,x,x0),y,y0)
    \\ 1.6266141799099942374112434696288952527 <----
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*(3/4))*(1-y)*(1+2*y+y^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
.........................................................................

  /* g = 0: (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z*t)*(1-y)*(1+(y*z)^2)*(1+y*z*t)
     BORDER
     x = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z*t)*(1-y)*(1+(y*z)^2)*(1+y*z*t) decreasing in t, so ti is maximized for t = 3/4 
            2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) 
     */
  L=2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) 
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Ly=Ty[3,1]
  Lz=Tz[6,1]
  
  T1=factor(polresultant(Ly,Lz,z))  
  /* All factors in T1 have no zeros in th eopen set */

  /* BORDER */
  /* y = 0: already discussed in greater generality (<= 2) */
  /* z = 3/4: already discussed in greater generality (<= 2) */
  /* y = 1: NULL */ 
  /* z = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4)) <= 2.0594 */ /* PARI ploth(y=0,1,my(z=1);2*(1+y)*(1+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+(y*z)^2)*(1+y*z*(3/4))) */ 
.........................................................................

  /* g = 0: (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z*t)*(1-y)*(1+(y*z)^2)*(1+y*z*t)
     BORDER
     z = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*t)*(1-y)*(1+y^2)*(1+y*t) */
  L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*t)*(1-y)*(1+y^2)*(1+y*t)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Lx=Tx[4,1]
  Ly=Ty[2,1]  
  Lt=Tt[6,1]  

  T1=factor(polresultant(Lx,Lt,t))
  T2=factor(polresultant(Ly,Lt,t))
  res1t=T1[3,1] 
  res2t=T2[4,1]
  
  T11=factor(polresultant(res1t,res2t,y))
  res1y=T11[3,1]

  R(x,y,t)=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*t)*(1-y)*(1+y^2)*(1+y*t)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VT=polrootsreal(subst(subst(Lt,x,VX[ix]),y,VY[iy]),[3/4,1]);

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

  search(2);
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */
  /* y = 0: already discussed in greater generality (<= 2) */
  /* t = 3/4: already discussed in greater generality (<= 2) */
  /* x = 1: 2*(1+y)*(1+y^2)*(1-y*t)*(1-y)*(1+y^2)*(1+y*t) decreasing in t, maximized for t = 3/4 
            2*(1+y)*(1+y^2)*(1-y*(3/4))*(1-y)*(1+y^2)*(1+y*(3/4)) <= 2.0594 */
  /* y = 1: NULL */
  /* t = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[3,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)) */ 
    /* y = 1: NULL */
.........................................................................

 /* g = 0: (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z*t)*(1-y)*(1+(y*z)^2)*(1+y*z*t)
     BORDER
    t = 1: (1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+(y*z)^2)*(1+y*z) */
  L=(1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+(y*z)^2)*(1+y*z)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Lx=Tx[4,1]
  Ly=Ty[2,1] /* Ly:= (8*z^6*y^7 + (-7*z^6 + 7*z^5)*y^6 + (-6*z^5 + 6*z^4)*y^5 + (-5*z^4 + 5*z^3)*y^4 - 4*z^3*y^3)*x^4 + ((7*z^6 - 7*z^5)*y^6 + (-6*z^6 + 12*z^5 - 6*z^4)*y^5 + (-5*z^5 + 10*z^4 - 5*z^3)*y^4 + (-4*z^4+ 8*z^3 - 4*z^2)*y^3 + (-3*z^3 + 3*z^2)*y^2)*x^3 + ((-6*z^5 + 6*z^4)*y^5 + (5*z^5 - 10*z^4 + 5*z^3)*y^4 + (4*z^4 - 8*z^3 + 4*z^2)*y^3 + (3*z^3 - 6*z^2 + 3*z)*y^2 + (2*z^2 - 2*z)*y)*x^2 + ((5*z^4 - 5*z^3)*y^4 + (-4*z^4 + 8*z^3 - 4*z^2)*y^3 + (-3*z^3 + 6*z^2 - 3*z)*y^2 + (-2*z^2 + 4*z - 2)*y + (-z +1))*x + (-4*z^3*y^3 + (3*z^3 - 3*z^2)*y^2 + (2*z^2 - 2*z)*y + (z - 1)); */ 
  Lz=Tz[5,1] /* Lz:= (6*z^5*y^5 + 5*z^4*y^4 + 4*z^3*y^3 + 3*z^2*y^2)*x^3 + (-5*z^4*y^4 - 4*z^3*y^3 - 3*z^2*y^2 - 2*z*y)*x^2 + (4*z^3*y^3 + 3*z^2*y^2 + 2*z*y + 1)*x + (-3*z^2*y^2 - 2*z*y - 1); */  

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z)) /* computed with Magma, Pari needs too much memory Factorization(Resultant(Ly,Lz,z)); */
  res1t=T1[4,1]          
  res2t=x*y - 3/4*x + 3/4
   
  T11=factor(polresultant(res1t,res2t,y)) 
  res1y=T11[2,1]

  R(x,y,z)=(1+x)*(1+x*y)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+(y*z)^2)*(1+y*z)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VZ=polrootsreal(subst(subst(Lz,x,VX[ix]),y,VY[iy]),[3/4,1]);

                 if(#Vz>0,
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
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */
  ....
  /* y = 0: already discussed in greater generality (<= 2) */
  ....
  /* z = 3/4: already discussed in greater generality (<= 2) */
  ....
  /* x = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z)*(1-y)*(1+(y*z)^2)*(1+y*z) */
    L=2*(1+y)*(1+(y*z)^2)*(1-y*z)*(1-y)*(1+(y*z)^2)*(1+y*z) 
    Ty=factor(deriv(L,y))
    Tz=factor(deriv(L,z))
    Ly=Ty[3,1]
    Lz=Tz[6,1]

    T1=factor(polresultant(Ly,Lz,z))
    /* All factors are strictly positive in the open set */
    
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* z = 3/4: already discussed in greater generality (<= 2) */ 
    /* y = 1: NULL */
    /* z = 1: 2*(1+y)*(1+(y*z)^2)*(1-y*z)*(1-y)*(1+(y*z)^2)*(1+y*z) <= 2 */ /* PARI ploth(y=0,1,my(z=1);2*(1+y)*(1+(y*z)^2)*(1-y*z)*(1-y)*(1+(y*z)^2)*(1+y*z)) */ 
  ....
  /* y = 1: NULL */
  ....
  /* z = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) 
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[3,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */ 
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)) */ 
    /* y = 1: NULL */
.........................................................................

/* x = 1: 2*(1+y)*(1-2*g*y*z+(y*z)^2)*(1-y*z*t)*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*t) decreasing in t, so it is maximized for t = 3/4
          2*(1+y)*(1-2*g*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4)) */

  L=2*(1+y)*(1-2*g*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z*(3/4))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tg=factor(deriv(L,g))
  /* All factors in Tg have no zeros in the open set */
  
  /* BORDER */
  /* y = 0: already discussed in greater generality (<= 2) */
  ....
  /* z = 3/4: already discussed in greater generality (<= 2) */
  ....
  /* g = 0: already discussed in greater generality (<= 2.0594) */
  ....
  /* y = 1: NULL */
  ....
  /* z = 1: 2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4)) */
    L=2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4))
    Ty=factor(deriv(L,y))
    Tg=factor(deriv(L,g))
    Ly=Ty[2,1]
    Lg=Tg[5,1]

    T1=factor(polresultant(Ly,Lg,g))
    /* All factors in T1 are strictly positive in the open set */
    
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* g = 0: 2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4)) <= 2.0594 */ /* PARI ploth(y=0,1,my(g=0);2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4))) */ 
    /* y = 1: NULL */
    /* g = 1: 2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(g=1);2*(1+y)*(1-2*g*y+y^2)*(1-y*(3/4))*(1-y)*(1+2*g*y+y^2)*(1+y*(3/4))) */ 
  ....
  /* g = 1: 2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4)) */
    L=2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4))
    Ty=factor(deriv(L,y))
    Tz=factor(deriv(L,z))
    /* All factors in Tz are strictly positive in the open set */
    
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* z = 3/4: already discussed in greater generality (<= 2) */ 
    /* y = 1: NULL */
    /* z = 1: 2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4)) <= 2 */ /* PARI ploth(y=0,1,my(z=1);2*(1+y)*(1-2*y*z+(y*z)^2)*(1-y*z*(3/4))*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*(3/4))) */ 
.........................................................................

  /* z = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*t)*(1-y)*(1+2*g*y+y^2)*(1+y*t) */

  /* INNER */

  R(x,y,t,g)= (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*t)*(1-y)*(1+2*g*y+y^2)*(1+y*t);

  {res1y=
  x^3 - 689/294*x^2 + 256/147*x - 125/294
  ;}
  
  {res2t=
  x^2*y^3 - 5/6*x^2*y^2 + 5/6*x*y^2 - 1/3*x*y - 1/6*x + 1/6
  ;}
  
  {res3g=
  x*y*t + 1/2*x - 1/2
  ;}
  
  {Lg=
  x^2*y^2 - x*y^2 - 4*x*y*g - x + 1
  ;}
  
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VT = polrootsreal(subst(subst(res3g,x,VX[ix]),y,VY[iy]),[3/4,1]);

                    if(#VT>0,
                       for(it=1,#VT,
                             VG=polrootsreal(subst(subst(subst(Lg,x,VX[ix]),y,VY[iy]),t,VT[it]),[0,1]);

                           if(#VG>0,
                              for(ig=1,#VG,
                                  tot+=1;
                                  value=R(VX[ix],VY[iy],VT[it],VG[ig]);
                                  if(value>thresold,
                                      print("[x,y,t,g]=",[VX[ix],VY[iy],VT[it],VG[ig]]);
                                      print("R(x,y,t,g)=",value)
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
  /* no stationary points */
.........................................................................

  /* z = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*t)*(1-y)*(1+2*g*y+y^2)*(1+y*t)
  BORDER
     t = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+2*g*y+y^2)*(1+y) */
  L= (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+2*g*y+y^2)*(1+y)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tg=factor(deriv(L,g))
  Lx=Tx[4,1]
  Ly=Ty[2,1]  
  Lg=Tg[7,1]  

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[3,1]*T1[4,1]
  res2g=T2[5,1]
  
  T11=factor(polresultant(res1g,res2g,y))
  res1y=T11[3,1]*T11[4,1]

  R(x,y,g)=(1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+2*g*y+y^2)*(1+y)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2g,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
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
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  ....
  /* y = 0: already discussed in greater generality (<= 2) */
  ....
  /* g = 0: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[3,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */ 
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) <= 2 */ /* PARI ploth(y=0,1,my(x=1);(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)) */ 
    /* y = 1: NULL */
  ....
  /* x = 1: 2*(1+y)*(1-2*g*y+y^2)*(1-y)*(1-y)*(1+2*g*y+y^2)*(1+y) decreasing in g, so it is maximized for g = 0 
            2*(1+y)*(1+y^2)*(1-y)*(1-y)*(1+y^2)*(1+y) <= 2 */
  ....
  /* y = 1: NULL */
  ....
  /* g = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[3,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */ 
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: already discussed in greater generality (<= 2.0594) */ 
    /* y = 1: NULL */
.........................................................................

  /* z = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*t)*(1-y)*(1+2*g*y+y^2)*(1+y*t) */
     BORDER
     g = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*t)*(1-y)*(1+y^2)*(1+y*t) */
  L= (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*t)*(1-y)*(1+y^2)*(1+y*t)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Lx=Tx[4,1]
  Ly=Ty[2,1]  
  Lt=Tt[6,1]  

  T1=factor(polresultant(Lx,Lt,t))
  T2=factor(polresultant(Ly,Lt,t))
  res1t=T1[3,1] 
  res2t=T2[4,1]
  
  T11=factor(polresultant(res1t,res2t,y))
  res1y=T11[3,1]

  R(x,y,t)=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y*t)*(1-y)*(1+y^2)*(1+y*t)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2t,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VT=polrootsreal(subst(subst(Lt,x,VX[ix]),y,VY[iy]),[3/4,1]);

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

  search(2);
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  ....
  /* y = 0: already discussed in greater generality (<= 2) */  
  ....
  /* t = 3/4: already discussed in greater generality (<= 2.0594) */  
  ....
  /* x = 1: already discussed in greater generality (<= 2.0594) */  
  ....
  /* y = 1: NULL */  
  ....
  /* t = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) */  
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[3,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */ 
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: already discussed in greater generality (<= 2.0594) */ 
    /* y = 1: NULL */
.........................................................................

  /* t = 1: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z) */

  R(x,y,z,g)= (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z);

  /* INNER */
  {res1y=
  x^3 - 689/294*x^2 + 256/147*x - 125/294
  ;}

  {res2z=
  x*y - 1/2*x + 1/2
  ;}

  {res3g=
  x^2*y^3*z^3 + 5/6*x^2*y^2*z^2 - 5/6*x*y^2*z^2 - 1/3*x*y*z + 1/6*x - 1/6
  ;}

  {Lg=
  x^2*y^2*z^2 - x*y^2*z^2 - 4*x*y*z*g - x + 1
  ;}
  
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2z,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VZ = polrootsreal(subst(subst(res3g,x,VX[ix]),y,VY[iy]),[3/4,1]);

                    if(#VZ>0,
                       for(iz=1,#VZ,
                             VG=polrootsreal(subst(subst(subst(Lg,x,VX[ix]),y,VY[iy]),z,VZ[iz]),[0,1]);

                           if(#VG>0,
                              for(ig=1,#VG,
                                  tot+=1;
                                  value=R(VX[ix],VY[iy],VZ[iz],VG[ig]);
                                  if(value>thresold,
                                      print("[x,y,z,g]=",[VX[ix],VY[iy],VZ[iz],VG[ig]]);
                                      print("R(x,y,z,g)=",value)
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
  /* no stationary points */
.........................................................................

  /* t = 1: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z)
  BORDER
     z = 1: (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+2*g*y+y^2)*(1+y) */
  L= (1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+2*g*y+y^2)*(1+y)
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tg=factor(deriv(L,g))
  Lx=Tx[4,1]
  Ly=Ty[2,1]  
  Lg=Tg[7,1]  

  T1=factor(polresultant(Lx,Lg,g))
  T2=factor(polresultant(Ly,Lg,g))
  res1g=T1[3,1]* T1[4,1]
  res2g=T2[5,1]
  
  T11=factor(polresultant(res1g,res2g,y))
  res1y=T11[3,1]*T11[4,1]

  R(x,y,g)=(1+x)*(1+x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+2*g*y+y^2)*(1+y)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2g,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
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
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  ....
  /* y = 0: already discussed in greater generality (<= 2) */
  ....
  /* g = 0: already discussed in greater generality (<= 2.0594) */
  ....
  /* x = 1: already discussed in greater generality (<= 2.0594) */
  ....
  /* y = 1: NULL */
  ....
  /* g = 1: (1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y) */
    L=(1+x)*(1+x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y^2)*(1+y)
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    Lx=Tx[4,1]
    Ly=Ty[3,1]

    T1=factor(polresultant(Lx,Ly,y))
    res1y=T1[3,1]
    
    u=polrootsreal(res1y,[0,1]) /* no roots */
    
    /* x = 0: already discussed in greater generality (<= 27/16) */ 
    /* y = 0: already discussed in greater generality (<= 2) */ 
    /* x = 1: already discussed in greater generality (<= 2.0594) */ 
    /* y = 1: NULL */
.........................................................................

  /* t = 1: (1+x)*(1+x*y)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*g*y*z+(y*z)^2)*(1+y*z)
  BORDER
     g = 1: (1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z) */
  L= (1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z)
  Tx=factor(deriv((1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z),x))             /* PARI needs too much memory to compute factor(deriv(L,x)) directly. This result coincides with factor(deriv(L,x)) up to factors whose zeros are also zeros for L, so that they are not important */
  Ty=factor(deriv(L,y))
  Tz=factor(deriv((1-2*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1+2*y*z+(y*z)^2)*(1+y*z),z)) /* PARI needs too much memory to compute factor(deriv(L,x)) directly. This result coincides with factor(deriv(L,z)) up to factors whose zeros are also zeros for L, so that they are not important */
  Lx=Tx[2,1]
  Ly=Ty[4,1]  
  Lz=Tz[4,1]  

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z))
  res1z=T1[3,1] 
  res2z=T2[4,1]
  
  T11=factor(polresultant(res1z,res2z,y))
  res1y=T11[3,1]

  R(x,y,z)=(1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z)
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2z,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                 VZ=polrootsreal(subst(subst(Lz,x,VX[ix]),y,VY[iy]),[3/4,1]);

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
  /* no stationary points */

  /* BORDER */
  /* x = 0: already discussed in greater generality (<= 27/16) */  
  /* y = 0: already discussed in greater generality (<= 2) */  
  /* z = 3/4: already discussed in greater generality (<= 2.0594) */  
  /* x = 1: already discussed in greater generality (<= 2.0594) */  
  /* y = 1: NULL */  
  /* z = 1: already discussed in greater generality (<= 2.0594) */  
.........................................................................

  /* g = 1: (1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z*t)*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*t) */

  /* INNER */

  R(x,y,z,t)= (1+x)*(1+x*y)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z*t)*(1-y)*(1+2*y*z+(y*z)^2)*(1+y*z*t);

  {res1y=
  x - 4/5
  ;}

  {res2z=
  x*y - 1/2*x + 1/2
  ;}

  {res3t=
  x*y*z + 1/2*x - 1/2
  ;}

  {Lt=
  x*y*z*t + 1/2*x - 1/2
  ;}
  
  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    print(#VX);
    if (#VX > 0,

        for(ix=1,#VX,
           VY = polrootsreal(subst(res2z,x,VX[ix]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VZ = polrootsreal(subst(subst(res3t,x,VX[ix]),y,VY[iy]),[3/4,1]);

                    if(#VZ>0,
                       for(iz=1,#VZ,
                             VT=polrootsreal(subst(subst(subst(Lt,x,VX[ix]),y,VY[iy]),z,VZ[iz]),[0,1]);

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
  /* no stationary points */
.........................................................................
