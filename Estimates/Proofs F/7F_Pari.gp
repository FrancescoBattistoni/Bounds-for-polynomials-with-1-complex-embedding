
    /* z = 1: 
      BORDER
       t = 1: (1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y*a)*(1-y)*(1+y)*(1+y^2)*(1+y*a) */

  L= (1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y*a)*(1-y)*(1+y)*(1+y^2)*(1+y*a);

/* INNER */

  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Ta=factor(deriv(L,a))
  Lx=Tx[6,1];
  Ly=Ty[1,1];
  La=Ta[7,1];

  T1=factor(polresultant(Ly,Lx,x))
  T2=factor(polresultant(La,Lx,x))
  res1x=T1[6,1];
  res2x=T2[1,1];

  T11=factor(polresultant(res1x,res2x,y))
  res1y=T11[5,1];

  R(x,y,a)=(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y*a)*(1-y)*(1+y)*(1+y^2)*(1+y*a)

  search(thresold)={
    tot=0;
    VA=polrootsreal(res1y,[0,1]);
    print(#VA);
    if (#VA > 0,

        for(ia=1,#VA,
           VY = polrootsreal(subst(res1x,a,VA[ia]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VX = polrootsreal(subst(subst(Lx,y,VY[iy]),a,VA[ia]),[0,1]);

                    if(#VX>0,
                       for(ix=1,#VX,

                           tot+=1;
                           value=R(VX[ix],VY[iy],VA[ia]);
                           if(value>thresold,
                               print("[x,y,a]=",[VX[ix],VY[iy],VA[ia]]);
                               print("R(x,y,a)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(1);
  /* no stationary points */

  /* BORDER */
  /* 0 borders have already been discussed in greater generality (it is <= 1.4093) */

  /* x = 1: (1+y)*(1-y)*(1+y^2)*(1-y*a)*(1-y)*(1+y)*(1+y^2)*(1+y*a) decreasing in a, maximized for a = 0 
            (1+y)*(1-y)*(1+y^2)*(1-y)*(1+y)*(1+y^2) <= 1 */

  /* y = 1: NULL */

  /* a = 1: (1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y) */
    L=(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y)
    
    Tx=factor(deriv(L,x))
    Ty=factor(deriv(L,y))
    /* All factors in Tx are strictly positive in the open set */

    /* x = 0: (1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y) <= 1.4093 */ /* PARI my(x=0);ploth(y=0,1,(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y)) */
    /* y = 0: the function is constant = 1 */
    /* x = 1: (1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y) <= 1      */ /* PARI my(x=1);ploth(y=0,1,(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y)) */
    /* y = 1: NULL */
.........................................................................

    /* z = 1: 
      BORDER
       a = 1: (1+x*y)*(1-x*y)*(1+(x*y*t)^2)*(1-x*y*t)*(1-y)*(1+y)*(1+(y*t)^2)*(1+y*t) */

  L= (1+x*y)*(1-x*y)*(1+(x*y*t)^2)*(1-x*y*t)*(1-y)*(1+y)*(1+(y*t)^2)*(1+y*t);

/* INNER */

  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Lx=Tx[6,1];
  Ly=Ty[1,1];
  Lt=Tt[6,1];

  T1=factor(polresultant(Ly,Lx,x))
  T2=factor(polresultant(Lt,Lx,x))
  res1x=T1[6,1];
  res2x=T2[3,1];

  T11=factor(polresultant(res1x,res2x,y))
  res1y=T11[5,1];

  R(x,y,t)=(1+x*y)*(1-x*y)*(1+(x*y*t)^2)*(1-x*y*t)*(1-y)*(1+y)*(1+(y*t)^2)*(1+y*t)

  search(thresold)={
    tot=0;
    VT=polrootsreal(res1y,[0,1]);
    print(#VT);
    if (#VT > 0,

        for(it=1,#VT,
           VY = polrootsreal(subst(res1x,t,VT[it]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VX = polrootsreal(subst(subst(Lx,y,VY[iy]),t,VT[it]),[0,1]);

                    if(#VX>0,
                       for(ix=1,#VX,

                           tot+=1;
                           value=R(VX[ix],VY[iy],VT[it]);
                           if(value>thresold,
                               print("[x,y,t]=",[VX[ix],VY[iy],VT[it]]);
                               print("R(x,y,t)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(1);
  /* no stationary points */

  /* BORDER */
  /* 0 borders have already been discussed in greater generality (it is <= 1.4093) */

  /* x = 1: already discussed in greater generality (it is <= 1.4093) */

  /* y = 1: NULL */

  /* t = 1: already discussed in greater generality (it is <= 1.4093) */
.........................................................................

  /*
   BORDER 
   t = 1:
   INNER
  */

  R(x,y,z,a)=(1+x*y)*(1-x*y*z)*(1+(x*y*z)^2)*(1-x*y*z*a)*(1-y)*(1+y*z)*(1+(y*z)^2)*(1+y*z*a);

  {res1y=
  x^19 - 5/2*x^18 + 1327/400*x^17 - 2241/500*x^16 + 59541/8000*x^15 - 77149/8000*x^14 + 
      59463/8000*x^13 - 49239/8000*x^12 + 23809/4000*x^11 - 14067/4000*x^10 - 14067/4000*x^9 + 
      23809/4000*x^8 - 49239/8000*x^7 + 59463/8000*x^6 - 77149/8000*x^5 + 59541/8000*x^4 - 
      2241/500*x^3 + 1327/400*x^2 - 5/2*x + 1
  ;}
  
  {res1z=
  x^15*y^3 + 34/5*x^14*y^4 - 21/5*x^14*y^3 + 12/5*x^14*y^2 + 79/5*x^13*y^5 - 514/25*x^13*y^4 + 
      588/25*x^13*y^3 - 258/25*x^13*y^2 + 48/25*x^13*y + 64/5*x^12*y^6 - 699/25*x^12*y^5 + 
      7664/125*x^12*y^4 - 1456/25*x^12*y^3 + 3938/125*x^12*y^2 - 1056/125*x^12*y + 64/125*x^12 - 
      64/25*x^11*y^6 + 5124/125*x^11*y^5 - 2138/25*x^11*y^4 + 12354/125*x^11*y^3 - 7748/125*x^11*y^2 +
      2399/125*x^11*y - 288/125*x^11 + 704/125*x^10*y^6 - 5348/125*x^10*y^5 + 12638/125*x^10*y^4 - 
      15134/125*x^10*y^3 + 2058/25*x^10*y^2 - 739/25*x^10*y + 556/125*x^10 - 2944/125*x^9*y^6 + 
      8634/125*x^9*y^5 - 16692/125*x^9*y^4 + 3877/25*x^9*y^3 - 12502/125*x^9*y^2 + 884/25*x^9*y - 
      672/125*x^9 + 704/125*x^8*y^6 - 8634/125*x^8*y^5 + 20672/125*x^8*y^4 - 24177/125*x^8*y^3 + 
      16192/125*x^8*y^2 - 5716/125*x^8*y + 784/125*x^8 - 64/25*x^7*y^6 + 5348/125*x^7*y^5 - 
      16692/125*x^7*y^4 + 24177/125*x^7*y^3 - 3672/25*x^7*y^2 + 282/5*x^7*y - 1088/125*x^7 + 
      64/5*x^6*y^6 - 5124/125*x^6*y^5 + 12638/125*x^6*y^4 - 3877/25*x^6*y^3 + 16192/125*x^6*y^2 - 
      282/5*x^6*y + 1288/125*x^6 + 699/25*x^5*y^5 - 2138/25*x^5*y^4 + 15134/125*x^5*y^3 - 
      12502/125*x^5*y^2 + 5716/125*x^5*y - 1088/125*x^5 - 79/5*x^4*y^5 + 7664/125*x^4*y^4 - 
      12354/125*x^4*y^3 + 2058/25*x^4*y^2 - 884/25*x^4*y + 784/125*x^4 - 514/25*x^3*y^4 + 
      1456/25*x^3*y^3 - 7748/125*x^3*y^2 + 739/25*x^3*y - 672/125*x^3 + 34/5*x^2*y^4 - 588/25*x^2*y^3 
      + 3938/125*x^2*y^2 - 2399/125*x^2*y + 556/125*x^2 + 21/5*x*y^3 - 258/25*x*y^2 + 1056/125*x*y - 
      288/125*x - y^3 + 12/5*y^2 - 48/25*y + 64/125
  ;}

  {res1a=
  x^5*y^4*z^3 + 3/5*x^4*y^4*z^3 + 4/5*x^4*y^3*z^3 - 4/5*x^4*y^3*z^2 + 2/5*x^3*y^3*z^3 - 
      2/5*x^3*y^3*z^2 - 3/5*x^3*y^2*z^2 + 3/5*x^3*y^2*z - 1/5*x^2*y^2*z^2 + 1/5*x^2*y^2*z + 
      2/5*x^2*y*z - 2/5*x^2*y - 1/5*x + 1/5
  ;}

  {La=
  x*y*z*a + 1/2*x - 1/2
  ;}

  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    if (#VX > 0,
       
        print(#VX);
        for(ix=1,#VX,
           VY = polrootsreal(subst(res1z,x,VX[ix]),[0,1]);
           if (#VY > 0,
      
              for(iy=1,#VY,
                    VZ = polrootsreal(subst(subst(res1a,x,VX[ix]),y,VY[iy]),[0,3/4]);
      
                    if(#VZ>0,
                       for(iz=1,#VZ, 
                             VA=polrootsreal(subst(subst(subst(La,x,VX[ix]),y,VY[iy]),z,VZ[iz]),[0,1]);
                             
                           if(#VA>0,
                              for(ia=1,#VA, 
      
                                  tot+=1;
                                  value=R(VX[ix],VY[iy],VZ[iz],VA[ia]);
                                  if(value>thresold,
                                      print("[x,y,z,a]=",[VX[ix],VY[iy],VZ[iz],VA[ia]]);
                                      print("R(x,y,z,a)=",value)
                                    )
                                  ); 
                              );
                           ); 
                       );
                  );
               );
            );
         );   
  print("number of stationary points: ",tot)
  };
  
  search(1);
  /* no stationary points */
.........................................................................

  /*
   BORDER 
   t = 1:
   BORDER
   a = 1: (1+x*y)*(1-x*y*z)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+y*z)*(1+(y*z)^2)*(1+y*z) 
  */

  L=(1+x*y)*(1-x*y*z)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+y*z)*(1+(y*z)^2)*(1+y*z)
  
/* INNER */

  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Lx=Tx[6,1];
  Ly=Ty[3,1];
  Lz=Tz[6,1];

  T1=factor(polresultant(Ly,Lx,x))
  T2=factor(polresultant(Lz,Lx,x))
  res1x=T1[5,1];
  res2x=T2[3,1];

  T11=factor(polresultant(res1x,res2x,y))
  res1y=T11[3,1]*T11[5,1];

  R(x,y,z)=(1+x*y)*(1-x*y*z)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+y*z)*(1+(y*z)^2)*(1+y*z)

  search(thresold)={
    tot=0;
    VZ=polrootsreal(res1y,[0,1]);
    print(#VZ);
    if (#VZ > 0,

        for(iz=1,#VZ,
           VY = polrootsreal(subst(res1x,z,VZ[iz]),[0,1]);
           if (#VY > 0,

              for(iy=1,#VY,
                    VX = polrootsreal(subst(subst(Lz,y,VY[iy]),z,VZ[iz]),[0,1]);

                    if(#VX>0,
                       for(ix=1,#VX,

                           tot+=1;
                           value=R(VX[ix],VY[iy],VZ[iz]);
                           if(value>thresold,
                               print("[x,y,z]=",[VX[ix],VY[iy],VZ[iz]]);
                               print("R(x,y,z)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(0.2);
  /* 1 stationary point with value <= 1 */

  /* BORDER */
  /* 0 borders have already been discussed in greater generality (it is <= 1.4093) */

  /* x = 1: already discussed in greater generality (it is <= 1.4093) */

  /* y = 1: NULL */

  /* z = 1: already discussed in greater generality (it is <= 1.4093) */
.........................................................................

  /*
   BORDER 
   a = 1:
   INNER
  */

  R(x,y,z,t)=(1+x*y)*(1-x*y*z)*(1+(x*y*z*t)^2)*(1-x*y*z*t)*(1-y)*(1+y*z)*(1+(y*z*t)^2)*(1+y*z*t);


  {res1y=
  x^10 - 3*x^9 + 463/100*x^8 - 1543/250*x^7 + 4049/500*x^6 - 1187/125*x^5 + 4049/500*x^4 - 
      1543/250*x^3 + 463/100*x^2 - 3*x + 1
  ;}

  {res2z=
  x*y - 1/2*x + 1/2
  ;}

  {res3t=
  x*y*z + 1/2*x - 1/2
  ;}

  {Lt=
  x^3*y^5*z^5*t^5 + 5/6*x^3*y^4*z^4*t^4 + 2/3*x^3*y^3*z^3*t^3 + 1/2*x^3*y^2*z^2*t^2 - 
      5/6*x^2*y^4*z^4*t^4 - 2/3*x^2*y^3*z^3*t^3 - 1/2*x^2*y^2*z^2*t^2 - 1/3*x^2*y*z*t + 
      2/3*x*y^3*z^3*t^3 + 1/2*x*y^2*z^2*t^2 + 1/3*x*y*z*t + 1/6*x - 1/2*y^2*z^2*t^2 - 1/3*y*z*t - 1/6
  ;}

  search(thresold)={
    tot=0;
    VX=polrootsreal(res1y,[0,1]);
    if (#VX > 0,
       
        print(#VX);
        for(ix=1,#VX,
           VY = polrootsreal(subst(res2z,x,VX[ix]),[0,1]);
           if (#VY > 0,
      
              for(iy=1,#VY,
                    VZ = polrootsreal(subst(subst(res3t,x,VX[ix]),y,VY[iy]),[0,3/4]);
      
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
                                  ); 
                              );
                           ); 
                       );
                  );
               );
            );
         );   
  print("number of stationary points: ",tot)
  };
  
  search(1);
  /* no stationary points */
.........................................................................
