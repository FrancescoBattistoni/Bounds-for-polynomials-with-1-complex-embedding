
/* IMPORTANT:

   Some computations in this file require data and functions which are stored in separated 
   files located into the directory "/Estimates/Proofs G/auxiliary_files_23G-B_Pari".
   To import these data we use the Pari instruction 
   read("ROOT/Estimates/Proofs G/auxiliary_files_23G-B_Pari/name of the file needed.gp")
   In order to make Pari execute correctly this command you have to substitute ROOT with 
   the path of the directory "/Estimates" in your PC.
*/

/*
     INNER 
*/

  R(y,z,t,a,b) = (1+(1/2)*y*z*t*a*b)*(1+y)*(1+y*z*t)*(1-y*z*t*a)*(1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-B_Pari/auxiliary_23G-B_Pari_Q1.gp")
  
  search(thresold)={
    tot=0;
    VY=polrootsreal(res1z,[0,1]);
    if (#VY > 0,
       
        print(#VY);
        for(iy=1,#VY,
           VZ = polrootsreal(subst(res2t,y,VY[iy]),[0,1]);
           if (#VZ > 0,
      
              for(iz=1,#VZ,
                    VT = polrootsreal(subst(subst(res1a,y,VY[iy]),z,VZ[iz]),[0,1]);
      
                    if(#VT>0,
                       for(it=1,#VT, 
                             VA=polrootsreal(subst(subst(subst(res1b,y,VY[iy]),z,VZ[iz]),t,VT[it]),[0,1]);
                             
                           if(#VA>0,
                              for(ia=1,#VA, 
                                    VB=polrootsreal(subst(subst(subst(subst(Ly,y,VY[iy]),z,VZ[iz]),t,VT[it]),a,VA[ia]),[0,1]);
      
                                  if(#VB>0,
                                       for(ib=1,#VB,
                                           tot+=1;
                                           value=R(VY[iy],VZ[iz],VT[it],VA[ia],VB[ib]);
                                           if(value>thresold,
                                               print("[y,z,t,a,b]=",[VY[iy],VZ[iz],VT[it],VA[ia],VB[ib]]);
                                               print("R(y,z,t,a,b)=",value)
                                             )
                                          )     
                                      );
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
  
  search(5);
  /* no stationary points */
.........................................................................

  /* BORDER
     y = 0: (1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b) =
            (1+z*t*a)*(1-z*t*a*b)*(1+b)          *(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b) <= (by 9A)
            2                                    *(1+a*b+(a*b)^2)*(1+a)*(1-a*b) =
            2                                    *(1+a*b)^2*(1+a)*(1-a*b) =
            2                                    *(1-(a*b)^2)*(1+a)*(1+a*b) <=
            2                                    *1          *2    *2 <= 8 */
.........................................................................

  /* BORDER
     z = 0: (1+y)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b) maximized for y = t = 1
            2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */
            
  L= 2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  La=Ta[2,1]
  Lb=Tb[2,1]

  T1=factor(polresultant(La,Lb,b))
  res1b=T1[3,1]

  u=polrootsreal(res1b,[0,1]) /* length 1 */
  a0=u[1]
  v=polrootsreal(subst(La,a,a0),[0,1]) /* length 1 */
  b0=v[1]
  subst(subst(L,a,a0),b,b0)  
  \\ 5.1002696253849125420793174817443455340  <----
  
  /* a = 0: 2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 4   */ /* PARI my(a=0);ploth(b=0,1,2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
  /* b = 0: 2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 4   */ /* PARI my(b=0);ploth(a=0,1,2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
  /* a = 1: 2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 5.3 */ /* PARI my(a=1);ploth(b=0,1,2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
  /* b = 1: 2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 5.3 */ /* PARI my(b=1);ploth(a=0,1,2*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
.........................................................................

  /* BORDER
     y = 1: (1+(1/2)*z*t*a*b)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER
  */   

  R(z,t,a,b) = (1+(1/2)*z*t*a*b)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-B_Pari/auxiliary_23G-B_Pari_Q2.gp")
  
  search(thresold)={
  tot=0;
  VZ=polrootsreal(res1t,[0,1]);
  print(#VZ);
  if (#VZ > 0,

      for(iz=1,#VZ,
         VT = polrootsreal(subst(res1a,z,VZ[iz]),[0,1]);
         if (#VT > 0,

            for(it=1,#VT,
                  VA = polrootsreal(subst(subst(res1b,z,VZ[iz]),t,VT[it]),[0,1]);

                  if(#VA>0,
                     for(ia=1,#VA,
                           VB=polrootsreal(subst(subst(subst(Lz,z,VZ[iz]),t,VT[it]),a,VA[ia]),[0,1]);

                         if(#VB>0,
                            for(ib=1,#VB,
                                tot+=1;
                                value=R(VZ[iz],VT[it],VA[ia],VB[ib]);
                                if(value>thresold,
                                    print("[z,t,a,b]=",[VZ[iz],VT[it],VA[ia],VB[ib]]);
                                    print("R(z,t,a,b)=",value)
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

  search(4);
  /* no stationary points */
.........................................................................

  /* y = 1:
  BORDER

  z = 1: (1+(1/2)*t*a*b)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b) */

  L= (1+(1/2)*t*a*b)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Tt=factor(deriv(L,t))
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  Lt=Tt[4,1]
  La=Ta[3,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(La,Lt,b))
  T2=factor(polresultant(Lb,Lt,b))
  res1b=T1[5,1]
  res2b=T2[4,1]

  T11=factor(polresultant(res1b,res2b,t))
  res1t=T11[4,1]*T11[5,1]

  R(t,a,b)=(1+(1/2)*t*a*b)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)

  search(thresold)={
    tot=0;
    VA=polrootsreal(res1t,[0,1]);
    print(#VA);
    if (#VA > 0,

        for(ia=1,#VA,
           VT = polrootsreal(subst(res1b,a,VA[ia]),[0,1]);
           if (#VT > 0,

              for(it=1,#VT,
                    VB = polrootsreal(subst(subst(Lt,t,VT[it]),a,VA[ia]),[0,1]);

                    if(#VB>0,
                       for(ib=1,#VB,

                           tot+=1;
                           value=R(VT[it],VA[ia],VB[ib]);
                           if(value>thresold,
                               print("[t,a,b]=",[VT[it],VA[ia],VB[ib]]);
                               print("R(t,a,b)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* 4 stationary points with value <= 4.8 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* t = 1: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */
    L=(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
    Ta=factor(deriv(L,a))
    Tb=factor(deriv(L,b))
    La=Ta[4,1]
    Lb=Tb[4,1]

    T1=factor(polresultant(La,Lb,b))
    res1b=T1[3,1]

    u=polrootsreal(res1b,[0,1]) /* no roots */

    /* a = 0: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8    */ /* PARI my(a=0);ploth(b=0,1,(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* b = 0: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 4.75 */ /* PARI my(b=0);ploth(a=0,1,(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* a = 1: NULL */
    /* b = 1: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8.23 */ /* PARI my(b=1);ploth(a=0,1,(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */

  /* a = 1: (1+(1/2)*t*b)*2*(1+t)*(1-t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) =
          4*(1+(1/2)*t*b)*(1-t*b)     *(1+t)*(1-t)*(1+t)*        (1+t*b+(t*b)^2)*(1-b)*(1+b) <=
          4*1                         *(32/27)          *        (1+b+b^2)*(1-b)*(1+b) <=
          4*1                         *(32/27)          *        (27/16) <= 8 */

  /* b = 1: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */
    L=(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
    
    Tt=factor(deriv(L,t))
    Ta=factor(deriv(L,a))
    Lt=Tt[4,1]
    La=Ta[3,1]

    T1=factor(polresultant(Lt,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* t = 0: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 4    */ /* PARI my(t=0);ploth(a=0,1,(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(t=0,1,(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* t = 1: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8.23 */ /* PARI my(t=1);ploth(a=0,1,(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */
.........................................................................

  /* y = 1:
  BORDER

  t = 1: (1+(1/2)*z*a*b)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */

  L= (1+(1/2)*z*a*b)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Tz=factor(deriv(L,z))
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  Lz=Tz[5,1]
  La=Ta[3,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(Lz,Lb,b))
  T2=factor(polresultant(La,Lb,b))
  res1b=T1[4,1]
  res2b=T2[9,1]

  T11=factor(polresultant(res1b,res2b,a))
  res1a=T11[5,1]*T11[6,1]

  R(z,a,b)=(1+(1/2)*z*a*b)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)

  search(thresold)={
    tot=0;
    VZ=polrootsreal(res1a,[0,1]);
    print(#VZ);
    if (#VZ > 0,

        for(iz=1,#VZ,
           VA = polrootsreal(subst(res1b,z,VZ[iz]),[0,1]);
           if (#VA > 0,

              for(ia=1,#VA,
                    VB = polrootsreal(subst(subst(Lb,z,VZ[iz]),a,VA[ia]),[0,1]);

                    if(#VB>0,
                       for(ib=1,#VB,

                           tot+=1;
                           value=R(VZ[iz],VA[ia],VB[ib]);
                           if(value>thresold,
                               print("[z,a,b]=",[VZ[iz],VA[ia],VB[ib]]);
                               print("R(z,a,b)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(6);
  /* 2 stationary points with value <= 6.3 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* z = 1: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */
    L=(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
    Ta=factor(deriv(L,a))
    Tb=factor(deriv(L,b))
    La=Ta[4,1]
    Lb=Tb[4,1]

    T1=factor(polresultant(La,Lb,b))
    res1b=T1[3,1]

    u=polrootsreal(res1b,[0,1]) /* no roots */

    /* a = 0: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8    */ /* PARI my(a=0);ploth(b=0,1,(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* b = 0: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 4.75 */ /* PARI my(b=0);ploth(a=0,1,(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* a = 1: NULL */
    /* b = 1: (1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8.23 */ /* PARI my(b=1);ploth(a=0,1,(1+(1/2)*a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */

  /* a = 1: (1+(1/2)*z*b)*2*(1+z)*(1-z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b) <= 
          4*(1+(1/2)*z*b)*(1-z*b)      *(1+z)*(1-z)*(1+z)    *(1+b+b^2)*(1-b)*(1+b) <=
          4*1                          *(32/27)              *(27/16) <= 8 */
  */

  /* b = 1: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 */
    L=(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2
    
    Tz=factor(deriv(L,z))
    Ta=factor(deriv(L,a))
    Lz=Tz[5,1]
    La=Ta[3,1]

    T1=factor(polresultant(Lz,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* z = 0: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 5.3  */ /* PARI my(z=0);ploth(a=0,1,(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(z=0,1,(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* z = 1: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.23 */ /* PARI my(z=1);ploth(a=0,1,(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */
.........................................................................

  /* y = 1:
  BORDER

  a = 1: (1+(1/2)*z*t*b)*2*(1+z*t)*(1-z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) */

  L= (1+(1/2)*z*t*b)*2*(1+z*t)*(1-z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
  Tz=factor(deriv(L,z))
  Tt=factor(deriv(L,t))
  Tb=factor(deriv(L,b))
  Lz=Tz[6,1]
  Lt=Tt[4,1]
  Lb=Tb[3,1]

  T1=factor(polresultant(Lt,Lz,b))
  T2=factor(polresultant(Lb,Lz,b))
  res1b=T1[6,1]
  res2b=T2[5,1]

  T11=factor(polresultant(res1b,res2b,t))
  res1t=T11[2,1]*T11[3,1]

  R(z,t,b)=(1+(1/2)*z*t*b)*2*(1+z*t)*(1-z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)

  search(thresold)={
    tot=0;
    VZ=polrootsreal(res1t,[0,1]);
    print(#VZ);
    if (#VZ > 0,

        for(iz=1,#VZ,
           VT = polrootsreal(subst(res1b,z,VZ[iz]),[0,1]);
           if (#VT > 0,

              for(it=1,#VT,
                    VB = polrootsreal(subst(subst(Lz,z,VZ[iz]),t,VT[it]),[0,1]);

                    if(#VB>0,
                       for(ib=1,#VB,

                           tot+=1;
                           value=R(VZ[iz],VT[it],VB[ib]);
                           if(value>thresold,
                               print("[z,t,b]=",[VZ[iz],VT[it],VB[ib]]);
                               print("R(z,t,b)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(6);
  /* no stationary points */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* z = 1: (1+(1/2)*t*b)*2*(1+t)*(1-t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) =
          4*(1+t)*(1-t)*(1+t)*(1+(1/2)*t*b)*(1-t*b)*(1+t*b+(t*b)^2)*(1-b)*(1+b) <=
          4*(1+t)*(1-t)*(1+t)*(1+t*b)*(1-t*b)*(1+t*b+(t*b)^2)*(1-b)*(1+b) <=
          4*(32/16)          *(27/16)                        *1 <= 8 */

  /* t = 1: (1+(1/2)*z*b)*2*(1+z)*(1-z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b) =
          4*(1+(1/2)*z*b)*(1-z*b)*(1+z)*(1-z)*(1+z)*(1+b+b^2)*(1-b)*(1+b) <=
          4*1                    *(32/27)    *(27/16) <= 8 */ 

  /* b = 1: NULL */
.........................................................................

  /* y = 1:
  BORDER

  b = 1: (1+(1/2)*z*t*a)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */

  L= (1+(1/2)*z*t*a)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
  Tz=factor(deriv(L,z))
  Tt=factor(deriv(L,t))
  Ta=factor(deriv(L,a))
  Lz=Tz[6,1]
  Lt=Tt[4,1]
  La=Ta[3,1]

  T1=factor(polresultant(Lt,Lz,a))
  T2=factor(polresultant(La,Lz,a))
  res1a=T1[4,1]
  res2a=T2[4,1]

  T11=factor(polresultant(res1a,res2a,t))
  res1t=T11[1,1]*T11[2,1]

  R(z,t,a)=(1+(1/2)*z*t*a)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2

  search(thresold)={
    tot=0;
    VZ=polrootsreal(res1t,[0,1]);
    print(#VZ);
    if (#VZ > 0,

        for(iz=1,#VZ,
           VT = polrootsreal(subst(res1a,z,VZ[iz]),[0,1]);
           if (#VT > 0,

              for(it=1,#VT,
                    VA = polrootsreal(subst(subst(Lz,z,VZ[iz]),t,VT[it]),[0,1]);

                    if(#VA>0,
                       for(ia=1,#VA,

                           tot+=1;
                           value=R(VZ[iz],VT[it],VA[ia]);
                           if(value>thresold,
                               print("[z,t,a]=",[VZ[iz],VT[it],VA[ia]]);
                               print("R(z,t,a)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* No stationary points */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* z = 1: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */
    L=(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
    Tt=factor(deriv(L,t))
    Ta=factor(deriv(L,a))
    Lt=Tt[4,1]
    La=Ta[3,1]

    T1=factor(polresultant(Lt,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* t = 0: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 4    */ /* PARI my(t=0);ploth(a=0,1,(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(t=0,1,(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* t = 1: (1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8.23 */ /* PARI my(t=1);ploth(a=0,1,(1+(1/2)*t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */

  /* t = 1: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 */
    L=(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2
    Tz=factor(deriv(L,z))
    Ta=factor(deriv(L,a))
    Lz=Tz[5,1]
    La=Ta[3,1]

    T1=factor(polresultant(Lz,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* z = 0: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 5.3  */ /* PARI my(z=0);ploth(a=0,1,(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(z=0,1,(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* z = 1: (1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.23 */ /* PARI my(z=1);ploth(a=0,1,(1+(1/2)*z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */

  /* a = 1: NULL */
.........................................................................

  /* BORDER
     z = 1: (1+(1/2)*y*t*a*b)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER
  */   

  R(y,t,a,b) = (1+(1/2)*y*t*a*b)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-B_Pari/auxiliary_23G-B_Pari_Q3.gp")
  
  search(thresold)={
  tot=0;
  VY=polrootsreal(res1t,[0,1]);
  print(#VY);
  if (#VY > 0,

      for(iy=1,#VY,
         VT = polrootsreal(subst(res1a,y,VY[iy]),[0,1]);
         if (#VT > 0,

            for(it=1,#VT,
                  VA = polrootsreal(subst(subst(res1b,y,VY[iy]),t,VT[it]),[0,1]);

                  if(#VA>0,
                     for(ia=1,#VA,
                           VB=polrootsreal(subst(subst(subst(Ly,y,VY[iy]),t,VT[it]),a,VA[ia]),[0,1]);

                         if(#VB>0,
                            for(ib=1,#VB,
                                tot+=1;
                                value=R(VY[iy],VT[it],VA[ia],VB[ib]);
                                if(value>thresold,
                                    print("[y,t,a,b]=",[VY[iy],VT[it],VA[ia],VB[ib]]);
                                    print("R(y,t,a,b)=",value)
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

  search(4);
  /* 2 stationary points with value <= 4.6 */
.........................................................................

  /* z = 1:
  BORDER

  t = 1: (1+(1/2)*y*a*b)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */

  L= (1+(1/2)*y*a*b)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Ty=factor(deriv(L,y))
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  Ly=Ty[6,1]
  La=Ta[5,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(La,Ly,b))
  T2=factor(polresultant(Lb,Ly,b))
  res1b=T1[3,1]
  res2b=T2[2,1]

  T11=factor(polresultant(res1b,res2b,a))
  res1a=T11[3,1]
  
  R(y,a,b)=(1+(1/2)*y*a*b)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1a,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VA = polrootsreal(subst(res1b,y,VY[iy]),[0,1]);
           if (#VA > 0,

              for(ia=1,#VA,
                    VB = polrootsreal(subst(subst(Ly,y,VY[iy]),a,VA[ia]),[0,1]);

                    if(#VB>0,
                       for(ib=1,#VB,

                           tot+=1;
                           value=R(VY[iy],VA[ia],VB[ib]);
                           if(value>thresold,
                               print("[y,a,b]=",[VY[iy],VA[ia],VB[ib]]);
                               print("R(y,a,b)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* 1 stationary point with value <= 4.6 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.23) */

  /* a = 1: (1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) */
    L=(1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b)
    
    Ty=factor(deriv(L,y))
    Tb=factor(deriv(L,b))
    Ly=Ty[5,1]
    Lb=Tb[4,1]

    T1=factor(polresultant(Ly,Lb,b))
    res1b=T1[1,1]

    u=polrootsreal(res1b,[0,1]) /* length 1 */
    y0=u[1]
    v=polrootsreal(subst(Lb,y,y0),[0,1]) /* length 1 */
    b0=u[1]
    subst(subst(L,y,y0),b,b0)
    \\ 4.2447729912516686734550220438570602195  <----
    

    /* y = 0: (1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) <= 4    */ /* PARI my(y=0);ploth(b=0,1,(1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b)) */
    /* b = 0: (1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) <= 4.75 */ /* PARI my(b=0);ploth(y=0,1,(1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b)) */
    /* y = 1: NULL */
    /* b = 1: NULL */

  /* b = 1: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 */
    L=(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2
    
    Ty=factor(deriv(L,y))
    Ta=factor(deriv(L,a))
    Ly=Ty[5,1]
    La=Ta[4,1]

    T1=factor(polresultant(Ly,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* y = 0: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 2.4  */ /* PARI my(y=0);ploth(a=0,1,(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(y=0,1,(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* y = 1: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.23 */ /* PARI my(y=1);ploth(a=0,1,(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */
.........................................................................

  /* z = 1:
  BORDER

  a = 1: (1+(1/2)*y*t*b)*(1+y)*(1+y*t)*(1-y*t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) */

  L= (1+(1/2)*y*t*b)*(1+y)*(1+y*t)*(1-y*t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Tb=factor(deriv(L,b))
  Ly=Ty[6,1]
  Lt=Tt[4,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(Lt,Ly,b))
  T2=factor(polresultant(Lb,Ly,b))
  res1b=T1[4,1]
  res2b=T2[2,1]

  T11=factor(polresultant(res1b,res2b,t))
  res1t=T11[3,1]*T11[4,1]
  
  R(y,t,b)=(1+(1/2)*y*t*b)*(1+y)*(1+y*t)*(1-y*t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1t,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VT = polrootsreal(subst(res1b,y,VY[iy]),[0,1]);
           if (#VT > 0,

              for(it=1,#VT,
                    VB = polrootsreal(subst(subst(Ly,y,VY[iy]),t,VT[it]),[0,1]);

                    if(#VB>0,
                       for(ib=1,#VB,

                           tot+=1;
                           value=R(VY[iy],VT[it],VB[ib]);
                           if(value>thresold,
                               print("[y,t,b]=",[VY[iy],VT[it],VB[ib]]);
                               print("R(y,t,b)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* 3 stationary point with value <= 4.6 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.23) */

  /* t = 1: (1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) =
          4*(1+(1/2)*y*b)*(1+y)*(1+y)*(1-y)*(1-b)*(1+b+b^2)*(1-b)*(1+b) <=
          4*(1+y)*  (1+y)*(1+y)*(1-y)*(1-b)*(1+b+b^2)*(1-b)*(1+b) <=
          4*(27/16)                  *1 = 6.75 */

  /* b = 1: NULL */
.........................................................................

  /* z = 1:
  BORDER

  b = 1: (1+(1/2)*y*t*a)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */

  L= (1+(1/2)*y*t*a)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Ta=factor(deriv(L,a))
  Ly=Ty[6,1]
  Lt=Tt[4,1]
  La=Ta[3,1]

  T1=factor(polresultant(Lt,Ly,a))
  T2=factor(polresultant(La,Ly,a))
  res1a=T1[5,1]
  res2a=T2[5,1]

  T11=factor(polresultant(res1a,res2a,t))
  res1t=T11[2,1]*T11[3,1]

  R(y,t,a)=(1+(1/2)*y*t*a)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1t,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VT = polrootsreal(subst(res1a,y,VY[iy]),[0,1]);
           if (#VT > 0,

              for(it=1,#VT,
                    VA = polrootsreal(subst(subst(Ly,y,VY[iy]),t,VT[it]),[0,1]);

                    if(#VA>0,
                       for(ia=1,#VA,

                           tot+=1;
                           value=R(VY[iy],VT[it],VA[ia]);
                           if(value>thresold,
                               print("[y,t,a]=",[VY[iy],VT[it],VA[ia]]);
                               print("R(y,t,a)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* No stationary points */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.23) */

  /* t = 1: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 */

    L=(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2
    Ty=factor(deriv(L,y))
    Ta=factor(deriv(L,a))
    Ly=Ty[5,1]
    La=Ta[4,1]

    T1=factor(polresultant(Ly,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* t = 0: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 2.4  */ /* PARI my(y=0);ploth(a=0,1,(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(y=0,1,(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* t = 1: (1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.23 */ /* PARI my(y=1);ploth(a=0,1,(1+(1/2)*y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */

  /* a = 1: NULL */
.........................................................................

  /* BORDER
     t = 1: (1+(1/2)*y*z*a*b)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER
  */   

  R(y,z,a,b) = (1+(1/2)*y*z*a*b)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-B_Pari/auxiliary_23G-B_Pari_Q4.gp")
  
  search(thresold)={
  tot=0;
  VY=polrootsreal(res1z,[0,1]);
  print(#VY);
  if (#VY > 0,

      for(iy=1,#VY,
         VZ = polrootsreal(subst(res1a,y,VY[iy]),[0,1]);
         if (#VZ > 0,

            for(iz=1,#VZ,
                  VA = polrootsreal(subst(subst(res1b,y,VY[iy]),z,VZ[iz]),[0,1]);

                  if(#VA>0,
                     for(ia=1,#VA,
                           VB=polrootsreal(subst(subst(subst(Ly,y,VY[iy]),z,VZ[iz]),a,VA[ia]),[0,1]);

                         if(#VB>0,
                            for(ib=1,#VB,
                                tot+=1;
                                value=R(VY[iy],VZ[iz],VA[ia],VB[ib]);
                                if(value>thresold,
                                    print("[y,z,a,b]=",[VY[iy],VZ[iz],VA[ia],VB[ib]]);
                                    print("R(y,z,a,b)=",value)
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

  search(4);
  /* no stationary points */
.........................................................................

  /* t = 1:
  BORDER

  a = 1: (1+(1/2)*y*z*b)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b) */

  L= (1+(1/2)*y*z*b)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b)
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tb=factor(deriv(L,b))
  Ly=Ty[6,1]
  Lz=Tz[5,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(Lz,Ly,b))
  T2=factor(polresultant(Lb,Ly,b))
  res1b=T1[4,1]
  res2b=T2[2,1]

  T11=factor(polresultant(res1b,res2b,z))
  res1z=T11[3,1]*T11[4,1]
  
  R(y,z,b)=(1+(1/2)*y*z*b)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b)

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1z,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VZ = polrootsreal(subst(res1b,y,VY[iy]),[0,1]);
           if (#VZ > 0,

              for(iz=1,#VZ,
                    VB = polrootsreal(subst(subst(Ly,y,VY[iy]),z,VZ[iz]),[0,1]);

                    if(#VB>0,
                       for(ib=1,#VB,

                           tot+=1;
                           value=R(VY[iy],VZ[iz],VB[ib]);
                           if(value>thresold,
                               print("[y,z,b]=",[VY[iy],VZ[iz],VB[ib]]);
                               print("R(y,z,b)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* 2 stationary point with value <= 4.8 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.23) */

  /* z = 1: already discussed in greater generality (it is <= 8.23) */

  /* b = 1: NULL */
.........................................................................

  /* t = 1:
  BORDER

  b = 1: (1+(1/2)*y*z*a)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 */

  L= (1+(1/2)*y*z*a)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Ta=factor(deriv(L,a))
  Ly=Ty[6,1]
  Lz=Tz[5,1]
  La=Ta[3,1]

  T1=factor(polresultant(Lz,Ly,a))
  T2=factor(polresultant(La,Ly,a))
  res1a=T1[5,1]
  res2a=T2[5,1]

  T11=factor(polresultant(res1a,res2a,z))
  res1t=T11[2,1]*T11[3,1]

  R(y,z,a)=(1+(1/2)*y*z*a)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2

  search(thresold)={
    tot=0;
    VY=polrootsreal(res1t,[0,1]);
    print(#VY);
    if (#VY > 0,

        for(iy=1,#VY,
           VZ = polrootsreal(subst(res1a,y,VY[iy]),[0,1]);
           if (#VZ > 0,

              for(iz=1,#VZ,
                    VA = polrootsreal(subst(subst(Ly,y,VY[iy]),z,VZ[iz]),[0,1]);

                    if(#VA>0,
                       for(ia=1,#VA,

                           tot+=1;
                           value=R(VY[iy],VZ[iz],VA[ia]);
                           if(value>thresold,
                               print("[y,z,a]=",[VY[iy],VZ[iz],VA[ia]]);
                               print("R(y,z,a)=",value)
                             )
                           );
                       );
                  );
               );
            );
         );
  print("number of stationary points: ",tot)
  };

  search(4);
  /* No stationary points */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.23) */

  /* z = 1: already discussed in greater generality (it is <= 8.23) */

  /* a = 1: NULL */
.........................................................................

/* BORDER
    a = 1: (1+(1/2)*y*z*t*b)*(1+y)*(1+y*z*t)*(1-y*z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
    
    CASE 1: reduction
*/

{res1b=
y^4*z^2*t^2 - 5/2*y^3*z^3*t^3 + 3/4*y^3*z^2*t^2 - 19/4*y^2*z^3*t^3 - 3/2*y^2*z^2*t^2 - 1/2*y^2 -
    2*y*z^3*t^3 - y*z^2*t^2 + y*z*t - 1/4*y + 5/4*z*t + 1/2
;}

{res2b=
y^10*z^8*t^6 - 5/2*y^9*z^9*t^7 + 3/2*y^9*z^8*t^7 + 9/4*y^9*z^8*t^6 - 17/2*y^8*z^9*t^7 +
    43/4*y^8*z^8*t^7 + 3/16*y^8*z^8*t^6 - 9/2*y^8*z^7*t^7 + 15/4*y^8*z^7*t^6 - 9/4*y^8*z^6*t^6 -
    3/2*y^8*z^6*t^4 - 337/32*y^7*z^9*t^7 + 799/32*y^7*z^8*t^7 - 181/64*y^7*z^8*t^6 -
    105/4*y^7*z^7*t^7 + 199/16*y^7*z^7*t^6 + 7/2*y^7*z^7*t^5 - 255/16*y^7*z^6*t^6 - 2*y^7*z^6*t^5 -
    3*y^7*z^6*t^4 - 363/64*y^6*z^9*t^7 + 413/16*y^6*z^8*t^7 - 75/32*y^6*z^8*t^6 - 869/16*y^6*z^7*t^7
    + 479/32*y^6*z^7*t^6 + 85/8*y^6*z^7*t^5 - 145/4*y^6*z^6*t^6 - 103/8*y^6*z^6*t^5 +
    1/32*y^6*z^6*t^4 + 21/4*y^6*z^5*t^5 - 37/8*y^6*z^5*t^4 + 21/8*y^6*z^4*t^4 + 3/4*y^6*z^4*t^2 -
    9/8*y^5*z^9*t^7 + 99/8*y^5*z^8*t^7 - 9/16*y^5*z^8*t^6 - 207/4*y^5*z^7*t^7 + 125/16*y^5*z^7*t^6 +
    181/16*y^5*z^7*t^5 - 73/2*y^5*z^6*t^6 - 819/32*y^5*z^6*t^5 + 205/64*y^5*z^6*t^4 +
    211/8*y^5*z^5*t^5 - 209/16*y^5*z^5*t^4 - 13/8*y^5*z^5*t^3 + 259/16*y^5*z^4*t^4 + 7/8*y^5*z^4*t^3
    + 21/16*y^5*z^4*t^2 + 9/4*y^4*z^8*t^7 - 93/4*y^4*z^7*t^7 + 3/2*y^4*z^7*t^6 + 319/64*y^4*z^7*t^5
    - 17*y^4*z^6*t^6 - 169/8*y^4*z^6*t^5 + 67/32*y^4*z^6*t^4 + 697/16*y^4*z^5*t^5 -
    399/32*y^4*z^5*t^4 - 69/16*y^4*z^5*t^3 + 59/2*y^4*z^4*t^4 + 5*y^4*z^4*t^3 - 1/8*y^4*z^4*t^2 -
    2*y^4*z^3*t^3 + 15/8*y^4*z^3*t^2 - y^4*z^2*t^2 - 1/8*y^4*z^2 - 4*y^3*z^7*t^7 + 3/4*y^3*z^7*t^5 -
    3*y^3*z^6*t^6 - 29/4*y^3*z^6*t^5 + 3/8*y^3*z^6*t^4 + 29*y^3*z^5*t^5 - 37/8*y^3*z^5*t^4 -
    121/32*y^3*z^5*t^3 + 83/4*y^3*z^4*t^4 + 261/32*y^3*z^4*t^3 - 71/64*y^3*z^4*t^2 -
    17/2*y^3*z^3*t^3 + 69/16*y^3*z^3*t^2 + 1/4*y^3*z^3*t - 85/16*y^3*z^2*t^2 - 1/8*y^3*z^2*t -
    3/16*y^3*z^2 - 3/4*y^2*z^6*t^5 + 27/4*y^2*z^5*t^5 - 1/2*y^2*z^5*t^4 - 81/64*y^2*z^5*t^3 +
    5*y^2*z^4*t^4 + 77/16*y^2*z^4*t^3 - 17/32*y^2*z^4*t^2 - 167/16*y^2*z^3*t^3 + 93/32*y^2*z^3*t^2 +
    9/16*y^2*z^3*t - 29/4*y^2*z^2*t^2 - 5/8*y^2*z^2*t + 1/32*y^2*z^2 + 1/4*y^2*z*t - 1/4*y^2*z +
    1/8*y^2 - 1/8*y*z^5*t^3 + 7/8*y*z^4*t^3 - 1/16*y*z^4*t^2 - 15/4*y*z^3*t^3 + 9/16*y*z^3*t^2 +
    3/8*y*z^3*t - 11/4*y*z^2*t^2 - 25/32*y*z^2*t + 7/64*y*z^2 + 7/8*y*z*t - 7/16*y*z + 9/16*y +
    5/64*z^3*t - 1/4*z^2*t + 1/32*z^2 + 11/16*z*t - 5/32*z + 1/2
;}

{res3b=
y^16*z^14*t^12 + 19/4*y^15*z^14*t^12 - y^15*z^13*t^12 + 73/8*y^14*z^14*t^12 - 17/4*y^14*z^13*t^12 +
    3/4*y^14*z^12*t^12 - 13/4*y^14*z^12*t^10 + 291/32*y^13*z^14*t^12 - 113/16*y^13*z^13*t^12 +
    35/16*y^13*z^12*t^12 - 207/16*y^13*z^12*t^10 + 2*y^13*z^11*t^10 + 1269/256*y^12*z^14*t^12 -
    369/64*y^12*z^13*t^12 + 81/64*y^12*z^12*t^12 - 1245/64*y^12*z^12*t^10 + 61/16*y^12*z^11*t^10 +
    17/16*y^12*z^10*t^10 + 25/8*y^12*z^10*t^8 + 1431/1024*y^11*z^14*t^12 - 297/128*y^11*z^13*t^12 -
    519/256*y^11*z^12*t^12 - 3401/256*y^11*z^12*t^10 - 285/64*y^11*z^11*t^10 + 989/64*y^11*z^10*t^10
    + 29/4*y^11*z^10*t^8 + 3/4*y^11*z^9*t^8 + 81/512*y^10*z^14*t^12 - 189/512*y^10*z^13*t^12 -
    203/64*y^10*z^12*t^12 - 1733/512*y^10*z^12*t^10 - 549/32*y^10*z^11*t^10 +
    6419/128*y^10*z^10*t^10 - 131/128*y^10*z^10*t^8 + 517/32*y^10*z^9*t^8 - 257/32*y^10*z^8*t^8 -
    13/16*y^10*z^8*t^6 - 51/32*y^9*z^12*t^12 + 363/1024*y^9*z^12*t^10 - 4361/256*y^9*z^11*t^10 +
    18149/256*y^9*z^10*t^10 - 4127/256*y^9*z^10*t^8 + 3445/64*y^9*z^9*t^8 - 899/16*y^9*z^8*t^8 +
    237/64*y^9*z^8*t^6 - 13/4*y^9*z^7*t^6 - 9/32*y^8*z^12*t^12 + 117/512*y^8*z^12*t^10 -
    3737/512*y^8*z^11*t^10 + 405/8*y^8*z^10*t^10 - 275/16*y^8*z^10*t^8 + 9517/128*y^8*z^9*t^8 -
    273/2*y^8*z^8*t^8 + 1253/64*y^8*z^8*t^6 - 199/8*y^8*z^7*t^6 + 309/32*y^8*z^6*t^6 -
    11/32*y^8*z^6*t^4 - 75/64*y^7*z^11*t^10 + 1151/64*y^7*z^10*t^10 - 3749/512*y^7*z^10*t^8 +
    1621/32*y^7*z^9*t^8 - 20135/128*y^7*z^8*t^8 + 467/16*y^7*z^8*t^6 - 3827/64*y^7*z^7*t^6 +
    3507/64*y^7*z^6*t^6 - 155/32*y^7*z^6*t^4 + 9/4*y^7*z^5*t^4 + 81/32*y^6*z^10*t^10 -
    291/256*y^6*z^10*t^8 + 4311/256*y^6*z^9*t^8 - 2987/32*y^6*z^8*t^8 + 4993/256*y^6*z^8*t^6 -
    4129/64*y^6*z^7*t^6 + 6975/64*y^6*z^6*t^6 - 837/64*y^6*z^6*t^4 + 53/4*y^6*z^5*t^4 -
    77/16*y^6*z^4*t^4 + 7/32*y^6*z^4*t^2 + 139/64*y^5*z^9*t^8 - 1751/64*y^5*z^8*t^8 +
    3075/512*y^5*z^8*t^6 - 4347/128*y^5*z^7*t^6 + 12633/128*y^5*z^6*t^6 - 455/32*y^5*z^6*t^4 +
    1637/64*y^5*z^5*t^4 - 1465/64*y^5*z^4*t^4 + 51/32*y^5*z^4*t^2 - 5/8*y^5*z^3*t^2 -
    99/32*y^4*z^8*t^8 + 173/256*y^4*z^8*t^6 - 2089/256*y^4*z^7*t^6 + 1341/32*y^4*z^6*t^6 -
    1833/256*y^4*z^6*t^4 + 1365/64*y^4*z^5*t^4 - 2319/64*y^4*z^4*t^4 + 13/4*y^4*z^4*t^2 -
    97/32*y^4*z^3*t^2 + 35/32*y^4*z^2*t^2 - 1/32*y^4*z^2 - 41/64*y^3*z^7*t^6 + 445/64*y^3*z^6*t^6 -
    1629/1024*y^3*z^6*t^4 + 1009/128*y^3*z^5*t^4 - 5963/256*y^3*z^4*t^4 + 689/256*y^3*z^4*t^2 -
    37/8*y^3*z^3*t^2 + 69/16*y^3*z^2*t^2 - 11/64*y^3*z^2 + 1/16*y^3*z + 3/32*y^2*z^6*t^6 -
    59/512*y^2*z^6*t^4 + 543/512*y^2*z^5*t^4 - 347/64*y^2*z^4*t^4 + 491/512*y^2*z^4*t^2 -
    173/64*y^2*z^3*t^2 + 647/128*y^2*z^2*t^2 - 35/128*y^2*z^2 + 1/4*y^2*z - 3/32*y^2 +
    1/64*y*z^5*t^4 - 7/64*y*z^4*t^4 + 127/1024*y*z^4*t^2 - 137/256*y*z^3*t^2 + 457/256*y*z^2*t^2 -
    41/256*y*z^2 + 9/32*y*z - 19/64*y + 1/512*z^4*t^2 - 5/512*z^3*t^2 + 1/32*z^2*t^2 - 1/32*z^2 +
    11/128*z - 7/32
;}

{res1ta = y^2 + 5/4*y - 3/4;}


res1bSy = lift(Mod(res1b,res1ta))
\\ res1bSy := (-59/32*t^3*z^3 - 39/32*t^2*z^2 + t*z + 3/8)*y + (-39/32*t^3*z^3 - 3/32*t^2*z^2 + 5/4*t*z + 1/8);
res2bSy = lift(Mod(res2b,res1ta))
\\ res2bSy := (-214955/32768*t^7*z^9 + (-189299/32768*t^7 - 221271/32768*t^6)*z^8 + (-224715/32768*t^7 - 200239/32768*t^6 + 19011/4096*t^5)*z^7 + (-219991/32768*t^6 + 1403/512*t^5 + 39241/8192*t^4)*z^6 + (24401/4096*t^5 + 28427/8192*t^4 - 2811/2048*t^3)*z^5 + (44325/8192*t^4 + 583/2048*t^3 - 2445/2048*t^2)*z^4 + (-173/64*t^3 - 141/512*t^2 + 1/4*t)*z^3 + (-549/256*t^2 - 37/128*t + 59/512)*z^2 + (9/16*t -1/8)*z + 13/32)*y + (91569/32768*t^7*z^9 + (103305/32768*t^7 + 96885/32768*t^6)*z^8 + (13329/32768*t^7 + 100221/32768*t^6 - 6201/4096*t^5)*z^7 + (36213/32768*t^6 - 1713/512*t^5 - 16107/8192*t^4)*z^6 + (13821/4096*t^5 - 22737/8192*t^4 - 447/2048*t^3)*z^5 + (14913/8192*t^4 + 4347/2048*t^3 + 663/2048*t^2)*z^4 + (-213/64*t^3 + 711/512*t^2 + 17/64*t)*z^3 + (-561/256*t^2 - 77/128*t + 7/512)*z^2 + (7/8*t - 11/32)*z + 19/32);
res3bSy = lift(Mod(res3b,res1ta))
\\ res3bSy := (-969625549/33554432*t^12*z^14 - 905307061/33554432*t^12*z^13 + (-844342797/33554432*t^12 + 2127689929/33554432*t^10)*z^12 + 1936884009/33554432*t^10*z^11 + (1762958361/33554432*t^10 - 422388523/8388608*t^8)*z^10 - 368019869/8388608*t^8*z^9 + (-330316115/8388608*t^8 + 40880943/2097152*t^6)*z^8 + 2009625/131072*t^6*z^7 + (8186429/524288*t^6 - 2225741/524288*t^4)*z^6 - 16345/8192*t^4*z^5 + (-19199/4096*t^4 + 22529/32768*t^2)*z^4 - 165/512*t^2*z^3 + (2565/2048*t^2 - 197/2048)*z^2 + 29/256*z- 23/128)*y + (429573687/33554432*t^12*z^14 + 401018991/33554432*t^12*z^13 + (373959927/33554432*t^12 - 942688491/33554432*t^10)*z^12 - 859088907/33554432*t^10*z^11 + (-773687259/33554432*t^10 + 186544689/8388608*t^8)*z^10 + 167096487/8388608*t^8*z^9 + (129050409/8388608*t^8 - 17261181/2097152*t^6)*z^8 - 1109979/131072*t^6*z^7 + (-1035399/524288*t^6 + 569271/524288*t^4)*z^6 + 25803/8192*t^4*z^5 +(-11811/4096*t^4 + 10765/32768*t^2)*z^4 - 149/128*t^2*z^3 + (3433/2048*t^2 - 265/2048)*z^2 + 55/256*z - 37/128);


{res1tSy=
z^36 + 3004665437910806708625953450590837755792427061711086707233/276107653310474645188402374584868\
    130068878034880232294400*z^35 + 541553364684672719751492715522248868498675577306180523000279274\
    18429/1466233264365928432490223233177248975569168036739063641083206041600*z^34 - 
    75515661786263228983891172301864552553078987529726713480368767805346199/13460021366879223010260\
    24928056714559572496257726460422514383146188800*z^33 - 
    4679711931727938291265662677219454634257767978726402972192128910503245799/203395878432841592155\
    04376690634797789095499005644290829106234209075200*z^32 + 
    6258992022395254180170307538709732127663682951753557887860887573965195673/207463796001498423998\
    1446422444749374487740898575717664568835889325670400*z^31 - 
    6726456159789171561413373095240439374704664347113492865232934159025506988987/627001694582306348\
    0832815854499686998451839160139946719585815132184248320*z^30 - 
    32683037212364411293199289470215958995254538824761698110903538834614134159945441/14389688890663\
    930688511312386076781661446970872521177721449445728362849894400*z^29 - 
    139523640962144504039302486615479378648601624816854157337299483770729666846912221/9784988445651\
    4728681876924225322115297839401933144008505856230952867379281920*z^28 - 
    2882104058266663718270589662591004840766314263407061085850656690436041717374189/402343274903432\
    272540612352900173171454931751369835561290527265431198105600*z^27 - 
    236114443974648518584658347869139784618234756794309646409048413962704049985785631/4119995135011\
    1464708158704936977732756985011340271161476149991980154686013440*z^26 - 
    5541226823726070712234885790335257312004695892364833243180557908012396314868071193/434888375362\
    287683030564107668098290212619564147306704470472137568299463475200*z^25 - 
    4115313304968061980702634173585419309868113415286281609765593058492211221504013501/244624711141\
    286821704692310563305288244598504832860021264640577382168448204800*z^24 - 
    11529258334415760804682229782460603010377031871732680671206525635874333217793484183/78279907565\
    2117829455015393802576922382715215465152068046849847622939034255360*z^23 - 
    21200292619117612578321323788075941391742304583805438267551430741832160474590022483/48924942228\
    2573643409384621126610576489197009665720042529281154764336896409600*z^22 - 
    1204113242445681441826958327580591883898507118257438289131913326125113129095784073/100358855852\
    835619160899409461868836202912207110916931800878185592684491571200*z^21 - 
    109945426131974486748261351493052882437698847823682373435960011935034687377325068833/1304665126\
    086863049091692323004294870637858692441920113411416412704898390425600*z^20 - 
    17059358264138269192844029630121213109105867157935043264959940793830409848063345139/19569976891\
    30294573637538484506442305956788038662880170117124619057347585638400*z^19 - 
    23867718456185322524181623107944400599654914243391926836054009314835118367768232111/16308314076\
    0857881136461540375536858829732336555240014176427051588112298803200*z^18 - 
    790200910022480017206019587014231809435662942557318294424531856764568824521513/1005858187258580\
    68135153088225043292863732937842458890322631816357799526400*z^17 - 
    235857586298032687561943864220373929956642387897818723559593679602396967561910931/1132521810839\
    290841225427363719005964095363448300277876225187858250779852800*z^16 - 
    49594311428544251182910540572087948108913241142262602552786756103426549880958757/15289044446330\
    426356543269410206580515287406552053751329040036086385528012800*z^15 - 
    1278837239875105550679598886413462076778638097573172739785800360679925263721850487/509634814877\
    6808785514423136735526838429135517351250443013345362128509337600*z^14 + 
    171868680813915594624262209327561362564132427596590826330146165772939581509524647/7644522223165\
    213178271634705103290257643703276026875664520018043192764006400*z^13 - 
    34383830750451204201142235379063787205165442148577717737632594463010193066117511/15926087964927\
    5274547325723022985213700910484917226576344167042566515916800*z^12 + 
    2276575930883786603639057951972262204594416308191626920169901273462132610043041/353913065887278\
    38788294605116218936377980107759383683632037120570336870400*z^11 - 
    28021822228518325558046359046134674347664508896646643166327844589312005931460683/19111305557913\
    0329456790867627582256441092581900671891613000451079819100160*z^10 + 
    16086587778500474512403540189428580806645431775861110164726463233215467691263829/15926087964927\
    5274547325723022985213700910484917226576344167042566515916800*z^9 - 
    7895756203795952262618150554673299850831811154094159916252022725618822900690099/119445659736956\
    455910494292267238910275682863687919932258125281924886937600*z^8 + 
    93515898796338490919388928803305162852346640717225865175638055083581260130619/11485159590091966\
    91447060502569604906496950612383845502481973864662374400*z^7 - 
    43988836623376990762596535186828499369433361435171178141360538729360476256191/93316921669497231\
    1800736658337803986528772372561874470766603765038179200*z^6 + 
    86875849837320401841715128145723534251732270398913453736405199519909128881863/18663384333899446\
    23601473316675607973057544745123748941533207530076358400*z^5 - 
    2801378313310867629475684041130810879009769576312701690150709960261039703163/933169216694972311\
    80073665833780398652877237256187447076660376503817920*z^4 + 
    177607211012991420621128694602686145548750822366855230301658656001251451707/2916153802171788474\
    3773020573056374579024136642558577211456367657443100*z^3 - 
    897389274808651526482796299547989500997217886990843110024192313440710366/4764957193091157637871\
    4085903686886567032903010716629430484260878175*z^2 + 
    1882712378009561604526440101470015909743761482515284374276674132313804172/270014240941832266146\
    046486787559023879853117060727566772744144976325*z + 
    623119996844017687425897507510805745100837372221260487005412794171427248/2700142409418322661460\
    46486787559023879853117060727566772744144976325
;}

{res1ySy=
z^10*t^8 + 561569/574409*z^9*t^8 - 432/661*z^9*t^7 + 397929/574409*z^8*t^8 - 359840/574409*z^8*t^7 -
    1463827/1723227*z^8*t^6 - 299600/574409*z^7*t^7 - 1567987/1723227*z^7*t^6 +
    775924/1723227*z^7*t^5 - 150081/574409*z^6*t^6 + 770612/1723227*z^6*t^5 + 37132/156657*z^6*t^4 +
    42140/156657*z^5*t^5 + 58660/156657*z^5*t^4 - 4512/52219*z^5*t^3 - 255548/1723227*z^4*t^4 -
    204304/1723227*z^4*t^3 - 41168/1723227*z^4*t^2 - 58240/1723227*z^3*t^3 - 148928/1723227*z^3*t^2
    + 5696/1723227*z^3*t + 54144/574409*z^2*t^2 + 1792/156657*z^2*t + 1216/1723227*z^2 +
    14848/1723227*z - 2048/156657
;}

/* we compute y from the modulo: */
w=polrootsreal(res1ta,[0,1]) /* length 1 */
y0=w[1]

/* we compute z from res1tSy: */
u=polrootsreal(res1tSy,[0,1]) /* length 1 */
z0=u[1]

/* we compute t from res1ySy: */
v=polrootsreal(subst(res1ySy,z,z0),[0,1]) /* length 2 */
t0=v[1]
/* Now we know the values for y0, z0, and t0. res1bSy depends only on y,z,t, hence in case y0,z0,t0 
is a stationary point res1bSy(y0,z0,t0) shoud be 0: substituting in res1bSy we get */
subst(subst(subst(res1bSy,y,y0),z,z0),t,t0)
\\ 0.72195309681775748322297413826554140170
/* it is not zero, so that the values we have foud are not consistent and the point is not a stationary point */

t0=v[2]
/* Now we know the values for y0, z0, and t0. res1bSy depends only on y,z,t, hence in case y0,z0,t0 
is a stationary point res1bSy(y0,z0,t0) shoud be 0: substituting in res1bSy we get */
subst(subst(subst(res1bSy,y,y0),z,z0),t,t0)
\\ 0.68532821524638713208588123229609351459
/* it is not zero, so that the values we have foud are not consistent and the point is not a stationary point */
.........................................................................

  /* BORDER
     a = 1: (1+(1/2)*y*z*t*b)*(1+y)*(1+y*z*t)*(1-y*z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
     INNER
     CASE 2
  */   

  R(y,z,t,b) = (1+(1/2)*y*z*t*b)*(1+y)*(1+y*z*t)*(1-y*z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-B_Pari/auxiliary_23G-B_Pari_Q5.gp")
  
  search(thresold)={
  tot=0;
  VZ=polrootsreal(res1y,[0,1]);
  print(#VZ);
  if (#VZ > 0,

      for(iz=1,#VZ,
         VY = polrootsreal(subst(res1tb,z,VZ[iz]),[0,1]);
         if (#VY > 0,

            for(iy=1,#VY,
                  VT = polrootsreal(subst(subst(res1b,y,VY[iy]),z,VZ[iz]),[0,1]);

                  if(#VT>0,
                     for(it=1,#VT,
                           VB=polrootsreal(subst(subst(subst(Ly,y,VY[iy]),z,VZ[iz]),t,VT[it]),[0,1]);

                         if(#VB>0,
                            for(ib=1,#VB,
                                tot+=1;
                                value=R(VY[iy],VZ[iz],VT[it],VB[ib]);
                                if(value>thresold,
                                    print("[y,z,t,b]=",[VY[iy],VZ[iz],VT[it],VB[ib]]);
                                    print("R(y,z,t,b)=",value)
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

  search(4);
  /* no stationary points */
.........................................................................
