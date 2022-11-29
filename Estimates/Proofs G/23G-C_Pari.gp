
/* IMPORTANT:

   Some computations in this file require data and functions which are stored in separated 
   files located into the directory "/Estimates/Proofs G/auxiliary_files_23G-C_Pari".
   To import these data we use the Pari instruction 
   read("ROOT/Estimates/Proofs G/auxiliary_files_23G-C_Pari/name of the file needed.gp")
   In order to make Pari execute correctly this command you have to substitute ROOT with 
   the path of the directory "/Estimates" in your PC.
*/

/*
     INNER 
*/

  R(y,z,t,a,b) = (1+y*z*t*a*b)*(1+y)*(1+y*z*t)*(1-y*z*t*a)*(1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-C_Pari/auxiliary_23G-C_Pari_Q1.gp")
  
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
                             VA=polrootsreal(subst(subst(subst(res2b,y,VY[iy]),z,VZ[iz]),t,VT[it]),[0,1]);
                             
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
  
  search(4);
  /* 1 stationary point with value <= 4.5 */
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
     y = 1: (1+z*t*a*b)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER
  */   

  R(z,t,a,b) = (1+z*t*a*b)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-C_Pari/auxiliary_23G-C_Pari_Q2.gp")
  
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
  /* 4 stationary points with value <= 4.5 */
.........................................................................

  /* y = 1:
  BORDER

  z = 1: (1+t*a*b)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b) */

  L= (1+t*a*b)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Tt=factor(deriv(L,t))
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  Lt=Tt[4,1]
  La=Ta[3,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(La,Lt,b))
  T2=factor(polresultant(Lb,Lt,b))
  res1b=T1[6,1]
  res2b=T2[6,1]

  T11=factor(polresultant(res1b,res2b,t))
  res1t=T11[4,1]*T11[5,1]

  R(t,a,b)=(1+t*a*b)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)

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
  /* 1 stationary point with value <= 4.9 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* t = 1: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */
    L=(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
    Ta=factor(deriv(L,a))
    Tb=factor(deriv(L,b))
    La=Ta[4,1]
    Lb=Tb[4,1]

    T1=factor(polresultant(La,Lb,b))
    res1b=T1[4,1]

    u=polrootsreal(res1b,[0,1]) /* no roots */

    /* a = 0: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8    */ /* PARI my(a=0);ploth(b=0,1,(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* b = 0: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 4.75 */ /* PARI my(b=0);ploth(a=0,1,(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* a = 1: NULL */
    /* b = 1: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8.79 */ /* PARI my(b=1);ploth(a=0,1,(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */

  /* a = 1: (1+t*b)*2*(1+t)*(1-t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) =
          4*(1+t*b)*(1-t*b)     *(1+t)*(1-t)*(1+t)*        (1+t*b+(t*b)^2)*(1-b)*(1+b) <=
          4*1                   *(32/27)          *        (1+b+b^2)*(1-b)*(1+b) <=
          4*1                   *(32/27)          *        (27/16) <= 8 */

  /* b = 1: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */
    L=(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
    
    Tt=factor(deriv(L,t))
    Ta=factor(deriv(L,a))
    Lt=Tt[5,1]
    La=Ta[4,1]

    T1=factor(polresultant(Lt,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* t = 0: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 4    */ /* PARI my(t=0);ploth(a=0,1,(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(t=0,1,(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* t = 1: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8.79 */ /* PARI my(t=1);ploth(a=0,1,(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */
.........................................................................

  /* y = 1:
  BORDER

  t = 1: (1+z*a*b)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */

  L= (1+z*a*b)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Tz=factor(deriv(L,z))
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  Lz=Tz[5,1]
  La=Ta[3,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(Lz,Lb,b))
  T2=factor(polresultant(La,Lb,b))
  res1b=T1[6,1]
  res2b=T2[10,1]

  T11=factor(polresultant(res1b,res2b,a))
  res1a=T11[7,1]*T11[8,1]

  R(z,a,b)=(1+z*a*b)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)

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
  /* 1 stationary points with value <= 6.2 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* z = 1: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */
    L=(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
    Ta=factor(deriv(L,a))
    Tb=factor(deriv(L,b))
    La=Ta[4,1]
    Lb=Tb[4,1]

    T1=factor(polresultant(La,Lb,b))
    res1b=T1[4,1]

    u=polrootsreal(res1b,[0,1]) /* no roots */

    /* a = 0: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8    */ /* PARI my(a=0);ploth(b=0,1,(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* b = 0: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 4.75 */ /* PARI my(b=0);ploth(a=0,1,(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */
    /* a = 1: NULL */
    /* b = 1: (1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) <= 8.79 */ /* PARI my(b=1);ploth(a=0,1,(1+a*b)*2*2*(1-a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)) */

  /* a = 1: (1+z*b)*2*(1+z)*(1-z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b) <= 
          4*(1+z*b)*(1-z*b)      *(1+z)*(1-z)*(1+z)    *(1+b+b^2)*(1-b)*(1+b) <=
          4*1                          *(32/27)              *(27/16) <= 8 */
  */

  /* b = 1: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 */
    L=(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2
    
    Tz=factor(deriv(L,z))
    Ta=factor(deriv(L,a))
    Lz=Tz[6,1]
    La=Ta[4,1]

    T1=factor(polresultant(Lz,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* z = 0: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 5.3  */ /* PARI my(z=0);ploth(a=0,1,(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(z=0,1,(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* z = 1: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.79 */ /* PARI my(z=1);ploth(a=0,1,(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */
.........................................................................

  /* y = 1:
  BORDER

  a = 1: (1+z*t*b)*2*(1+z*t)*(1-z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) */

  L= (1+z*t*b)*2*(1+z*t)*(1-z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
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

  R(z,t,b)=(1+z*t*b)*2*(1+z*t)*(1-z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)

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

  /* z = 1: (1+t*b)*2*(1+t)*(1-t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) =
          4*(1+t)*(1-t)*(1+t)*(1+t*b)*(1-t*b)*(1+t*b+(t*b)^2)*(1-b)*(1+b) <=
          4*(32/16)          *(27/16)                        *1 <= 8 */

  /* t = 1: (1+z*b)*2*(1+z)*(1-z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b) =
          4*(1+z*b)*(1-z*b)*(1+z)*(1-z)*(1+z)*(1+b+b^2)*(1-b)*(1+b) <=
          4*1                    *(32/27)    *(27/16) <= 8 */ 

  /* b = 1: NULL */
.........................................................................

  /* y = 1:
  BORDER

  b = 1: (1+z*t*a)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */

  L= (1+z*t*a)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
  Tz=factor(deriv(L,z))
  Tt=factor(deriv(L,t))
  Ta=factor(deriv(L,a))
  Lz=Tz[7,1]
  Lt=Tt[5,1]
  La=Ta[4,1]

  T1=factor(polresultant(Lt,Lz,a))
  T2=factor(polresultant(La,Lz,a))
  res1a=T1[4,1]
  res2a=T2[4,1]

  T11=factor(polresultant(res1a,res2a,t))
  res1t=T11[2,1]*T11[3,1]

  R(z,t,a)=(1+z*t*a)*2*(1+z*t)*(1-z*t*a)*(1+z*t*a)*(1-z*t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2

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

  /* z = 1: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */
    L=(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
    Tt=factor(deriv(L,t))
    Ta=factor(deriv(L,a))
    Lt=Tt[5,1]
    La=Ta[4,1]

    T1=factor(polresultant(Lt,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* t = 0: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 4    */ /* PARI my(t=0);ploth(a=0,1,(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(t=0,1,(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* t = 1: (1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 <= 8.79 */ /* PARI my(t=1);ploth(a=0,1,(1+t*a)*2*(1+t)*(1-t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */

  /* t = 1: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 */
    L=(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2
    Tz=factor(deriv(L,z))
    Ta=factor(deriv(L,a))
    Lz=Tz[6,1]
    La=Ta[4,1]

    T1=factor(polresultant(Lz,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* z = 0: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 5.3  */ /* PARI my(z=0);ploth(a=0,1,(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(z=0,1,(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* z = 1: (1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.79 */ /* PARI my(z=1);ploth(a=0,1,(1+z*a)*2*(1+z)*(1-z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */

  /* a = 1: NULL */
.........................................................................

  /* BORDER
     z = 1: (1+y*t*a*b)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER
  */   

  R(y,t,a,b) = (1+y*t*a*b)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a*b)*(1+t*a*b+(t*a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-C_Pari/auxiliary_23G-C_Pari_Q3.gp")
  
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
  /* no stationary points */
.........................................................................

  /* z = 1:
  BORDER

  t = 1: (1+y*a*b)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b) */

  L= (1+y*a*b)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
  Ty=factor(deriv(L,y))
  Ta=factor(deriv(L,a))
  Tb=factor(deriv(L,b))
  Ly=Ty[6,1]
  La=Ta[5,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(La,Ly,b))
  T2=factor(polresultant(Lb,Ly,b))
  res1b=T1[4,1]
  res2b=T2[4,1]

  T11=factor(polresultant(res1b,res2b,a))
  res1a=T11[2,1]
  
  R(y,a,b)=(1+y*a*b)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)

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
  /* 1 stationary point with value <= 4.8 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.79) */

  /* a = 1: (1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) */
    L=(1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b)
    
    Ty=factor(deriv(L,y))
    Tb=factor(deriv(L,b))
    Ly=Ty[5,1]
    Lb=Tb[4,1]

    T1=factor(polresultant(Ly,Lb,b))
    res1b=T1[3,1]

    u=polrootsreal(res1b,[0,1]) /* length 1 */
    y0=u[1]
    v=polrootsreal(subst(Lb,y,y0),[0,1]) /* length 1 */
    b0=u[1]
    subst(subst(L,y,y0),b,b0)
    \\ 4.3981457012621605547027187596761429002  <----
    

    /* y = 0: (1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) <= 4    */ /* PARI my(y=0);ploth(b=0,1,(1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b)) */
    /* b = 0: (1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) <= 4.75 */ /* PARI my(b=0);ploth(y=0,1,(1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b)) */
    /* y = 1: NULL */
    /* b = 1: NULL */

  /* b = 1: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 */
    L=(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2
    
    Ty=factor(deriv(L,y))
    Ta=factor(deriv(L,a))
    Ly=Ty[5,1]
    La=Ta[4,1]

    T1=factor(polresultant(Ly,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* y = 0: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 2.4  */ /* PARI my(y=0);ploth(a=0,1,(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(y=0,1,(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* y = 1: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.79 */ /* PARI my(y=1);ploth(a=0,1,(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */
.........................................................................

  /* z = 1:
  BORDER

  a = 1: (1+y*t*b)*(1+y)*(1+y*t)*(1-y*t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b) */

  L= (1+y*t*b)*(1+y)*(1+y*t)*(1-y*t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
  Ty=factor(deriv(L,y))
  Tt=factor(deriv(L,t))
  Tb=factor(deriv(L,b))
  Ly=Ty[6,1]
  Lt=Tt[4,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(Lt,Ly,b))
  T2=factor(polresultant(Lb,Ly,b))
  res1b=T1[5,1]
  res2b=T2[5,1]

  T11=factor(polresultant(res1b,res2b,t))
  res1t=T11[3,1]*T11[4,1]
  
  R(y,t,b)=(1+y*t*b)*(1+y)*(1+y*t)*(1-y*t)*(1+t)*(1-t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)

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
  /* 3 stationary point with value <= 4.7 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.79) */

  /* t = 1: (1+y*b)*(1+y)*(1+y)*(1-y)*2*(1-b)*(1+b+b^2)*2*(1-b)*(1+b) =
          4*(1+y*b)*(1+y)*(1+y)*(1-y)*(1-b)*(1+b+b^2)*(1-b)*(1+b) <=
          4*(1+y)*  (1+y)*(1+y)*(1-y)*(1-b)*(1+b+b^2)*(1-b)*(1+b) <=
          4*(27/16)                  *1 = 6.75 */

  /* b = 1: NULL */
.........................................................................

  /* z = 1:
  BORDER

  b = 1: (1+y*t*a)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2 */

  L= (1+y*t*a)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2
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
  res1t=T11[3,1]*T11[4,1]

  R(y,t,a)=(1+y*t*a)*(1+y)*(1+y*t)*(1-y*t*a)*(1+t*a)*(1-t*a)*(1+t*a+(t*a)^2)*(1+a)*(1-a)*2

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

  /* y = 1: already discussed in greater generality (it is <= 8.79) */

  /* t = 1: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 */

    L=(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2
    Ty=factor(deriv(L,y))
    Ta=factor(deriv(L,a))
    Ly=Ty[5,1]
    La=Ta[4,1]

    T1=factor(polresultant(Ly,La,a))
    res1a=T1[3,1]

    u=polrootsreal(res1a,[0,1]) /* no roots */

    /* t = 0: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 2.4  */ /* PARI my(y=0);ploth(a=0,1,(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 0: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8    */ /* PARI my(a=0);ploth(y=0,1,(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* t = 1: (1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2 <= 8.79 */ /* PARI my(y=1);ploth(a=0,1,(1+y*a)*(1+y)*(1+y)*(1-y*a)*(1+a)*(1-a)*(1+a+a^2)*(1+a)*(1-a)*2) */
    /* a = 1: NULL */

  /* a = 1: NULL */
.........................................................................

  /* BORDER
     t = 1: (1+y*z*a*b)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER
  */   

  R(y,z,a,b) = (1+y*z*a*b)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-C_Pari/auxiliary_23G-C_Pari_Q4.gp")
  
  search(thresold)={
  tot=0;
  VY=polrootsreal(res1z,[0,1]);
  print(#VY);
  if (#VY > 0,

      for(iy=1,#VY,
         VZ = polrootsreal(subst(res1a,y,VY[iy]),[0,1]);
         if (#VZ > 0,

            for(iz=1,#VZ,
                  VA = polrootsreal(subst(subst(res2b,y,VY[iy]),z,VZ[iz]),[0,1]);

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

  a = 1: (1+y*z*b)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b) */

  L= (1+y*z*b)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b)
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Tb=factor(deriv(L,b))
  Ly=Ty[6,1]
  Lz=Tz[5,1]
  Lb=Tb[5,1]

  T1=factor(polresultant(Lz,Ly,b))
  T2=factor(polresultant(Lb,Ly,b))
  res1b=T1[5,1]
  res2b=T2[5,1]

  T11=factor(polresultant(res1b,res2b,z))
  res1z=T11[3,1]*T11[4,1]
  
  R(y,z,b)=(1+y*z*b)*(1+y)*(1+y*z)*(1-y*z)*(1+z)*(1-z*b)*(1+b+b^2)*2*(1-b)*(1+b)

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
  /* 2 stationary point with value <= 4.9 */

  /* BORDER */
  /* 0-borders have already been discussed in greater generality (function is <= 8) */

  /* y = 1: already discussed in greater generality (it is <= 8.79) */

  /* z = 1: already discussed in greater generality (it is <= 8.79) */

  /* b = 1: NULL */
.........................................................................

  /* t = 1:
  BORDER

  b = 1: (1+y*z*a)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2 */

  L= (1+y*z*a)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2
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
  res1t=T11[3,1]

  R(y,z,a)=(1+y*z*a)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a)*(1+a+a^2)*(1+a)*(1-a)*2

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

  /* y = 1: already discussed in greater generality (it is <= 8.79) */

  /* z = 1: already discussed in greater generality (it is <= 8.79) */

  /* a = 1: NULL */
.........................................................................

/* BORDER
    a = 1: (1+y*z*t*b)*(1+y)*(1+y*z*t)*(1-y*z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
    
    CASE 1: substitution y - 1/2
*/

{res1Yta=
z^45 - 15592479/536350*z^44 + 4134050761/13408750*z^43 - 66022292961/67043750*z^42 - 
    8895760287987/1340875000*z^41 + 160998800484151/2681750000*z^40 - 6892571572641/167609375*z^39 -
    5933040684021203/5363500000*z^38 + 26260799267739171/8581600000*z^37 + 
    922269199876469321/85816000000*z^36 - 4534505232221212709/85816000000*z^35 - 
    2523081278702161477/42908000000*z^34 + 11897966889438313187/21454000000*z^33 + 
    655755968707040417/5363500000*z^32 - 45300867253128375023/10727000000*z^31 + 
    1003798637543614753/1340875000*z^30 + 135359405879597516903/5363500000*z^29 - 
    41730138031133496357/5363500000*z^28 - 1554456010609342977/12620000*z^27 + 
    85657540858780825083/2681750000*z^26 + 16704212791161179583/33521875*z^25 - 
    7332997003432585473/167609375*z^24 - 278837556254069724003/167609375*z^23 - 
    69779087445225024/265625*z^22 + 753755465531023625208/167609375*z^21 + 
    67045224411932080488/33521875*z^20 - 1599253035132761835408/167609375*z^19 - 
    244757505829933309312/33521875*z^18 + 2541853929936425059008/167609375*z^17 + 
    586495712115334275968/33521875*z^16 - 2788509546213073482624/167609375*z^15 - 
    4895491701784534219008/167609375*z^14 + 1687114145939762686976/167609375*z^13 + 
    5726394849177038045184/167609375*z^12 + 39699492349074440192/33521875*z^11 - 
    919713957668696948736/33521875*z^10 - 1398380684193903149056/167609375*z^9 + 
    2403924869666830286848/167609375*z^8 + 1262116921370981957632/167609375*z^7 - 
    731630613312688881664/167609375*z^6 - 110227554441472507904/33521875*z^5 + 
    99450939844597907456/167609375*z^4 + 110526879985165139968/167609375*z^3 - 
    5395312353589854208/167609375*z^2 - 8076080031782141952/167609375*z + 
    729583139634020352/167609375
;}

{res1Yb=
z^6*t^5 - 8/5*z^5*t^5 + 22/15*z^5*t^4 - 24/5*z^4*t^5 - 34/15*z^4*t^4 - 16/5*z^4*t^3 - 32/5*z^3*t^4 +
    68/15*z^3*t^3 - 64/15*z^3*t^2 + 32/3*z^2*t^3 + 16/3*z^2*t^2 + 32/15*z^2*t + 128/15*z*t^2 - 
    16/5*z*t + 32/15*z - 128/15*t - 32/15
;}

{LYb=
z^2*t^4*b^5 - 2/3*z^2*t^4*b^3 + 5/6*z^2*t^3*b^4 - 1/2*z^2*t^3*b^2 + 2/3*z^2*t^2*b^3 - 1/3*z^2*t^2*b 
    + 5/6*z*t^3*b^4 - 1/2*z*t^3*b^2 + 2/3*z*t^2*b^3 - 1/3*z*t^2*b + 1/2*z*t*b^2 - 1/6*z*t - 
    4/3*t^2*b^3 + 2/3*t^2*b - t*b^2 + 1/3*t - 2/3*b
;}

u=polrootsreal(res1Yta,[0,1]) /* length 3 */
z0=u[1]
v=polrootsreal(subst(res1Yb,z,z0),[0,1]) /* no roots */
z0=u[2]
v=polrootsreal(subst(res1Yb,z,z0),[0,1]) /* no roots */
z0=u[3]
v=polrootsreal(subst(res1Yb,z,z0),[0,1]) /* no roots */
.........................................................................

  /* BORDER
     a = 1: (1+y*z*t*b)*(1+y)*(1+y*z*t)*(1-y*z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b)
     INNER
     CASE 2
  */   

  R(y,z,t,b) = (1+y*z*t*b)*(1+y)*(1+y*z*t)*(1-y*z*t)*(1+z*t)*(1-z*t*b)*(1+t*b+(t*b)^2)*2*(1-b)*(1+b);

  /* The polynomials are quite large, we read them from a file */
  read("ROOT/Estimates/Proofs G/auxiliary_files_23G-C_Pari/auxiliary_23G-C_Pari_Q5.gp")
  
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
