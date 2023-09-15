/* a = 1, Interior */

  {rest=
  z^13 - 1/3*z^12 + 14/405*z^11 - 14/1215*z^10 - 55907/54675*z^9 + 
      55907/164025*z^8 + 41996/273375*z^7 - 41996/820125*z^6 + 378133/1366875*z^5 
      - 378133/4100625*z^4 - 1802/18225*z^3 + 1802/54675*z^2 + 1/50625*z - 
      1/151875
  ;}
  
  {res1y=
  z^4*t^4 - 2/3*z^3*t^5 - 2/3*z^3*t^3 - 4/3*z^2*t^6 + 5/9*z^2*t^4 + 1/9*z^2*t^2 + 
      10/9*z*t^5 + 4/9*z*t^3 - 2/9*z*t - 1/9*t^4 - 1/3*t^2 + 1/9
  ;}
  
  {res1g=
  y*z^2*t^2 - 1/2*y*z*t + 1/2*z*t + 1/3*t^2 - 1/3
  ;}
  
  {Ly=
  y^2*z^2*t - 4/3*y*z*t*g + 2/3*y*z + 1/3*t - 2/3*g
  ;}
  
  R(y,z,t,g)=(1+y*z*t)*(1-2*y*z*g+(y*z)^2)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*(1-2*t*g+t^2)*2*sqrt(1-g^2);
  
  v=polrootsreal(rest,[0,1]);  /* Length 3 */
  z0=v[1];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[2];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[3];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 2 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* Length 1 */
  y0=u[1];
  ug=polrootsreal(subst(subst(subst(Ly,t,t0),z,z0),y,y0),[1/2,1]);  /* no roots */
  
  t0=w[2];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* Length 1 */
  y0=u[1];
  ug=polrootsreal(subst(subst(subst(Ly,t,t0),z,z0),y,y0),[1/2,1]);  /* no roots */
....................................................  

/* a = 1, g = 1/2 */
 /* inner */

{resz=
y^8 - 4*y^7 + 68/9*y^6 - 92/9*y^5 + 50/27*y^4 - 92/9*y^3 + 68/9*y^2 - 4*y + 1
;}
   
{res1t=
y*z + 1/2*y - 1/2
;}
   
{Ly=
y^2*z^2*t - 2/3*y*z*t + 2/3*y*z + 1/3*t - 1/3
;}
   
R(y,z,t)=(1+y*z*t)*(1-y*z+(y*z)^2)*(1+z+z^2)*(1-z*t)*(1+t+t^2)*(1-t+t^2)*2*sqrt(1-(1/2)^2)

  search(thresold)={
    tot=0;
    VY = polrootsreal(resz,[0,1]);
    if (#VY > 0,

       for(iy=1,#VY,
             VZ = polrootsreal(subst(res1t,y,VY[iy]),[0,1]);

             if(#VZ>0,
                for(iz=1,#VZ,
                      VT=polrootsreal(subst(subst(Ly,z,VZ[iz]),y,VY[iy]),[0,1]);     

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

  search(0);   
  
  /* 1 critical point, value <= 3.1 */
....................................................    

  /* a = 1, g = 1/2, z=1 */
L=(1+y*t)*(1-y+y^2)*3*(1-t)*(1+t+t^2)*(1-t+t^2)   /* times 2*sqrt(1-(1/2)^2) */

Ty = factor(deriv(L,y))
Tt = factor(deriv(L,t))
Ly=Ty[4,1];
Lt=Tt[1,1];
T=factor(polresultant(Ly,Lt,t))
res=T[2,1];
v=polrootsreal(res,[0,1]) /* length one */
y0=v[1];
u=polrootsreal(subst(Ly,y,y0),[0,1]) /* no roots */

/* y= 0:  3*(1-t)*(1+t+t^2)*(1-t+t^2)*2*sqrt(1-(1/2)^2)       <= 3*2*sqrt(3/4) <= 5.2 */
/* t= 0:  (1-y+y^2)*3*2*sqrt(1-(1/2)^2)                       <= 3*2*sqrt(3/4) <= 5.2 */
/* y= 1:  (1+t)*3*(1-t)*(1+t+t^2)*(1-t+t^2)*2*sqrt(1-(1/2)^2) <= 3*2*sqrt(3/4) <= 5.2 */
/* t= 1:  NULL */                      
....................................................    

/* a = 0, Interior */
  
  {rest=
  z^4 - 22/35*z^2 + 1/25
  ;}
  
  {res1y=
  z^2*t^2 - 8/3*z*t^3 + 2/3*z*t + 5/3*t^2 - 2/3
  ;}
  
  {res1g=
  y*z^2*t^2 - 1/2*y*z*t + 1/2*z*t + 1/3*t^2 - 1/3
  ;}

  {Ly=
  y^2*z^2*t - 4/3*y*z*t*g + 2/3*y*z + 1/3*t - 2/3*g
  ;}
  
  R(y,z,t,g)=(1+y*z*t)*(1-2*y*z*g+(y*z)^2)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*2*sqrt(1-g^2);
  v=polrootsreal(rest,[0,1]);  /* Length 2 */
  z0=v[1];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[2];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* Length 1 */
  y0=u[1];
  ug=polrootsreal(subst(subst(subst(Ly,t,t0),z,z0),y,y0),[1/2,1]);  /* No roots */
....................................................  
  
/* a = 0, g = 1/2 */

{resz=
y^4 - 2*y^3 - 2/3*y^2 - 2*y + 1
;}

{res1t=
y*z + 1/2*y - 1/2
;}

{Ly=
y^2*z^2*t - 2/3*y*z*t + 2/3*y*z + 1/3*t - 1/3
;}

  R(y,z,t)=(1+y*z*t)*(1-y*z+(y*z)^2)*(1+z+z^2)*(1-z*t)*(1+t+t^2)*2*sqrt(1-(1/2)^2)

u=polrootsreal(resz,[0,1]) /* length one */
y0=u[1];
v=polrootsreal(subst(res1t,y,y0),[0,1]) /* length one */
z0=v[1];
w=polrootsreal(subst(subst(Ly,z,z0),y,y0),[0,1]) /* length one */
t0=w[1];
R(y0,z0,t0) \\ 4.0246896762418487391946859391920716504
....................................................  

/* a = 0, g = 1/2, y = 1:
L= (1+z*t)*(1-z+z^2)*(1+z+z^2)*(1-z*t)*(1+t+t^2) /*  times 2*sqrt(1-(1/2)^2) */

Tz = factor(deriv(L,z))
Tt = factor(deriv(L,t))
Lz=Tz[3,1];
Lt=Tt[3,1];
T=factor(polresultant(Lz,Lt,t))
res=T[3,1];
v=polrootsreal(res,[0,1]) /* length two */
z0=v[1];
u=polrootsreal(subst(Lz,z,z0),[0,1]) /* length one */
t0=u[1];
my(t=t0,z=z0);eval(L)*2*sqrt(1-(1/2)^2)   \\  4.9573594494899125484118281554952721815

z0=v[2];
u=polrootsreal(subst(Lz,z,z0),[0,1]) /* length one */
t0=u[1];
my(t=t0,z=z0);eval(L)*2*sqrt(1-(1/2)^2)   \\  4.9509092222637217628497875385756421717

/* z= 0: (1+t+t^2)*2*sqrt(1-(1/2)^2)                         <= 5.2  */
/* t= 0: (1-z+z^2)*(1+z+z^2)*2*sqrt(1-(1/2)^2)               <= 5.2  */
/* z= 1: (1+t)*3*(1-t)*(1+t+t^2)*2*sqrt(1-(1/2)^2)           <= 6.85 */
/* t= 1: (1+z)*(1-z+z^2)*(1+z+z^2)*(1-z)*3*2*sqrt(1+(1/2)^2) <= 5.2  */                      
....................................................  

/* a = 0, g = 1/2, z = 1 */

L=(1+y*t)*(1-y+y^2)*3*(1-t)*(1+t+t^2) /* times  2*sqrt(1-(1/2)^2)  */
Ty = factor(deriv(L,y))
Tt = factor(deriv(L,t))
Ly=Ty[3,1];
Lt=Tt[1,1];
T=factor(polresultant(Ly,Lt,t))
res=T[2,1];
v=polrootsreal(res,[0,1]) /* length one */
y0=v[1];
u=polrootsreal(subst(Ly,y,y0),[0,1]) /* length one */
t0=u[1];
my(y=y0,t=t0);eval(L)*2*sqrt(1-(1/2)^2)   \\  4.3090261350217697335779317754610631178

/* y= 0:  3*(1-t)*(1+t+t^2)*2*sqrt(1-(1/2)^2)       <= 5.2 */
/* t= 0:  (1-y+y^2)*3*2*sqrt(1-(1/2)^2)             <= 5.2 */
/* y= 1:  (1+t)*3*(1-t)*(1+t+t^2)*2*sqrt(1-(1/2)^2) <= 5.2 */
/* t= 1:  NULL */                      
....................................................  

/* a = 0, y = 1  */
  
  {rest=
  z^10 + 13/33*z^8 - 86/825*z^6 - 2/15*z^4 - 1/275*z^2 + 3/275
  ;}
  
  {res1g=
  z^6*t^4 - 32/27*z^4*t^6 + 4/27*z^4*t^4 - 17/27*z^4*t^2 + 25/27*z^2*t^4 - 
      2/9*z^2*t^2 + 2/27*z^2 - 5/27*t^2 + 2/27
  ;}
  
  {Lt=
  z^2*t^3 + 3/2*z^2*t^2*g + 1/2*z^2*t - 1/2*t - 1/2*g
  ;}

  R(z,t,g)=(1+z*t)*(1-2*z*g+z^2)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*2*sqrt(1-g^2);
  
  v=polrootsreal(rest,[0,1]);  /* Length 2 */
  z0=v[1];
  w=polrootsreal(subst(res1g,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[2];
  w=polrootsreal(subst(res1g,z,z0),[0,1]);  /* Length 3 */
  t0=w[1];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* No roots */
  t0=w[2];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* No roots */
  t0=w[3];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[1/2,1]);  /* No roots */
....................................................
/* a = 0, y = 1, z = 1:

L=4*(1-g^2)*(1-t^2)*(1+2*t*g+t^2)  /* times 2*sqrt(1-g^2) */

Tt = factor(deriv(L,t))
Tg = factor(deriv(L,g)*(1-g^2)-g*L)
Lt=Tt[3,1];
Lg=Tg[5,1];
T=factor(polresultant(Lt,Lg,g))
res=T[4,1];
v=polrootsreal(res,[0,1]) /* length one */
t0=v[1];
u=polrootsreal(subst(Lt,t,t0),[1/2,1]) /* No roots */

/* t= 0:  4*(1-g^2)*2*sqrt(1-g^2) <= 5.2 */
/* g= 1/2:  4*(1-(1/2)^2)*(1-t^2)*(1+t+t^2)*2*sqrt(1-(1/2)^2) <= 6.85 */
/* t= 1:  NULL */
/* g= 1:  NULL */                      
....................................................
