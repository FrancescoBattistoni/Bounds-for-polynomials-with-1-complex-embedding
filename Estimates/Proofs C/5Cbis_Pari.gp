
/* z = 1 t = 1 */

{resy= x - 5/9
;}

{res2g= x^3*y^4 - 3/4*x^3*y^2 + 1/2*x^2*y^2 - 1/4*x^2 - 3/4*x*y^2 + 1/2*x - 1/4
;}

{Lx= x^4*y^4 + 4/5*x^3*y^4 - 8/5*x^3*y^3*g - 6/5*x^2*y^3*g + 4/5*x*y*g + 2/5*y*g - 1/5
;}
   
v=polrootsreal(resy,[0,1]);   /* length one */
x0=v[1];
u=polrootsreal(subst(res2g,x,x0),[0,1])   /* no roots */
.........................................................................

/* t = 1 */

{resy=
z^8 - 583/65*z^7 + 97469/5850*z^6 - 6637/650*z^5 + 138298/26325*z^4 - 
    311266/26325*z^3 + 10417/810*z^2 - 1145/234*z + 70/117
;}


{res2x=
y^2*z^3 + 5/3*y^2*z^2 - 2/3*y*z + 2/3*y - 1/3*z + 1/3
;}


{res3g=
x^2*y^3*z^3 - 5/6*x^2*y^2*z^2 + 5/6*x*y^2*z^2 - 1/3*x*y*z - 1/6*x + 1/6
;}


{Lg=
x^2*y^2*z^2 - x*y^2*z^2 - 4*x*y*z*g - x + 1
;}


v=polrootsreal(resy,[0,1]); /* length two */
z0=v[1];
u=polrootsreal(subst(res2x,z,z0),[0,1])  /* no roots */

z0=v[2];
u=polrootsreal(subst(res2x,z,z0),[0,1])  /* no roots */
/* no stationary points */
.........................................................................

/* t = 1 g = 1/2 */

{resy= x^3 - 3/5*x^2 + 3/5*x - 4/5
;}

{res2z= x*y + 1/2*x - 1/2
;}

{Lz= x^3*y^3*z^3 - 1/2*x^3 + 1/2
;}
   
v=polrootsreal(resy,[0,1]);   /* length one */
x0=v[1];
u=polrootsreal(subst(res2z,x,x0),[0,1])   /* length one */
y0=u[1];
w=polrootsreal(subst(subst(Lz,y,y0),x,x0),[0,1])   /* no roots */
/* no stationary points */
.........................................................................
