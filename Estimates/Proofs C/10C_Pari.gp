/* y = 1 */

{resx=
z^5 - 2/3*z^4 + 43/54*z^3 - 1/3*z^2 + 7/54*z - 1/27
;}


{Lz=
x^2*z^3 - 1/6*x^2*z^2 - 1/2*x^2*z - 1/6*x*z^2 + 1/6*x - 1/2*z + 1/6
;}

R(x,z)=2*(1-x^2)*(1-x*z)*(1+x*z)^2*(1-z)*(1+z)^2;

v=polrootsreal(resx,[0,1]);  /* Length 1 */

z0=v[1];

w=polrootsreal(subst(Lz,z,z0),[0,1]);  /* Length 2 */

x0=w[1];

R(x0,z0)  /* = 2.44  > 64/27 */

x0=w[2];

R(x0,z0)  /* < 1.16 */


/* Interior */

{resz=
y^5 + 7/8*y^4 + 3/4*y^3 + 5/64*y^2 - 3/64*y - 9/128
;}


{res2x=
y^2*z^3 + 1/3*y^2*z^2 - 2*y^2*z - 2*y*z + 2/3*y - z + 1/3
;}


{Lz=
x^2*y^3*z^3 - 1/6*x^2*y^2*z^2 - 1/2*x^2*y*z - 1/6*x*y^2*z^2 + 1/6*x - 1/2*y*z + 
    1/6
;}


v=polrootsreal(resz,[0,1]);  /* Length 1 */

y0=v[1];

w=polrootsreal(subst(res2x,y,y0),[0,1]);  /* Length 1 */

z0=w[1];

u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[0,1]);

/* This has no roots */

