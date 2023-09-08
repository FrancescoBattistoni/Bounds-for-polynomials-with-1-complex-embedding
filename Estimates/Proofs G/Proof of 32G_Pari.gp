/* Study of inequality 32G */
\*

+' -  -             <= 13.7  when g in [1/2,1]
   -' -'
      + 2R^{1/2}
*/

PROOF

The function is (1-2*x*g+x^2)*(1+x*y)*(1+x*y*z)*(1+2*y*g+y^2)*(1+2*y*z*g+(y*z)^2)*(1-z)*2*sqrt(1-g^2). We have (x,y,z,g) in [0,1]^3 x [1/2,1].

* The function is strictly increasing in y, so that we can assume y = 1, getting

  <= (1-2*x*g+x^2)*(1+x)*(1+x*z)*2*(1+g)*(1+2*z*g+z^2)*(1-z)*2*sqrt(1-g^2).

==================
BORDER
==================

* CASE x = 0:
  The function with x = 0 becomes
  2*(1+g)*(1+2*z*g+z^2)*(1-z)*2*sqrt(1-g^2).
  The product (1+2*z*g+z^2)*(1-z) <= (1+2*z+z^2)*(1-z) <= 32/27, thus the function is
  <= 2*(1+g)*(32/27)*2*sqrt(1-g^2) <= 6.2

* CASE z = 0:
  The function with z = 0 becomes
  (1-2*x*g+x^2)*(1+x)*2*(1+g)*2*sqrt(1-g^2).
  The product (1-2*x*g+x^2)*(1+x) <= (1-x+x^2)*(1+x) <= 2   (recall we are assuming g >= 1/2), thus the function is
  <= 8*(1+g)*sqrt(1-g^2) <= 10.4

* CASE g = 1/2:
  The function with g = 1/2 becomes
  (1-x+x^2)*(1+x)*(1+x*z)*2*(3/2)*(1+z+z^2)*(1-z)*sqrt(3)
   = (1+x^3)*(1+x*z)*2*(3/2)*(1-z^3)*sqrt(3). This increases in x, so that
  <= 2*(1+z)*2*(3/2)*(1-z^3)*sqrt(3) <= 13.7.

* CASE x = 1:
  The function with x = 1 becomes
  2*(1-g)*2*(1+z)*2*(1+g)*(1+2*z*g+z^2)*(1-z)*2*sqrt(1-g^2)
  = 16*(1+2*z*g+z^2)*(1-z^2)*(1-g^2)^(3/2)
  <= 16
  PROOF (in PARI code)
  */ BORDER /*
  * x = 1, z = 0 already computed, <= 10.4
  
  * x = 1, g = 1/2 already computed, <= 13.7

  * x = 1, z = 1
    the function is 0

  * x = 1, g = 1
    the function is 0

  */ INNER /*
  */ PARI code /*
  L=16*(1+2*z*g+z^2)*(1-z^2);   \\ The function is L*(1-g^2)^(3/2)
  T=factor(deriv(L,z))
  Lz=T[1,1];
  T=factor(deriv(L,g)*(1-g^2)-3*L*g)
  Lg=T[3,1];
  T=factor(polresultant(Lz,Lg,z))
  res=T[4,1];  \\ 28*g^2-1               \\ there is also the factor g*(g-1) but we already know that in this case the value is <= 16
  \\ the solution is g= 1/sqrt(28) which is < 1/2
  \\ hence there are no stationary point in the open region
  END PROOF

* CASE z = 1:
  The function with z = 1  is NULL

* CASE g = 1:
  The function with g = 1 is NULL.



==================
INNER POINTS
==================
 L = (1-2*x*g+x^2)*(1+x)*(1+x*z)*2*(1+g)*(1+2*z*g+z^2)*(1-z)*2;   /* times sqrt(1-g^2) */
 T=factor(deriv(L,x))
 Lx=T[4,1];
 T=factor(deriv(L,z))
 Lz=T[3,1];
 T=factor(deriv(L,g)*(1-g^2) - g*L)
 Lg=T[5,1];
 T=factor(polresultant(Lx,Lz,z))
 res1z=T[3,1]*T[4,1];
 T=factor(polresultant(Lg,Lz,z))
 res2z=T[4,1];
 T=factor(polresultant(res1z,res2z,g))
 res1g=T[5,1]*T[6,1]*T[8,1]    \\ there is also other factors which are positive
 rx=polrootsreal(res1g,[0,1])                        /* four roots in the interval */

 x0=rx[1]   \\ 0.18439122973211267142107717449201977566
 subst(res2z,x,x0)  \\  -11.234630347711829197385072643455418169*g^6 + 74.888927235134765921983908169974940100*g^5 -183.95526021914602117638468314855452513*g^4 + 203.75597759790171355509857076995014793*g^3 - 101.02456702427942950742704118419712765*g^2 + 20.360767914461713229419075861073117609*g - 1.0250042058368676324708594940013248040
 rg=polrootsreal(subst(res2z,x,x0),[1/2,1])  \\ no roots in the interval

 x0=rx[2]   \\ 0.71207308530663467079235999224906767950
 subst(res2z,x,x0)  \\  -647.01107947000113263342188941994932413*g^6 + 2029.6669768609924850979075100288525322*g^5 -2137.5190050877161454854401944192392443*g^4 + 735.73477492605917124253937817030221768*g^3 + 60.707226689367426562201596598391430566*g^2 - 45.471118496744012572279848330899550067*g + 3.8956431508330297962172411308868211244
 rg=polrootsreal(subst(res2z,x,x0),[1/2,1])  \\ no roots in the interval

 x0=rx[3]   \\ 0.71428571428571428571428571428571428571
 subst(res2z,x,x0)  \\  -653.06122448979591836734693877551020408*g^6 + 2047.1470220741357767596834652228238234*g^5 -2153.6978639852442434699827452847028024*g^4 + 739.63807597174646618330797541840559631*g^3 + 61.781196610255930777142177154076957731*g^2 - 45.695382026196567756632015571743066240*g + 3.8914397912434444831660277605419510579
 rg=polrootsreal(subst(res2z,x,x0),[1/2,1])  \\ no roots in the interval

 x0=rx[4]   \\ 0.73588989435406735522369819838596156591
 subst(res2z,x,x0)  \\  -714.12875771394480225799889889351583802*g^6 + 2223.5497525959419830450811092294517860*g^5 -2317.0478226533840411131647428710220929*g^4 + 779.29585124248159422979733024905781328*g^3 + 72.344750017294448172323491931065936610*g^2 - 47.849215601412493709943736108271094360*g + 3.8374785049378095205918439242420103937
 rg=polrootsreal(subst(res2z,x,x0),[1/2,1])  \\ no roots in the interval

 /* NO STATIONARY POINTS */
END PROOF
