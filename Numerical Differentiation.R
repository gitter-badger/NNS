Newton <- function(f, x0=x0,N=N,tol=tol,h=h) {
  h=h
  i = 1
  x1 = x0
  p=numeric(N)

  while (i<=N) {
    df.dx = (f(x0+h)-f(x0))/h
    x1 = (x0 - (f(x0)/df.dx))
    p[i] = x1
    i = i + 1

    if(abs(x1-x0) < tol) break
    #if((x1/x0)<(1.1) && (x1/x0)>(.9)) break
    x0 = x1

  }
  return(p[1:(i-1)])
}


Numerical.differentiation = function(point,h,tol=1e-10,print.trace="FALSE"){

  options(digits=20)
  Bs = numeric()
  Bl = numeric()
  Bu = numeric()
  ### Step 1 initalize the boundaries for B

  ### Initial step size h
  f.x = f(point)
  f.x.h = f(point - h)

  ### Y = mX + B
  Y = f.x

  m = (f.x - f.x.h)/h

  mX = ((f.x - f.x.h)/h)*point

  B = f.x - ((f.x - f.x.h)/h)*point



  ### Initial interval for B given inputted h-step-value

  f.x.h.lower = f(point - h)
  f.x.h.upper = f(point + h)


  B1 = f.x - ((f.x - f.x.h.lower)/h)*point
  B2 = f.x - ((f.x.h.upper-f.x)/h)*point



  low.B = min(c(B1,B2))
  high.B = max(c(B1,B2))

  lower.B = low.B
  upper.B = high.B


  ## Return value at f(x) if lower.B and upper.B are identical to 20 digits
  if(lower.B==upper.B){

    par(mfrow=c(1,2))
    plot(f, xlim = c(point-(100*h), point+(100*h)),col='blue',ylab='f(x)')
    points(point,f.x, pch=19,col='red')
    plot(f, xlim = c(point-1, point+1),col='blue',ylab='f(x)')
    points(point,f.x, pch=19,col='red')
    return(c("Misspecified Scenario"))}


  new.B = mean(c(lower.B,upper.B))


  i=1
  while(i>=1L){
    Bl[i] = lower.B
    Bu[i] = upper.B
    new.f = function(x) -f.x + ((f.x - f(point - x))/x)*point + new.B

    ###  SOLVE FOR h, we just need the negative or positive sign from the tested B
    ###  This can be better served with an exact solution

    inferred.h = uniroot(new.f, c(-2*h,2*h))$root

    if(print.trace=="TRUE") {print(c("h"=inferred.h,"Lower B" = lower.B,"Upper B" = upper.B))}


    ## NARROW THE RANGE OF B BASED ON SIGN OF INFERRED.H
    if(B1==high.B){
      if(sign(inferred.h) < 0) {
        lower.B = new.B
        upper.B = upper.B
      }

      else {
        upper.B = new.B
        lower.B = lower.B
      }}  else{
        if(sign(inferred.h) < 0) {
          lower.B = lower.B
          upper.B = new.B
        }
        else {
          upper.B = upper.B
          lower.B = new.B
        }}



    new.B = mean(c(lower.B,upper.B))
    Bs[i]=new.B
    i = i+1


    ## Stop when the inferred h is within the tolerance level

    if(abs(inferred.h)<tol) {


      final.B = mean(c(upper.B,lower.B))
      slope = solve(point, f.x-final.B)
print(slope)
      par(mfrow=c(1,3))

      print(max(f((point-(100*h)):(point+(100*h)))))

      plot(f, xlim = c(point-(100*h), point+(100*h)),col='blue',ylab='f(x)',lwd=2,
           ylim = c(min(c(min(c(B1,B2)),min(na.omit(f((point-(100*h)):(point+(100*h))))))),max(na.omit(f((point-(100*h)):(point+(100*h)))))),
           main='f(x) and initial y-intercept range')
      abline(h=0,v=0)
      points(point,f.x, pch=19,col='red')
      points(x=rep(0,2),y=c(B1,B2),col='green',pch=19)
      segments(0,B1,point-h,f.x.h.lower,col='green')
      segments(0,B2,point+h,f.x.h.upper,col='green')


      plot(f, xlim = c(point-1, point+1),col='blue',ylab='f(x)',lwd=3,main='f(x) narrowed range and secant lines')
      points(point,f.x, pch=19,col='red')
      points(point-h,f.x.h.lower,col='green',pch=19)
      points(point+h,f.x.h.upper,col='green',pch=19)
      points(point,f.x, pch=19,col='red')
      segments(0,B1,point-h,f.x.h.lower,col='green')
      segments(0,B2,point+h,f.x.h.upper,col='green')

      plot(Bs,ylim=c(min(c(Bl,Bu)),max(c(Bl,Bu))),ylab="y-inetercept",col='green',pch=19,
           main='Iterated range of y-intercept')
      points(Bl,col='red',ylab='')
      points(Bu,col='blue',ylab='')



      return(as.matrix(c("Value of f(x) at point"=f(point),
                         "Final y-intercept (B)" = final.B,
                         "DERIVATIVE"=slope,
                         "Inferred h" = inferred.h,
                         "iterations"=i,
                         Finite.step(point,h)[1:2],
                         "Averaged Finite Step Initial h "=Finite.step(point,h)[3],
                         "Inferred h"=Finite.step(point,inferred.h)[1:2],
                         "Inferred h Averaged Finite Step"=Finite.step(point,inferred.h)[3])))
    }
  }

}






Finite.step = function(point,h){

  f.x = f(point)
  f.x.h.min = f(point - h)
  f.x.h.pos = f(point + h)

  neg.step = (f.x - f.x.h.min)/h
  pos.step = (f.x.h.pos - f.x)/h

  return((c("f(x-h)"=neg.step,"f(x+h)"=pos.step,mean(c(neg.step,pos.step)))))

}

