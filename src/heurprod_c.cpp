#include <Rcpp.h>
using namespace Rcpp;

//' Heuristic Product Within Distance (Spatially Balanced Sampling).
//'
//' The function \code{hpwd} provides a fast selection of spatially balanced samples: it is a heuristic and a fast implemention of the algorithm \code{\link{pwd}}.
//' It generates samples approximately with the same probabilities of the \code{\link{pwd}} but with a significantly smaller number of steps.
//' In fact, this algorithm randomly selects a sample of size \eqn{n} exactly with \eqn{n} steps, updating at each step the selection probability of not-selected units, depending on their distance from the units, that were already selected in the previous steps.
//' Note that the inclusion probabilities are all constant and equal to \eqn{nsamp/N}, where \eqn{nsamp} is sample size and \eqn{N} is population size.
//'
//' @param dis A distance matrix NxN that specifies how far are all the pairs of units in the population.
//' @param nsamp Sample size.
//' @param nrepl Number of samples drawn.
//' @return Return a matrix 2x\code{nrepl} with \code{nrepl} samples drawn. In particular, the element \eqn{a_{ij}}{a_ij} is the j-th unit of the population drawn in the i-th sample.
//' @references
//' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
//'
//' \insertRef{fast_selection}{Spbsampling}
//' @examples
//' ##Example 1
//' ##Draw 50 samples of dimension 10 without constant probabilities
//' dis<-as.matrix(dist(cbind(lucas_abruzzo$x,lucas_abruzzo$y))) #distance matrix
//' drawn_samples<-hpwd(dis,10,50)
//' \donttest{
//' ##Example 2
//' ##Draw 100 samples of dimension 15 with constant probabilities equal to nsamp/N
//' #with N=population size
//' dis<-as.matrix(dist(cbind(lucas_abruzzo$x,lucas_abruzzo$y))) #distance matrix
//' vec<-rep(1,nrow(dis)) #vector of constraints
//' stand_dist<-stprod(dis,vec,1e-15,1000) #standardized matrix
//' drawn_samples<-hpwd(stand_dist,15,100)
//' }
//' @export
// [[Rcpp::export]]

NumericMatrix hpwd (NumericMatrix dis, int nsamp, int nrepl)
{
  NumericMatrix selez(nsamp*nrepl,2);
  int npop=dis.nrow();
  NumericVector r(1);
  NumericVector c(npop);
  NumericVector psel(npop);
  double drawn;
  for(int cc=1;cc<=nrepl;cc++)
  {
    for(int i=0;i<npop;i++)
    {
      psel(i)=1.0/npop;
    }
    for(int j=1;j<=nsamp;j++)
    {
      selez((cc-1)*nsamp+j-1,0)=cc;
      r=runif(1);
      for(int i=0;i<npop;i++)
      {
        c(i)=0;
      }
      drawn=0;
      c(0)=psel(0);
      if(c(0)>r(0))
      {
        drawn=0;
      }
      else
      {
      for(int i=1;i<npop;i++)
      {
        c(i)=c(i-1)+psel(i);
        if(c(i)>r(0))
        {
          drawn=i;
          break;
        }
      }
      }
      selez((cc-1)*nsamp+j-1,1)=drawn+1;
      for(int i=0;i<npop;i++)
      {
        psel(i)=psel(i)*dis(selez((cc-1)*nsamp+j-1,1)-1,i);
      }
      psel=psel/sum(psel);
    }
  }
  return selez;
}
