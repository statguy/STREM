/*
Taken from:
http://www.cplusplus.com/forum/beginner/49408/#msg268694
*/

#include <Rcpp.h>
using namespace Rcpp;

double PerpDot(double a[2], double b[2]) { return a[1] * b[0] - a[0] * b[1]; }

int doIntersect(double A1[2], double A2[2], double B1[2], double B2[2]) {
  double a[2] = { A2[0] - A1[0], A2[1] - A1[1] };
  double b[2] = { B2[0] - B1[0], B2[1] - B1[1] };
  double f = PerpDot(a, b);
  if (f == 0) return false;
  
  double c[2] = { B2[0] - A2[0], B2[1] - A2[1] };
  double aa = PerpDot(a,c);
  double bb = PerpDot(b,c);
  
  if (f < 0) {
    if (aa > 0) return false;
    if (bb > 0) return false;
    if (aa < f) return false;
    if (bb < f) return false;
  }
  else {
    if (aa < 0) return false;
    if (bb < 0) return false;
    if (aa > f) return false;
    if (bb > f) return false;
  }
  
  return true;
}

// [[Rcpp::export]]
SEXP internal_countIntersections(Rcpp::NumericMatrix track, Rcpp::NumericMatrix surveyRoute) {
  int nTracks = track.nrow();
  int nSurveyRoutes = surveyRoute.nrow();
  int count = 0;
  
  if (nTracks == 1) return(wrap(0));
    
  for (int s = 1; s < nSurveyRoutes; s++) {
    double s1[2] = { surveyRoute(s-1, 0), surveyRoute(s-1, 1) };
    double s2[2] = { surveyRoute(s, 0), surveyRoute(s, 1) };
    
    for (int t = 1; t < nTracks; t++) {
      double t1[2] = { track(t-1, 0), track(t-1, 1) };
      double t2[2] = { track(t, 0), track(t, 1) };
      
      if (doIntersect(s1, s2, t1, t2) == true) count++;
    }
  }
  
  return(wrap(count));
}
