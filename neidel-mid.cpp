#include <iostream>
#include <cmath>


using namespace std;
double function(const double &x1,const double &x2){
  return x1*exp((7+x1)/2.)+pow(x2,2)*exp((7+x1)/2.)-8.*x2*exp((7+x1)/2.)+23.*exp((7+x1)/2.);
}

int main () {
   const int  n = 2;
   const double eps = 0.0001;
   double x[n + 5][n];
   double s;

   x[0][0] = -1; x[0][1] = 5;
   x[1][0] = 0;  x[1][1] = 1;
   x[2][0] = 1;  x[2][1] = 1;

    
   int h, l, g;
   long int kol = 0;

   do {
      kol++;
      h = 0;
      l = 0;
      g = 0;
      double max = function(x[0][0], x[0][1]);
      double min = function(x[0][0], x[0][1]);
      for (int k = 1; k <= n; k ++) {
         if (function(x[k][0], x[k][1]) < min) {
            min = function(x[k][0], x[k][1]);
            l = k;
         }
         if (function(x[k][0], x[k][1]) > max) {
            max = function(x[k][0], x[k][1]);
            h = k;
         }
      }

       max = function(x[0][0], x[0][1]);

      for (int k = 1; k <= n; k ++) {
         if ((function(x[k][0], x[k][1]) > max) && (k != h)) {
            max = function(x[k][0], x[k][1]);
            g = k;
         }
      }

      x[n + 1][0] = (x[0][0] + x[1][0] + x[2][0] - x[h][0]) / n;
      x[n + 1][1] = (x[0][1] + x[1][1] + x[2][1] - x[h][1]) / n;


      for (int i = 0; i < n; i++){
         x[n + 2][i] = 2 * x[n + 1][i] - x[h][i];
      }

      if (function(x[n + 2][0], x[n + 2][1]) < function(x[l][0], x[l][1])) {
         x[n + 3][0] = x[n + 1][0] + 2 * (x[n + 2][0] - x[n + 1][0]);
         x[n + 3][1] = x[n + 1][1] + 2 * (x[n + 2][1] - x[n + 1][1]);

      if (function(x[n + 3][0], x[n + 3][1]) < function(x[l][0], x[l][1])) {
         x[h][0] = x[n + 3][0];
         x[h][1] = x[n + 3][1];
      }else{
         x[h][0] = x[n + 2][0];
         x[h][1] = x[n + 2][1];
      }
   }
   else
      if (function(x[n + 2][0], x[n + 2][1]) > function(x[g][0], x[g][1])) {
         if (!(function(x[n + 2][0], x[n + 2][1]) > function(x[h][0], x[h][1]))) {
            x[h][0] = x[n + 2][0];
            x[h][1] = x[n + 2][1];
      }
         x[n + 4][0] = x[n + 1][0] + 0.5 * (x[h][0] - x[n + 1][0]);
         x[n + 4][1] = x[n + 1][1] + 0.5 * (x[h][1] - x[n + 1][1]);
            if (!(function(x[n + 4][0], x[n + 4][1]) > function(x[h][0], x[h][1]))) {
               x[h][0] = x[n + 4][0];
               x[h][1] = x[n + 4][1];
            }
   else {
      for (int k = 0; k < n; k++) {
         x[k][0] = x[k][0] + 0.5 * (x[k][0] - x[l][0]);
         x[k][1] = x[k][1] + 0.5 * (x[k][1] - x[l][1]);
      }
   }
   }
   else {
      x[h][0] = x[n + 2][0];
      x[h][1] = x[n + 2][1];
   }

      float s1 = 0, s2 = 0;
      for (int k = 0; k < n + 1; k++) {
         s1 += function(x[k][0], x[k][1]);
         s2 += function(x[k][0], x[k][1]) * function(x[k][0], x[k][1]);
      }
      s = s2 - s1 * s1 / (n + 1);
      s /= (n + 1);

      min = function(x[0][0], x[0][1]);
      for (int k = 1; k <= n; k++) {
         if (function(x[k][0], x[k][1]) < min) {
            min = function(x[k][0], x[k][1]);
            l = k;
         }
      }

   } while (s > eps);

   cout << "x1 = " << x[n + 1][0] << endl;
   cout << "x2 = " << x[n + 1][1] << endl;
   cout << "minF = " << function(x[n + 1][0], x[n + 1][1]) << endl;
   cout << "iterations = " << kol << endl;

}
