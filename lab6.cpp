#include<iostream>
#include<cmath>
#include<fstream>


double function(const double &x1,const double &x2){
  return x1*exp((7+x1)/2.)+pow(x2,2)*exp((7+x1)/2.)-8.*x2*exp((7+x1)/2.)+23.*exp((7+x1)/2.);
}
double function1(const double &x){
  return x*exp((7+x)/2.)+pow(x,2)*exp((7+x)/2.)-8.*x*exp((7+x)/2.)+23.*exp((7+x)/2.);
}

double partial_f1_x1(const double  &x1,const double  &x2){
  return (x1*exp((7+x1)/2.))/2. + (pow(x2,2) * exp((7+x1)/2.))/2. -4*x2*exp((7+x1)/2.) + (25*exp((7.+ x1) /2.))/2.;
}

double partial_f1_x2(const double &x1,const double &x2){
  return 2 * x2 * exp((7.+x1)/2.) - 8 * exp((7.+x1)/2);
}

double partial_f2_x1(const double &x1,const double &x2){
  return x1*exp((7+x1)/2.)/4. + (pow(x2,2) * exp((7+x1)/2.))/4.-2*x2*exp((7+x1)/2.)+(27*exp((7+x1)/2.))/4.;
}

double partial_f2_x2(const double &x1,const double &x2){
  return 2 * exp((7+x1)/2.);
}

double partial_f_x1_x2(const double  &x1,const double  &x2){
  return x2*exp((7+x1)/2.) - 4*exp((7+x1)/2.);
}


double **inverse_matrix(double** const A, int N){
  double temp;
  double** UnitMatrix = new double* [N];
  double** Inverse_Matrix = new double* [N];
  for (int i = 0; i < N; i++){
		UnitMatrix[i] = new double[N];
    Inverse_Matrix[i] = new double[N];
  }
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      Inverse_Matrix[i][j] = A[i][j]; // <---
        UnitMatrix[i][j] = 0.0;
          if (i == j){
            UnitMatrix[i][j] = 1.0;
          }
        }
      }
    //-------------
    for (int k = 0; k < N; k++) {
        temp = Inverse_Matrix[k][k];
        for (int j = 0; j < N; j++){
            Inverse_Matrix[k][j] /= temp;
            UnitMatrix[k][j] /= temp;
        }
        for (int i = k + 1; i < N; i++){
            temp =  Inverse_Matrix[i][k];
            for (int j = 0; j < N; j++){
                Inverse_Matrix[i][j] -= Inverse_Matrix[k][j] * temp;
                UnitMatrix[i][j] -=  UnitMatrix[k][j] * temp;
            }
        }
    }
    for (int k = N - 1; k > 0; k--){
        for (int i = k - 1; i >= 0; i--){
            temp = Inverse_Matrix[i][k];
            for (int j = 0; j < N; j++){
                Inverse_Matrix[i][j] -= Inverse_Matrix[k][j] * temp;
                UnitMatrix[i][j] -= UnitMatrix[k][j] * temp;
            }
        }
    }
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			  Inverse_Matrix[i][j] =  UnitMatrix[i][j];
      }
    }


	for (int i = 0; i < N; i++){
		delete[]  UnitMatrix[i];
  }
	delete[]  UnitMatrix;
  return Inverse_Matrix;
}


double **matrix_gesse(const double  &x1,const double  &x2){
  const int size = 2;
  double **gesse =  new double *[size];
  double **inv_gesse = new double *[size];
  for(int i = 0; i < size; i++){
        gesse[i] = new double[size];
        inv_gesse[i] = new double[size];
      }

  gesse[0][0] = partial_f2_x1(x1,x2);
  gesse[1][1] = partial_f2_x2(x1,x2);
  gesse[1][0] = partial_f_x1_x2(x1,x2);
  gesse[0][1] = partial_f_x1_x2(x1,x2);

  return gesse;

}

double *swanna_x1(const double &x1,const double &x2){
    double h= 0.001;
    double *a_x1 = new double[2];
    a_x1[0] = x1;

    a_x1[1] = a_x1[0] + h;
    if (function(a_x1[1], x2) > function(a_x1[0], x2))
    {
        h *=-1;
        a_x1[1] = a_x1[0] + h;
    }
    while (function(a_x1[1], x2) < function(a_x1[0], x2))
    {
        a_x1[0] = a_x1[1];
        a_x1[1] = a_x1[0] + h;
        h *= 2;
    }
  return a_x1;
}

double *swanna_x2(const double &x1,const double &x2){
    double h = 0.001;
    double *a_x2 = new double[2];
    a_x2[0] = x2;
    a_x2[1] = a_x2[0] + h;
    if (function(x1,a_x2[1]) > function(x1,a_x2[0]))
    {
        h *=-1;
        a_x2[1] = a_x2[0] + h;
    }
     while (function(x1,a_x2[1]) < function(x1,a_x2[0]))
    {
        a_x2[0] = a_x2[1];
        a_x2[1] = a_x2[0] + h;
        h *= 2;
    }

    return a_x2;
}

double *method_uniform(double a, double b){
    double N=1000;
    double aa = a, bb = b; 
    double dn = abs(b - a)/N;
    double temp=0;
    double min=0;
    double *x = new double [2];
    x[0] = a;
    x[1] = b;
    min = function(a, b);
    for(double i = a; i < b; i += dn){
      for(double j = b; j > a; j -= dn){
        temp = function(i,j);
        if(temp < min){
          min = temp;
          x[0] = i;
          x[1] = j;
        }
      }
    }
    return x;
  
}

int main(){

  std::ofstream out;


  double x1 = -1;
  double x2 = 5;
  const int size = 2;
  const double eps = 0.00001;


  double **inv_gesse = new double *[size];
  double **gesse =  new double *[size];
  double *grad =  new double [size];
  double *P  = new double [size];
  double *uniform = new double [size];

  for(int i = 0; i < size; i++){
        gesse[i] = new double[size];
        inv_gesse[i] = new double[size];
      }
  gesse = matrix_gesse(x1,x2);
  inv_gesse = inverse_matrix(gesse,size);
  grad[0] = partial_f1_x1(x1,x2);
  grad[1] = partial_f1_x2(x1,x2);


double x[size];
x[0] = x1; // -1
x[1] = x2; // 5
double det_grad = sqrt(grad[0] * grad[0] + grad[1] * grad[1]);
double temp[size];
double *swanna_1 = new double[size];
double *swanna_2 = new double[size];
int n = 0;





out.open("data.dat");
while(det_grad >= eps){
  
  for(int i = 0; i < size; i++){
     P[i] = 0;
        for(int j = 0; j < size; j++){  
          P[i] +=  grad[j] * inv_gesse[i][j]; //направление спуска
        }
      }

  

  for(int i = 0; i < size; i++){
    temp[i] = x[i] - 1*P[i];
    x[i] = temp[i];
  }

  gesse = matrix_gesse(x[0],x[1]);
  inv_gesse = inverse_matrix(gesse,size);
  grad[0] = partial_f1_x1(x[0],x[1]);
  grad[1] = partial_f1_x2(x[0],x[1]);
  det_grad = sqrt(grad[0] * grad[0] + grad[1] * grad[1]);
  swanna_1 = swanna_x1(x[0],x[1]);
  swanna_2 = swanna_x2(x[0],x[1]);
  uniform = method_uniform(x[0],x[1]);
  out << x[0] << " " << x[1] << std::endl;
  n++;
}
  out.close();
  std::cout << "=======N-R======" << std::endl;
  std::cout <<" x1 = " << x[0] << " x2 = " << x[1] << " \n Iteration: " << n  << std::endl;
  std::cout <<" f(x1,x2) = " << function(x[0],x[1]) << std::endl;
  std::cout << "--Swanna--" << std::endl;
  std::cout <<"[" <<   swanna_1[0] << " : " << swanna_1[1] << "]" << std::endl;
  std::cout <<"[" <<   swanna_2[0] << " : " << swanna_2[1] << "]" << std::endl;
  std::cout << "===Method uniform===" << std::endl;
  std::cout << uniform[0] << " " << uniform[1] << std::endl;






for (int i = 0; i < size; i++){
		delete[]  gesse[i];
    delete[]  inv_gesse[i];
  }
	delete[]  gesse;
  delete[]  inv_gesse;
  delete[]  grad;
  delete[]  P;

  return 0;
}
