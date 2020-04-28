#include<iostream>
#include<cmath>



double function(const double &x1,const double &x2){
  return x1*exp((7+x1)/2.)+pow(x2,2)*exp((7+x1)/2.)-8.*x2*exp((7+x1)/2.)+23.*exp((7+x1)/2.);
}

double fRand(double fMin, double fMax)
{  
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double ramdom_search(const double  &x1,const double  &x2){
    double eps = 0.0000001;
    double lambda = 1;
    double alfa = 2;
    int N  = 10;
    int err_N = 3*N;

    double x[2];

    int g = 1000;   
    double game_x1[g] = {0};
    double game_x2[g] = {0};
    double game_f[g] = {0};

    double temp_x1 = x1;
    double temp_x2 = x2;
    double min = function(temp_x1,temp_x2);
    double unit_vector[2] = {0};
    int i= 0;
    for(int j = 0; j < g; j++){ //Небольшая модификация, 10 раудов отбора 
        while(lambda > eps && i < 3*N){
            unit_vector[0] = fRand(-N,N) * lambda;
            unit_vector[1] = fRand(-N,N) * lambda;

            x[0] = temp_x1 +  unit_vector[0];
             x[1] = temp_x2 +  unit_vector[1];

            if(min > function(x[0],x[1])){
                temp_x1 = x[0];
                temp_x2 = x[1];
                min = function(x[0],x[1]);
                lambda /= alfa;
                i = 0;
            }else{
                i++;
            }
        }
        game_x1[j] = x[0];
        game_x2[j] = x[1];
        game_f[j] = min;
    }
    for(int i = 0; i < g; i++){
        for(int j = 0; j < g; j++){
             if(game_f[i] > game_f[j]){
                 min = game_f[j];
                 x[0] =  game_x1[j];
                 x[1] =  game_x2[j];
             }
        }
    }

    std::cout << x[0] << " " << x[1] << " " << function(x[0],x[1]) << std::endl;
    return 0;
}


int main(){
    srand(time(NULL));
    double x1 = -1, x2 = 5;
    ramdom_search(x1,x2);
    return 0;
}