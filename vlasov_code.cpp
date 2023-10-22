#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

float p0= 0.5;
float x0 = 1.0;
double x_min = -2.0;
double x_max = 2.0;
double p_min = -2.0;
double p_max = 2.0;
const int n_x = 10;
const int n_p = 10;
double epsilon = 1.0;
double d_t = 0.01;


double x_values[n_x];
double p_values[n_p];

double force[n_x];

// da las condiciones de borde en negativo
int b_c(int value)
{
    if (value < 0)
    {
        value = n_x + value;
    }
    return value;
}


// ---------------------------- waterbag array-----------------------------------------------------
void water_bag(double function[n_x][n_p])
{
    for (int j = 0; j < n_p; j++) {
        for (int i = 0; i < n_x; i++) {
            if (abs(x_values[i]) < x0) {
                if (abs(p_values[j]) < p0) {
                    function[i][j] = 1.0 / (4.0 * x0 * p0);
                } else {
                    function[i][j] = 0.0;
                }
            }
        }
    } 

}
// ----------------------------Creation of x and p values array-----------------------------------------------------
void init_xvalues()
{
    for (int i = 0; i < n_x; i++)
    {
        x_values[i] = x_min + (x_max - x_min) * i / (n_x - 1);
    }
}
void init_pvalues()
{
    for (int i = 0; i < n_p; i++)
    {
        p_values[i] = p_min + ((p_max - p_min) * i) / (n_p - 1);
    }
}



//--------------------------------FUNCIONES DEL CUBIC SPLINE-------------------------------------------------
double Alpha(double x, int i)
{
    return (x_values[(i + 1) % n_x] - x) / (x_values[(i + 1) % n_x] - x_values[i]);
}
double Beta(double x, int i)
{
    return 1 - Alpha(x, i);
}
double Gamma(double x, int i)
{
double alpha = Alpha(x, i);
    return ((pow(alpha, 3) - alpha) / 6) * pow((x_values[(i + 1) % n_x] - x_values[i]), 2);
}
double Delta(double x, int i)
{
    double beta = Beta(x, i);
    return ((pow(beta, 3) - beta) / 6) * pow((x_values[(i + 1) % n_x] - x_values[i]), 2);
}


//-------second derivate-------------------


double sec_derivate_array(double func_array[n_x][n_p],int i,int j)
{
    return -(1/560)*func_array[b_c(i-4)][j] + (8/315)*func_array[b_c(i-3)][j] - (1/5)*func_array[b_c(i-2)][j] + (8/5)*func_array[b_c(i-1)][j] - (205/72)*func_array[(i)%n_x][j] + (8/5)*func_array[(i+1)%n_x][j] - (1/5)*func_array[(i+2)%n_x][j] + (8/315)*func_array[(i+3)%n_x][j] - (1/560)*func_array[(i+4)%n_x][j];

}

double sec_derivate_array_plus(double func_array[n_x][n_p],int i,int j)
{
    return -(1/560)*func_array[b_c(i+1-4)][j] + (8/315)*func_array[b_c(i+1-3)][j] - (1/5)*func_array[b_c(i+1-2)][j] + (8/5)*func_array[i+1-1][j] - (205/72)*func_array[(i+1)%n_x][j] + (8/5)*func_array[(i+1+1)%n_x][j] - (1/5)*func_array[(i+1+2)%n_x][j] + (8/315)*func_array[(i+1+3)%n_x][j] - (1/560)*func_array[(i+1+4)%n_x][j];
}

//--------------------------------FUNCIONES DEL CUBIC SPLINE-------------------------------------------------

double Func_csp_array_x(double func_array[n_x][n_p],double x,int i,int j)
{
    return Alpha(x,i)*func_array[i][j]+ Beta(x,i)*func_array[(i+1)%n_x][j]+ Gamma(x,i)*sec_derivate_array(func_array,i,j) + Delta(x,i)*sec_derivate_array_plus(func_array,i,j)   ;
}
double Func_csp_array_p(double func_array[n_x][n_p],double p,int i,int j)
{
    return  Alpha(p,j)*func_array[i][j] + Beta(p,j)*func_array[i][(j+1)%n_p] + Gamma(p,j)*sec_derivate_array(func_array,i,j) + Delta(p,j)*sec_derivate_array_plus(func_array,i,j);

}


//------------------------------FALTA LA FUNCION DE LA INTEGRAL----------------------------------------------
//------------DEBE SER UNA INTEGRACION NUMERICA, CASO RAPIDO, TRAPECIOS, CASO LENTO, CUADRATURA--------------


// double integration()
// {
//     double h = (x_max - x_min)/ n_x    
//     for 

// }




// def integration_array_eval(functions_array): # requiere igualar a una variable
//     integral = np.zeros(n_x)#.astype('O')
//     for k in range(n_x):
//         for j in range(n_p):
//             for i in range(n_x):
//                 integral[k] = sp.integrate(-epsilon*functions_array[i][j].subs({x: x1, p: p1}), (x1, x_min,x_values[k]), (p1, p_min,p_max))+sp.integrate(-epsilon*-functions_array[i][j].subs({x: x1, p: p1}), (x1,x_values[k], x_max ), (p1,p_min, p_max))
//     return integral



//-----------------------------------FUNCIONES DEL LOOP--------------------------------------------------------

void func_one(double func_array[n_x][n_p], double func_one[n_x][n_p])
{
    for (int j = 0; j < n_p; j++) {
        for (int i = 0; i < n_x; i++) {
            func_one[i][j] = Func_csp_array_x(func_array,x_values[i]-((p_values[j]*d_t)/2)    ,i,j);
        }
    }
}
//- ((p_values[j] * d_t)/2)

void func_two(double func_array[n_x][n_p], double func_one[n_x][n_p],double force[n_x])
{
    for (int j = 0; j < n_p; j++) {
        for (int i = 0; i < n_x; i++) {
            func_one[i][j] = Func_csp_array_x(func_array,p_values[j] - (force[i] * d_t),i,j);
        }
    }
}



double main_func[n_x][n_p];
double func_one_array[n_x][n_p];

int main()
{ 

    init_xvalues();//inicia los valores de x
    init_pvalues();//inicia los valores de p 

   
    water_bag(main_func);// crea la funcion inicial 
    func_one(main_func,func_one_array);// realiza el primer paso, first advetion in x.

    printf("main_func \n");
     for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_p; j++) {
            //cout << main_func[(i+1)%n_x][j];
            printf("%f (%d %d) ", main_func[i][j],i,j);
        }
        printf("\n");
    }   

    printf("func_one_array \n");
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_p; j++) {
            //cout << main_func[(i+1)%n_x][j];
            printf("%f (%d %d) ", func_one_array[i][j],i,j);
        }
        printf("\n");
    }

    

}





    // for (int j = 0; j < n_p; j++) {
    //     for (int i = 0; i < n_x; i++) {
    //         printf("%i%i ", i,j);
    //     }
    //     printf("\n");
    // }    

    // //printf("%f", function[4][4]);    
    // // for (int j = 0; j < n_p; j++) {
    // //     for (int i = 0; i < n_x; i++) {
    // //         printf("%f ",  sec_derivate_array_plus(function,i,j));
    // //     }
    // //     printf("\n");
    // // } 


        // int func_array[n_x][n_p];
    // for (int i = 0; i < n_x; i++) {
    //     for (int j = 0; j < n_p; j++) {
    //         func_array[i][j] = i * 10 + j;
    //     }
    // }
    // printf("%d", func_array[2][1]);
//  //   Imprimir la matriz func_array en forma de matriz
    // for (int i = 0; i < n_x; i++) {
    //     for (int j = 0; j < n_p; j++) {
    //         printf("%d ", func_array[i][j]);
    //     }
    //     printf("\n");
    // }






    // int func_array[n_x][n_p];
    // for (int i = 0; i < n_x; i++) {
    //     for (int j = 0; j < n_p; j++) {
    //         func_array[i][j] = i * 10 + j;
    //     }
    // }

    // cout << b_c(1-3) << endl;
    // printf("%d", func_array[-3%n_x][1]);
