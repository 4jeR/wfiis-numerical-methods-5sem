#include <fstream>
#include <cmath>
#include <iostream>
#include <array>
#include <vector>
#include <functional>



constexpr int nx = 150;
constexpr int ny = 100;
constexpr double eps = 1.0;
constexpr double delta = 0.1;
constexpr double xmax = delta * nx;
constexpr double ymax = delta * ny;
constexpr double TOL = 0.00000001;


std::array<std::array<double, ny+1>, nx+1> InitializeMatrix();
double CalcRo1(double x, double y);
double CalcRo2(double x, double y);
double CalcRo(double x, double y);

void GlobalPoisson(const char* filepath, const char* filepath2, double wG, const std::array<std::array<double, ny+1>, nx+1>& p);
void  LocalPoisson(const char* filepath, double wL, const std::array<std::array<double, ny+1>, nx+1>& p);

double CalcSum(const std::array<std::array<double, ny+1>, nx+1>& V, const std::array<std::array<double, ny+1>, nx+1>& p);
void PrintError(const std::array<std::array<double, ny+1>, nx+1>& V, const char* filepath);

int main(){

    
    std::array<std::array<double, ny+1>, nx+1> p = InitializeMatrix();

    GlobalPoisson("global-06.txt","map-global-06.txt", 0.6, p);
    GlobalPoisson("global-10.txt","map-global-10.txt", 1.0, p);

    LocalPoisson("local-10.txt", 1.0, p);
    LocalPoisson("local-14.txt", 1.4, p);
    LocalPoisson("local-18.txt", 1.8, p);
    LocalPoisson("local-19.txt", 1.9, p);



    return 0;
}


std::array<std::array<double, ny+1>, nx+1> InitializeMatrix() {
    std::array<std::array<double, ny+1>, nx+1> result;
    for (int i = 0; i < nx+1; i++) {
        for (int j = 0; j < ny+1; j++) {
            result[i][j] = CalcRo(i * delta, j * delta);
        }
    }
    return result;
}




double CalcRo1(double x, double y){
    double delta_x = 0.1 * xmax;
    double delta_y = 0.1 * ymax;
    double val1 = ((x - 0.35 * xmax) * (x - 0.35 * xmax)) / (delta_x * delta_x);
    double val2 = ((y - 0.5  * ymax) * (y - 0.5  * ymax)) / (delta_y * delta_y);
    
    return std::exp(-1.0 * val1 - val2);
}

double CalcRo2(double x, double y){
    double delta_x = 0.1 * xmax;
    double delta_y = 0.1 * ymax;
    double val1 = ((x - 0.65 * xmax) * (x - 0.65 * xmax)) / (delta_x * delta_x);
    double val2 = ((y - 0.5  * ymax) * (y - 0.5  * ymax)) / (delta_y * delta_y);
    
    return -std::exp(-1.0 * val1 - val2);
}


double CalcRo(double x, double y){
    return CalcRo1(x, y) + CalcRo2(x, y);
}

void GlobalPoisson(const char* filepath, const char* filepath2, double wG, const std::array<std::array<double, ny+1>, nx+1>& p){
    std::array<std::array<double, ny+1>, nx+1> Vs, Vn;
    
    for(auto& arr : Vn)
        for(auto& el : arr)
            el = 0.0;

    for(auto& arr : Vs)
        for(auto& el : arr)
            el = 0.0;

    
    // initial V values
    for(int i = 0; i <= nx; ++i){
        Vn[i][0] = Vs[i][0] = 10.0;
        Vn[i][ny] = Vs[i][ny] = 0.0;
    }


    double S_prev = 0.0;
    double S_curr = CalcSum(Vn, p);
    double value;
    long int iter = 0;

    std::ofstream file;
    file.open(filepath);
    // ITERATION OF 3 STEPS
    do{
        ++iter;
        S_prev = S_curr;




        // wzor 9 - STEP1
        for(int i = 1; i <= nx; ++i){
            for(int j = 1; j <= ny; ++j){
                Vn[i][j] = 0.25 * (Vs[i+1][j] + Vs[i-1][j] +
                                Vs[i][j+1] + Vs[i][j-1] + 
                                delta * delta * CalcRo(i, j) / eps);
            }
        }        
        // wzor 10, 11 - STEP2
        for(int j = 1; j <= ny-1; ++j){
            Vn[0][j] = Vn[1][j];
            Vn[nx][j] = Vn[nx-1][j];
        }
        //wzor 12 - STEP3
        for(int i = 0; i <= nx; ++i){
            for(int j = 0; j <= ny; ++j){
                Vs[i][j] = (1.0 - wG) * Vs[i][j] + wG * Vn[i][j];
            }
        } 



        S_curr = CalcSum(Vn, p);

        value = std::fabs((S_curr - S_prev) / (S_prev));
        file << iter << " " << S_curr << "\n";


    }while(value > TOL);
    
    PrintError(Vn, filepath2);
    file.close();
}

void LocalPoisson(const char* filepath, double wL, const std::array<std::array<double, ny+1>, nx+1>& p){
    std::array<std::array<double, ny+1>, nx+1> V;
    
    for(auto& arr : V)
        for(auto& el : arr)
            el = 0.0;

    

    // initial V values
    for(int i = 0; i <= nx; ++i){
        V[i][0] = 10.0;
        V[i][ny] = 0.0;
    }
    double S_prev = 0.0;
    double S_curr = CalcSum(V, p);
    double value;
    long int iter = 0;

    std::ofstream file;
    file.open(filepath);
    do{
        ++iter;

        // ITERATION OF 2 STEPS
        for(int i = 1; i <= nx; ++i){
            for(int j = 1; j <= ny; ++j){
                V[i][j] = (1.0 - wL) * V[i][j] + 0.25 * wL * 
                (V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + delta*delta * CalcRo(i,j) / eps);
            }
        }

        for(int j = 1; j <= ny-1; ++j){
            V[0][j] = V[1][j];
            V[nx][j] = V[nx-1][j];
        }


            
        S_prev = S_curr;
        S_curr = CalcSum(V, p);
        value = std::fabs((S_curr - S_prev) / (S_prev));
        file << iter << " " << S_curr << "\n";


    }while(value >= TOL);
    

    file.close();
}


double CalcSum(const std::array<std::array<double, ny+1>, nx+1>& V, const std::array<std::array<double, ny+1>, nx+1>& p){
    double sum = 0.0;
    for(int i = 0; i <= nx-1; ++i){
        for(int j = 0; j <= ny-1; ++j){
            double first = std::pow(( V[i+1][j] - V[i][j] ) / delta, 2.0);
            double second = std::pow(( V[i][j+1] - V[i][j] ) / delta, 2.0);;
            sum += delta*delta * (0.5 * first + 0.5 * second - p[i][j] * V[i][j]);
            }                                                                       
    }
    return sum;
}


void PrintError(const std::array<std::array<double, ny+1>, nx+1>& V, const char* filepath)
{
    std::ofstream file;
    file.open(filepath);

    for(int i = 1; i < nx; ++i){
        for(int j = 1; j < ny; ++j){
            double up = V[i+1][j] - 2.0 * V[i][j] + V[i-1][j] + V[i][j+1] - 2.0 * V[i][j] + V[i][j-1];
            double down = delta * delta;

            file << i * delta << " " << j * delta << " " << V[i][j] << " " << up / down + CalcRo(i, j) / eps << "\n";
        }
    }

    file.close();
}