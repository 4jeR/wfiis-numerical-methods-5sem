#include <fstream>
#include <cmath>
#include <iostream>
#include <array>

constexpr int nx = 128;
constexpr int ny = 128;
constexpr double delta = 0.2;
constexpr double xmax = delta * nx;
constexpr double ymax = delta * ny;
constexpr double TOL = 0.00000001;
unsigned iter = 0;

using matrix = std::array<std::array<double, ny+1>, nx+1>;

void SolvePoisson(int k, matrix& V, const char* filepath, const char* filepath2);
void WriteToFile(unsigned iteration, double S_value, std::ofstream& file);
void WriteToFile2(int k, matrix V, std::ofstream& file);
void SetBoundaryConds(matrix& V);
void ThickenMatrix(int k, matrix& V);
matrix& NextV(int k, matrix& V);
double CalcS(int k, const matrix& V);
double CalcValue(double S, double S_prev);


int main(){

    matrix V_matrix;

    for(int i = 0; i < ny+1; ++i){
        for(int j = 0; j < nx+1; ++j){
            V_matrix[i][j] = 0.0;
        }
    }

    SetBoundaryConds(V_matrix);
    SolvePoisson(16, V_matrix, "k16.txt", "map16.txt");
    SolvePoisson( 8, V_matrix,  "k8.txt",  "map8.txt");
    SolvePoisson( 4, V_matrix,  "k4.txt",  "map4.txt");
    SolvePoisson( 2, V_matrix,  "k2.txt",  "map2.txt");
    SolvePoisson( 1, V_matrix,  "k1.txt",  "map1.txt");

    return 0;
}

void SolvePoisson(int k, matrix& V, const char* filepath, const char* filepath2){
    double value;
    double S_now = CalcS(k, V);
    double S_prev = S_now;

    std::ofstream file, file2;
    file.open(filepath);
    file2.open(filepath2);
    
    do{
        V = NextV(k, V); // nastepna wartosc V
        S_prev = S_now;
        S_now = CalcS(k, V); // obliczenie sumy
        WriteToFile(iter, S_now, file); // iter - S(iter)
        value = CalcValue(S_now, S_prev); // warunek stopu
        ++iter;
    }while(value >= TOL);
    
    WriteToFile2(k, V, file2); // x - y - potencjal V 

    file.close();
    file2.close();

    ThickenMatrix(k, V); // zageszczenie siatki
}

void SetBoundaryConds(matrix& V){
    for(int j = 0; j < ny+1; ++j){
        double y = delta * j;
        V[0][j] = 1.0 * std::sin(M_PI * y / ymax);
    }

    for(int i = 0; i < nx+1; ++i){
        double x = delta * i;
        V[i][ymax] = -1.0 * std::sin(2.0 * M_PI * x / xmax);  
    }

    for(int j = 0; j < ny+1; ++j){
        double y = delta * j;
        V[xmax][j] = 1.0 * std::sin(M_PI * y / ymax); 
    }
    
    for(int i = 0; i < nx+1; ++i){
        double x = delta * i;
        V[i][0] = 1.0 * std::sin(2.0 * M_PI * x / xmax);  
    }
}


matrix& NextV(int k, matrix& V){   
    for(int i = k; i <= nx-1; i += k){
        for(int j = k; j<= ny-1; j += k){
            V[i][j] = 0.25 * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);    
        }
    }

    return V;
}


double CalcS(int k, const matrix& V){
    double sum = 0.0;

    for(int i = 0; i <= nx-1; i += k){
        for(int j = 0; j <= ny-1; j += k){
            double preBracket = std::pow(k * delta, 2.0) / 2.0;
            double down = 2.0 * k * delta;
            
            double one = (V[i+k][j] - V[i][j]) / down;
            double two = (V[i+k][j+k] - V[i][j+k]) / down;
            double three = (V[i][j+k] - V[i][j]) / down;
            double four = (V[i+k][j+k] - V[i+k][j]) / down;
            
            sum += preBracket * (std::pow(one+two, 2.0) + std::pow(three+four, 2.0)); 
        }
    }

    return sum;
}

double CalcValue(double S, double S_prev){
    return std::fabs( (S/S_prev) - 1.0 );
}


void WriteToFile(unsigned iteration, double S_value, std::ofstream& file){
    file << iteration << "\t" << S_value << "\n";
}



void WriteToFile2(int k, matrix V, std::ofstream& file){
    for(int i = 0; i < nx+1; i += k){
        for(int j = 0; j < ny+1; j += k){
            file << delta * i << "\t" << delta * j << "\t" << V[i][j] << "\n";
        }
    }
}

void ThickenMatrix(int k, matrix& V){   
    if(k != 1){
        for(int i=0; i <= nx-k; i += k){
            for(int j=0; j <= ny-k; j += k){
                V[i + k/2][j+k/2] = 0.25 * (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                V[i+k][j+k/2] = 0.5 * (V[i+k][j] + V[i+k][j+k]);
                V[i+k/2][j+k] = 0.5 * (V[i][j+k] + V[i+k][j+k]);
                V[i+k/2][j] = 0.5 * (V[i][j] + V[i+k][j]);
                V[i][j+k/2] = 0.5 * (V[i][j] + V[i][j+k]);
            }
        }
    }
}