#include <fstream>
#include <cmath>
#include <iostream>
#include <array>

constexpr double delta = 0.01;
constexpr double p = 1.0;
constexpr double u = 1.0;
constexpr int n_x = 200;
constexpr int n_y = 90;
constexpr int i_1 = 50;
constexpr int j_1 = 55;
constexpr int IT_MAX = 20000;

using matrix = std::array<std::array<double, n_y+1>, n_x+1>;
using y_nodes = std::array<double, n_y+1>;
using x_nodes = std::array<double, n_x+1>;

void Relaxation(double Q_in, const char* filepath, const x_nodes& x, const y_nodes& y);
void SetBoundaryConditionsFi(double Qin, matrix& fi, const x_nodes& x, const y_nodes& y);
void SetBoundaryConditionsZeta(double Qin, const matrix& fi, matrix& zeta, const x_nodes& x, const y_nodes& y);
void ErrorControl(int iteration, const matrix& fi, const matrix& zeta);
void WriteResultsToFile(const char* filepath, const matrix& fi,const matrix& zeta, const x_nodes& x, const y_nodes& y);

double CalcFi(int i, int j, const matrix& fi, const matrix& zeta);
double CalcZeta(int i, int j, const matrix& fi, const matrix& zeta, double omega);
double CalcQ_out(double Qin, const y_nodes& y);

bool Inside(int i, int j);

int main(){
    x_nodes xnodes;
    y_nodes ynodes;
    
    for(int i = 0; i <= n_x; ++i)
        xnodes[i] = delta * i;
    
    for(int j = 0; j <= n_y; ++j)
        ynodes[j] = delta * j;
    
    Relaxation(-1000.0, "q-1000.txt", xnodes, ynodes);
    Relaxation(-4000.0, "q-4000.txt", xnodes, ynodes);
    Relaxation( 4000.0, "q+4000.txt", xnodes, ynodes);

    return 0;
}

void Relaxation(double Q_in, const char* filepath, const x_nodes& x, const y_nodes& y){
    double omega = 1.0;
    matrix fi,zeta;
    
    for(auto& arr : fi)
        for(auto& el : arr)
            el = 0.0;
        
    for(auto& arr : zeta)
        for(auto& el : arr)
            el = 0.0;
        
    SetBoundaryConditionsFi(Q_in, fi, x, y);

    for(int it = 1; it <= IT_MAX; ++it){
        omega = (it >= 2000) ? 1.0 : 0.0;
        
        for(int i = 1; i <= n_x-1; ++i){
            for(int j = 1; j <= n_y-1; ++j){
                if(Inside(i, j)){
                    fi[i][j] = CalcFi(i, j, fi, zeta);
                    zeta[i][j] = CalcZeta(i, j, fi, zeta, omega);
                }
            }
        }
        SetBoundaryConditionsZeta(Q_in, fi, zeta, x, y);
        // ErrorControl(it, fi, zeta);
    }

    WriteResultsToFile(filepath, fi, zeta, x, y);
}

bool Inside(int i, int j){
    return ( (0 < i   && i < n_x) && (0 < j   && j < n_y) 
        && 
          !( (0 <= i && i <= i_1) && (0 <= j && j <= j_1) ) 
    );
}

void SetBoundaryConditionsFi(double Qin, matrix& fi, const x_nodes& x, const y_nodes& y){
    // A
    for(int j = j_1; j <= n_y; ++j){
        fi[0][j] = Qin * (std::pow(y[j], 3.0)/3.0 - std::pow(y[j],2.0) * (y[j_1] + y[n_y]) / 2.0 + y[j]*y[j_1]*y[n_y]) / (2.0 * u);
    }
    // C
    double Q_out = CalcQ_out(Qin, y);
    for(int j = 0; j <= n_y; ++j){
        fi[n_x][j] =   Q_out * (std::pow(y[j], 3.0)/3.0 - std::pow(y[j], 2.0)*y[n_y]/2.0) / (2.0 * u)
                     +
                       Qin * std::pow(y[j_1], 2.0) * (-1.0 * y[j_1] + 3.0 * y[n_y]) / (12.0 * u);
    }
    // B
    for(int i = 1; i <= n_x-1; ++i){
        fi[i][n_y] = fi[0][n_y];
    }
    // D
    for(int i = i_1; i <= n_x-1; ++i){
        fi[i][0] = fi[0][j_1];
    }
    // E
    for(int j = 1; j <= j_1; ++j){
        fi[i_1][j] = fi[0][j_1];
    }
    // F
    for(int i = 1; i <= i_1; ++i){
        fi[i][j_1] = fi[0][j_1];
    }
}

void SetBoundaryConditionsZeta(double Qin, const matrix& fi, matrix& zeta, const x_nodes& x, const y_nodes& y){
    // A
    for(int j = j_1; j <= n_y; ++j){
        zeta[0][j] = Qin * (2.0 * y[j] - y[j_1] - y[n_y]) / (2.0 * u);
    }
    // C
    double Q_out = CalcQ_out(Qin, y);
    for(int j = 0; j <= n_y; ++j){
        zeta[n_x][j] = Q_out * (2.0 * y[j] - y[n_y]) / (2.0 * u);
    }
    // B
    for(int i = 1; i <= n_x-1; ++i){
        zeta[i][n_y] = 2.0 * (fi[i][n_y-1] - fi[i][n_y]) / (delta * delta);
    }
    // D
    for(int i = i_1+1; i <= n_x-1; ++i){
        zeta[i][0] = 2.0 * (fi[i][1] - fi[i][0]) / (delta * delta);
    }
    // E
    for(int j = 1; j <= j_1-1; ++j){
        zeta[i_1][j] = 2.0 * (fi[i_1+1][j] - fi[i_1][j]) / (delta * delta);
    }
    // F
    for(int i = 1; i <= i_1; ++i){
        zeta[i][j_1] = 2.0 * (fi[i][j_1+1] - fi[i][j_1]) / (delta * delta);
    }
    // E/F
        zeta[i_1][j_1] = 0.5 * (zeta[i_1-1][j_1] + zeta[i_1][j_1-1]);

}

double CalcFi(int i, int j, const matrix& fi, const matrix& zeta){
    return 0.25 * (   fi[i+1][j] 
                    + fi[i-1][j] 
                    + fi[i][j+1]
                    + fi[i][j-1]
                    - delta * delta * zeta[i][j]
    );
}

double CalcZeta(int i, int j, const matrix& fi, const matrix& zeta, double omega){
    double first = 0.25 * (   zeta[i+1][j] 
                            + zeta[i-1][j] 
                            + zeta[i][j+1] 
                            + zeta[i][j-1]
    );
    double second = omega * p * (  (fi[i][j+1] - fi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j])
                                 - (fi[i+1][j] - fi[i-1][j]) * (zeta[i][j+1] - zeta[i][j-1])
                                ) / (16.0 * u);
    return first - second;
}

[[maybe_unused]]
void ErrorControl(int iteration, const matrix& fi, const matrix& zeta){
    double sum = 0.0;
    double j2 = j_1 + 2;
    
    for(int i = 1; i <= n_x-1; ++i){
        sum += (  fi[i+1][j2] 
                + fi[i-1][j2] 
                + fi[i][j2+1] 
                + fi[i][j2-1] 
                - 4.0 * fi[i][j2] 
                - delta*delta*zeta[i][j2]);
    }

    std::cout << iteration << " -> " << sum << "\n";
}

void WriteResultsToFile(const char* filepath, const matrix& fi, const matrix& zeta, const x_nodes& x, const y_nodes& y){
    std::ofstream file;
    double u_xy, v_xy;

    file.open(filepath);

    for(int i = 0; i <= n_x; ++i){
        for(int j = 0; j <= n_y; ++j){
            if(Inside(i, j)){
                u_xy =  (fi[i][j+1] - fi[i][j-1])/(2.0 * delta);
                v_xy = -(fi[i+1][j] - fi[i-1][j])/(2.0 * delta);  
            }
            else{
                u_xy = 0.0;
                v_xy = 0.0;
            }
            file << x[i] << " " << y[j] << " " << fi[i][j] << " " << zeta[i][j] << " " << u_xy << " " << v_xy << "\n";
        }
        file << "\n";
    }

    file.close();
}

double CalcQ_out(double Qin, const y_nodes& y){
    double yny3 = std::pow(y[n_y], 3.0);
    return Qin * (yny3 - std::pow(y[j_1], 3.0) - 3.0 * y[j_1] * std::pow(y[n_y], 2.0) + 3.0 * std::pow(y[j_1], 2.0) * y[n_y] ) / yny3;
}

