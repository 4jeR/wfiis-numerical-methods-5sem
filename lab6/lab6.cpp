#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>
#include "mgmres.h"


constexpr double delta = 0.1;
double CalcRo1(double x, double y, const double xmax, const double ymax);
double CalcRo2(double x, double y, const double xmax, const double ymax);
double CalcRo(double x, double y, const double xmax, const double ymax);

constexpr std::size_t ij_to_l(std::size_t i, std::size_t j, std::size_t nx);
constexpr std::size_t l_to_j(std::size_t l, std::size_t nx);
constexpr std::size_t l_to_i(std::size_t l, std::size_t j, std::size_t nx);

void CSR_Algorithm(const char* fp_a, const char* fp_b, const char* fp_v, 
                   std::vector<double>& a, std::vector<int>& ia,
                   std::vector<int>& ja, std::vector<double>& b, std::vector<double>& v,
                   std::size_t nx, std::size_t ny, int& nz_num,
                   double V1, double V2, double V3, double V4, 
                   double eps_1, double eps_2, bool isRoZero);

void SolvePoisson(const char* filepath, std::size_t nx, std::size_t ny, int nz_num,
                  std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& a,
                  std::vector<double>& b, std::vector<double>& V);


int main(){
    std::vector<double> a;
    std::vector<int> ia;
    std::vector<int> ja;
    std::vector<double> b;
    std::vector<double> V;

    double eps1, eps2, V1, V2, V3, V4;
    std::size_t nx, ny;
    
    int nz_num = 0;


    eps1 = 1.0;
    eps2 = 1.0;
    V1 =  10.0;
    V2 = -10.0;
    V3 =  10.0;
    V4 = -10.0;
    nx = 4;
    ny = 4; 
    CSR_Algorithm("a1.txt","b1.txt", "v1.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, true);
    

    nx = 50;
    ny = 50;
    CSR_Algorithm("a2.txt","b2.txt", "v2.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, true);
    

    nx = 100;
    ny = 100;
    CSR_Algorithm("a3.txt","b3.txt", "v3.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, true);
    

    nx = 200;
    ny = 200;
    CSR_Algorithm("a4.txt","b4.txt", "v4.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, true);
    

    V1 = 0;
    V2 = 0;
    V3 = 0;
    V4 = 0;
    nx = 100;
    ny = 100;
    CSR_Algorithm("a5.txt","b5.txt", "v5.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, false);
    

    eps1 = 1.0;
    eps2 = 2.0;
    CSR_Algorithm("a6.txt","b6.txt", "v6.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, false);
    

    eps1 = 1.0;
    eps2 = 10.0;
    CSR_Algorithm("a7.txt","b7.txt", "v6.txt", a, ia, ja, b, V, nx, ny, nz_num, V1, V2, V3, V4, eps1, eps2, false);
    
    return 0;
}

void CSR_Algorithm(const char* fp_a, const char* fp_b,const char* fp_v, 
                   std::vector<double>& a, std::vector<int>& ia,
                   std::vector<int>& ja, std::vector<double>& b, std::vector<double>& v,
                   std::size_t nx, std::size_t ny, int& nz_num,
                   double v1, double v2, double v3, double v4, 
                   double eps_1, double eps_2, bool isRoZero)
{
    int k = -1;
    const std::size_t N = (nx + 1) * (ny + 1);

    a.clear();
    a.resize(5 * N);

    ia.clear();
    ia.resize(N+1);
    
    ja.clear();
    ja.resize(5 * N);

    b.clear();
    b.resize(N);


    std::vector<double> eps_vec;

    for(std::size_t l = 0; l < N+nx+1; ++l){
        std::size_t i = l_to_i(l, l_to_j(l, nx), nx);

        eps_vec.push_back((i <= nx / 2 ? eps_1 : eps_2));
    }

    for(std::size_t l = 0; l < N; ++l){
        std::size_t j_idx = l_to_j(l, nx);
        std::size_t i_idx = l_to_i(l, j_idx, nx);
        std::size_t brzeg = 0;
        double vb = 0.0;

        if(i_idx == 0){
            brzeg = 1;
            vb = v1;
        }
        else if(i_idx == nx){
            brzeg = 1;
            vb = v3;
        }

        if(j_idx == ny){
            brzeg = 1;
            vb = v2;
        }
        else if(j_idx == 0){
            brzeg = 1;
            vb = v4;
        }

        b[l] = (isRoZero) ? 0.0: -1.0 * CalcRo(delta * i_idx, delta * j_idx, delta * nx, delta * ny);

        if(brzeg == 1){
            b[l] = vb;
        }

        ia[l] = -1;
        // lewa skrajna przekatna
        if(l-nx-1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0){
                ia[l] = k;
            }

            a[k] = eps_vec[l] / (delta*delta);
            ja[k] = l - nx - 1;
        }
        // poddiagonala
        if(l - 1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0){
                ia[l] = k;
            }

            a[k] = eps_vec[l] / (delta*delta);
            ja[k] = l - 1;
        }
        // diagonala
        k++;
        if(ia[l] < 0){
            ia[l] = k;
        }

        if(brzeg == 0){
            a[k] = -(2.0 * eps_vec[l] * eps_vec[l+1] + eps_vec[l+nx+1]) / (delta * delta);
        }
        else{
            a[k] = 1;
        }
        ja[k] = l;
        // naddiagonala
        if(l < N && brzeg == 0){
            k++;
            a[k] = eps_vec[l+1] / (delta * delta);
            ja[k] = l+1;
        }

        // prawa skrajna przekatna
        if(l < N - nx - 1 && brzeg == 0){
            k++;
            a[k] = eps_vec[l+nx+1] / (delta * delta);
            ja[k] = l + nx + 1;
        }

    } // end of l loop

    nz_num = k + 1;
    ia[N] = nz_num;

    std::ofstream file;
    file.open(fp_a);
    for(std::size_t l = 0; l < N; ++l){
        std::size_t j_idx = l_to_j(l, nx);
        std::size_t i_idx = l_to_i(l, j_idx, nx);
        
        file << l << " " << i_idx << " " << j_idx << " " << a[l] << "\n";
    }
    file.close();

    file.open(fp_b);
    for(std::size_t l = 0; l < N; ++l){
        std::size_t j_idx = l_to_j(l, nx);
        std::size_t i_idx = l_to_i(l, j_idx, nx);
        
        file << l << " " << i_idx << " " << j_idx << " " << b[l] << "\n";
    }
    file.close();

    // SOLVE
    SolvePoisson(fp_v, nx, ny, nz_num, ia, ja, a, b, v);
}




double CalcRo1(double x, double y, const double xmax, const double ymax){
    double d = 0.1 * xmax;

    double first = (x - 0.25 * xmax) * (x - 0.25 * xmax) / (d*d); 
    double second = (y - 0.5 * ymax) * (y - 0.5 * ymax) / (d*d);

    return std::exp(- first - second);
}

double CalcRo2(double x, double y, const double xmax, const double ymax){
    double d = 0.1 * xmax;

    double first = (x - 0.25 * xmax)*(x - 0.25 * xmax) / (d*d); 
    double second = (y - 0.5 * ymax)*(y - 0.5 * ymax) / (d*d);

    return -1.0 * std::exp(- first - second);
}


double CalcRo(double x, double y, const double xmax, const double ymax){
    return CalcRo1(x, y, xmax, ymax) + CalcRo2(x, y, xmax, ymax);
}

constexpr std::size_t ij_to_l(std::size_t i, std::size_t j, std::size_t nx){
    return i + j * (nx+1);
}

constexpr std::size_t l_to_j(std::size_t l, std::size_t nx){
    return l / (nx+1);
}

constexpr std::size_t l_to_i(std::size_t l, std::size_t j, std::size_t nx){
    return l - j * (nx+1);
}




void SolvePoisson(const char* filepath, std::size_t nx, std::size_t ny, int nz_num, std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& a, std::vector<double>& b, std::vector<double>& V){
    const std::size_t N = (nx + 1) * (ny + 1);
    
    V.clear();
    V.resize(N);

    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1e-8;
    double tol_rel = 1e-8;  

    pmgmres_ilu_cr(N, nz_num, ia.data(), ja.data(), a.data(), V.data(), b.data(), itr_max, mr, tol_abs, tol_rel);
    
    std::ofstream file;
    file.open(filepath);
    for(std::size_t l = 0; l < N; l++){
        std::size_t j_idx = l_to_j(l, nx);
        std::size_t i_idx = l_to_i(l, j_idx, nx);
        
        file << delta * i_idx << " " << delta * j_idx << " " << V[l] << "\n";
    }
    file.close();
}