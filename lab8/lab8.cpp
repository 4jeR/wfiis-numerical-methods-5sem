#include <fstream>
#include <cmath>
#include <iostream>
#include <array>

constexpr int nx = 400;
constexpr int ny = 90;
constexpr int i_1 = 200;
constexpr int i_2 = 210;
constexpr int j_1 = 50;
constexpr double delta_triangle = 0.01;
constexpr double delta = 10.0 * delta_triangle;
constexpr double x_A = 0.45;
constexpr double y_A = 0.45;
constexpr int IT_MAX = 200;

using matrix = std::array<std::array<double, ny+1>, nx+1>;  

double t(int n);
double x(int i);
double y(int j);
void FillPsi(matrix& psi);
void FillV(matrix& vx, matrix& vy,const matrix& psi);
double CalcVmax(const matrix& vx, const matrix& vy);
double CalcDt(double vmax);
void FillU0(matrix& u0);
void PicardIteration(matrix& u0, matrix& u1, matrix& vx, matrix& vy, double dt, double D, const char* filepath);
void AdvectionDiffusion(double D, const char* filepath);

int main(){
    AdvectionDiffusion(0.1, "WynikD01.dat");
    AdvectionDiffusion(0.0, "WynikD0.dat");

    return 0;
}

void AdvectionDiffusion(double D, const char* filepath){
    matrix u0, u1, psi, vx, vy;
    FillPsi(psi);
    FillV(vx, vy, psi);
    double vmax = CalcVmax(vx, vy);
    double dt = CalcDt(vmax);
    FillU0(u0);
    PicardIteration(u0, u1, vx, vy, dt, D, filepath);
}

double CalcDt(double vmax){
    return delta_triangle / (4.0 * vmax);
}

double t(int n){
    return delta * n;
}

double x(int i){
    return delta_triangle *i;
}

double y(int j){
    return delta_triangle*j;
}

void FillPsi(matrix& psi){
    std::ifstream file;
    file.open("psi.dat");

    if(!file)
        std::cout << "Couldn't load file" << std::endl;
    else{
        int i, j;
        while (file >> i){
            file >> j;
            file >> psi[i][j];
        }
    }

    file.close();
}

void FillV(matrix& vx, matrix& vy,const matrix& psi){
    for(int i = 1; i <= nx-1; ++i){
        for(int j = 1; j <= ny-1; ++j){
            vx[i][j] =  (psi[i][j+1]-psi[i][j-1]) / (2.0 * delta_triangle);
            vy[i][j] = -(psi[i+1][j]-psi[i-1][j]) / (2.0 * delta_triangle);
        }
    }
    // zastawka
    for(int i = i_1; i <= i_2; ++i){
        for(int j = 0; j <= j_1; ++j){
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    // dolny i gorny brzeg
    for(int i = 1; i <= nx-1; ++i){
        vx[i][0] = 0.0;
        vx[i][ny] = 0.0;
        vy[i][0] = 0.0;
        vy[i][ny] = 0.0;
    }
    // lewy i prawy brzeg
    for(int j = 0; j <= ny; ++j){
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx-1][j];
    }
}

double CalcVmax(const matrix& vx, const matrix& vy){
    double vmax = 0.0;

    for(int i = 0; i <= nx; ++i){
        for(int j = 0; j <= ny; ++j){
            double vx2 = vx[i][j] * vx[i][j];
            double vy2 = vy[i][j] * vy[i][j];
            vmax = (std::sqrt(vx2+vy2) > vmax) ? std::sqrt(vx2+vy2) : vmax;
        }
    }

    return vmax;
}

void FillU0(matrix& u0){
    for(int i = 0; i <= nx; ++i){
        for(int j = 0; j <= ny; ++j){
            double up = (x(i)-x_A)*(x(i)-x_A) + (y(j)-y_A)*(y(j)-y_A);
            double down = 2.0*delta*delta;
            u0[i][j] = std::exp(-up/down) / (2.0 * M_PI * delta * delta);
        }
    }
}

void PicardIteration(matrix& u0, matrix& u1, matrix& vx, matrix& vy, double dt, double D, const char* filepath){
    std::ofstream file;
    file.open(filepath);

    for( int it=1; it<=IT_MAX; it++ ){
        for( int i=0; i<=nx; i++ ){
            for( int j=0; j<=ny; j++ ){
                u1[i][j]=u0[i][j];
            }
        }
        for( int k=1; k<=20; k++ ){
            for( int i=0 ;i<=nx; i++ ){
                for( int j=1 ;j<=ny-1; j++){
                    if( i<i_1 || i>i_2 || j>j_1 ){
                        if(i==0){
                            u1[i][j] = ( 1./( 1+( 2*D*dt / pow(delta, 2)) ) ) * ( u0[i][j] - (dt/2.) * vx[i][j] *
                            ( ( (u0[i+1][j] - u0[nx][j])/(2.*delta) ) + (u1[i+1][j] - u1[nx][j])/(2.*delta) ) - (dt / 2) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.*delta) + (u1[i][j+1] - u1[i][j-1])/(2.*delta) ) + dt/2. * D * 
                            ( ( u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(delta,2) + ( u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2) )
                            );
                        }
                        else if(i == 0 || i==nx){
                            u1[i][j] = ( 1./( 1+( 2*D*dt / pow(delta, 2)) ) ) * ( u0[i][j] - (dt/2.) * vx[i][j] *
                            ( ( (u0[0][j] - u0[i-1][j])/(2.*delta) ) + (u1[0][j] - u1[i-1][j])/(2.*delta) ) - (dt / 2) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.*delta) + (u1[i][j+1] - u1[i][j-1])/(2.*delta) ) + dt/2. * D * 
                            ( ( u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(delta,2) + ( u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2) )
                            );
                        }
                        else{
                            u1[i][j] = ( 1./( 1+( 2*D*dt / pow(delta, 2)) ) ) * ( u0[i][j] - (dt/2.) * vx[i][j] *
                            ( ( (u0[i+1][j] - u0[i-1][j])/(2.*delta) ) + (u1[i+1][j] - u1[i-1][j])/(2.*delta) ) - (dt / 2) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.*delta) + (u1[i][j+1] - u1[i][j-1])/(2.*delta) ) + dt/2. * D * 
                            ( ( u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(delta,2) + ( u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2) )
                            );
                        }
                    }
                }
            }
        }        
        
        for( int i=0; i<=nx; i++ ){
            for( int j=0; j<=ny; j++ ){
                u0[i][j]=u1[i][j];
            }
        }

        double c=0.;
        double xsr=0.;
        for( int i=0; i<=nx; i++ ){
            for( int j=0; j<=ny; j++ ){
                c += u0[i][j];
                xsr += x(i) * u0[i][j];
            }
        }
        file << it*dt << " " << c* pow(delta, 2) << " " << xsr*pow(delta, 2) << std::endl;
    }
    file.close();
    // mapy rozkladu?
}