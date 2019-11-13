#include <fstream>
#include <cmath>
#include <utility>
#include <functional>


const double dt_zero = 1.0;
const double S = 0.75;
const double p = 2;
const double t_max = 40;
const double alfa = 5;
const double x_zero = 0.01;
const double v_zero = 0.0;


void TimeControl(std::function< std::pair<double, double> (double, std::pair<double,double>) > schemat_num, double TOL, const char* filepath);
std::pair<double, double> RK2(double dt, std::pair<double, double> xv_curr);
std::pair<double, double> Trapezy(double dt, std::pair<double, double> xv_curr);
double f(std::pair<double, double> xv_pair);
double g(std::pair<double, double> xv_pair);
double a11();
double a12(double dt);
double a21(double dt, std::pair<double, double>xv_next);
double a22(double dt, double xnext);
double F(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair);
double G(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair);
double CalcDx(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair);
double CalcDv(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair);

// ===================================== MAIN ===================================== 


int main(){

    TimeControl(RK2, std::pow(10.0, -2.0), "RK2_a.txt");
    TimeControl(RK2, std::pow(10.0, -5.0), "RK2_b.txt");

    TimeControl(Trapezy, std::pow(10.0, -2.0), "Trapez_a.txt");
    TimeControl(Trapezy, std::pow(10.0, -5.0), "Trapez_b.txt");


    return 0;
}

// ===================================== MAIN ===================================== 

void TimeControl(std::function< std::pair<double, double> (double, std::pair<double,double>) > schemat_num, double TOL, const char* filepath){
    double t = 0.0;
    double dt = dt_zero;
    std::pair<double, double> xv_curr(x_zero, v_zero);
    std::pair<double, double> xv_next1;
    std::pair<double, double> xv_next2;
    
    double Ex;
    double Ev;

    std::ofstream file;
    file.open(filepath);
    do{
        xv_next1 = schemat_num(dt, xv_curr);
        xv_next2 = schemat_num(dt, xv_next1);

        xv_next2 = schemat_num(2.0 * dt, xv_curr);

        Ex = (xv_next1.first  - xv_next2.first ) / (std::pow(2.0, p) - 1.0);
        Ev = (xv_next1.second - xv_next2.second) / (std::pow(2.0, p) - 1.0);
        
        if(std::max(std::abs(Ex), std::abs(Ev)) < TOL){
            t = t + 2.0 * dt;
            xv_curr = xv_next2;
            file << t << " " << dt << " " << xv_curr.first << " " << xv_curr.second << "\n";
        }
        double up = S * TOL;
        double down = std::max(std::abs(Ex), std::abs(Ev));
        dt = dt * std::pow(up / down, 1.0 / (p + 1.0));
    
    }while(t < t_max);
    file.close();
}


std::pair<double, double> Trapezy(double dt, std::pair<double, double> xv_curr){
    double sigma = std::pow(10.0, -10.0);
    double dx = 100.0;
    double dv = 110.0;
    std::pair<double, double> xv_next(xv_curr);

    while( std::abs(dx) > sigma || std::abs(dv) > sigma ){
        dx = CalcDx(dt, xv_curr, xv_next);
        dv = CalcDv(dt, xv_curr, xv_next);

        xv_next.first += dx;
        xv_next.second += dv;
    }

    return xv_next;
}

std::pair<double, double> RK2(double dt, std::pair<double, double> xv_curr){
    double thisx, thisv;
    thisx = xv_curr.second;
    thisv = alfa * (1.0 - xv_curr.first * xv_curr.first) * xv_curr.second - xv_curr.first;
    std::pair<double, double> k1 (thisx, thisv);
    thisx = xv_curr.second + dt * k1.second;
    thisv = alfa * (1.0 - std::pow(xv_curr.first + dt * k1.first, 2.0) * (xv_curr.second + dt * k1.second) - (xv_curr.first - dt * k1.first));
    std::pair<double, double> k2 (thisx, thisv);

    thisx = xv_curr.first + 0.5 * dt * (k1.first + k2.first);
    thisv = xv_curr.second + 0.5 * dt * (k1.second + k2.second);
    return std::pair<double, double>(thisx, thisv);
}

double f(std::pair<double, double> xv_pair){
    return xv_pair.second;
}

double g(std::pair<double, double> xv_pair){
    return alfa * (1.0 - xv_pair.first * xv_pair.first) * xv_pair.second - xv_pair.first;
}

double a11(){
    return 1.0;
}

double a12(double dt){
    return -0.5 * dt;
}

double a21(double dt, std::pair<double, double>xv_next){
    return -0.5 * dt * (-2.0 * alfa * xv_next.first * xv_next.second - 1.0);
}

double a22(double dt, double xnext){
    return 1.0 - 0.5 * dt * alfa * (1.0 - xnext * xnext);
}

double F(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair){
    return xv_next_pair.first - xv_pair.first - 0.5 * dt * (f(xv_pair) + f(xv_next_pair));
}

double G(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair){
    return xv_next_pair.second - xv_pair.second - 0.5 * dt * (g(xv_pair) + g(xv_next_pair));
}

double CalcDx(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair){
    double up = -1.0 * F(dt, xv_pair, xv_next_pair) * a22(dt, xv_next_pair.first) + G(dt, xv_pair, xv_next_pair) * a12(dt);
    double down = a11() * a22(dt, xv_next_pair.first) - a12(dt) * a21(dt, xv_next_pair);

    return up / down;
}

double CalcDv(double dt, std::pair<double, double> xv_pair, std::pair<double, double> xv_next_pair){
    double up = -1.0 * a11() * G(dt, xv_pair, xv_next_pair) + a21(dt, xv_next_pair) * F(dt, xv_pair, xv_next_pair); ;
    double down = a11() * a22(dt, xv_next_pair.first) - a12(dt) * a21(dt, xv_next_pair);

    return up / down;
}