#include <fstream>
#include <cmath>

// task 1 functions
void Picard(double beta, double N, double gamma, double t_max, double dt, double curr_u, double TOL, const char* filepath);
void Newton(double beta, double N, double gamma, double t_max, double dt, double curr_u, double TOL, const char* filepath);
double NextUPicard(double un, double dt, double gamma, double beta, double iter, double N);
double NextUNewton(double un, double prev_un, double dt, double gamma, double beta, double iter, double N);


// task 2 functions
void RK2_Method(double beta, double N, double gamma, double t_max, double dt, double currU1, double currU2, double TOL, const char* filepath);
double DeltaU1(double m11, double m12, double m21, double m22, double F1, double F2);
double DeltaU2(double m11, double m12, double m21, double m22, double F1, double F2);
double f(double t, double u, double beta, double N, double gamma);
double Korektor(double un, double dt, double b1, double b2, double f_val1, double f_val2);
double F1(double U1, double  U2, double un, double dt, double a11, double a12, double alfa, double beta);
double F2(double U1, double  U2, double un, double dt, double a21, double a22, double alfa, double beta);


int main(){
    // warunki poczatkowe

    double beta = 0.001;
    double N = 500;
    double gamma = 0.1;
    double t_max = 100.0;
    double dt = 0.1;

    double curr_u = 1.0;
    double TOL = std::pow(10.0, -6.0);


    Picard(beta, N, gamma, t_max, dt, curr_u, TOL, "picard.txt");
    Newton(beta, N, gamma, t_max, dt, curr_u, TOL, "newton.txt");
    RK2_Method(beta, N, gamma, t_max, dt, curr_u, curr_u, TOL, "rk2.txt");

    return 0;
}


// ========================== task 1 functions ==========================


void Picard(double beta, double N, double gamma, double t_max, double dt, double curr_u, double TOL, const char* filepath){
    double iter = 0.0;
    double diff = 0.0;

    std::ofstream file;
    file.open(filepath);

    for(double t = 0.0; t < t_max; t += dt){
        file << t << " " << curr_u << " " << N - curr_u << "\n";
        curr_u = NextUPicard(curr_u, dt, gamma,  beta, iter, N);
        diff = std::abs(NextUPicard(curr_u, dt, gamma, beta, iter, N) - curr_u);
        if(diff < TOL || iter > 20.0)
            break;
    }

    file.close();
}

void Newton(double beta, double N, double gamma, double t_max, double dt, double curr_u, double TOL, const char* filepath){
    double iter = 0.0;
    double diff = 0.0;
    double prev_u = curr_u;

    std::ofstream file;
    file.open(filepath);

    for(double t = 0.0; t < t_max; t += dt){
        file << t << " " << curr_u << " " << N - curr_u << "\n";
        curr_u = NextUNewton(curr_u, prev_u, dt, gamma,  beta, iter, N);
        prev_u = curr_u;
        diff = std::abs(NextUNewton(curr_u, prev_u, dt, gamma, beta, iter, N) - curr_u);
        if(diff < TOL || iter > 20.0)
            break;
    }

    file.close();
}


double NextUPicard(double un, double dt, double gamma, double beta, double iter, double N){
    double alfa = beta * N - gamma;
    double pow_un_iter = std::pow(un, iter);
    double value1 =  alfa * un - beta * un * un;
    double value2 = alfa * pow_un_iter - beta * pow_un_iter * pow_un_iter ; 
    return un + dt * (value1 + value2) / 2.0;    
}

double NextUNewton(double un, double prev_un, double dt, double gamma, double beta, double iter, double N){
    double alfa = beta * N - gamma;
    double un_to_iter = std::pow(un, iter);
    double value1 = alfa * un - beta * un * un;
    double value2 = alfa * un_to_iter - beta * un_to_iter * un_to_iter;

    double licznik = un - prev_un - dt * (value1 + value2) / 2.0;
    double mianownik = 1.0 - dt * (alfa - 2.0 * beta * un_to_iter) / 2.0;
    
    return un - (licznik / mianownik);
}


// ========================== task 2 functions ==========================


void RK2_Method(double beta, double N, double gamma, double t_max, double dt, double currU1, double currU2, double TOL, const char* filepath){
    double sqrt3 = std::sqrt(3.0);
    double alfa = beta * N - gamma;

    double a11 = 0.25;
    double a12 = 0.25 - sqrt3 / 6.0;
    double a21 = 0.25 + sqrt3 / 6.0;
    double a22 = 0.25;

    double b1 = 0.5;
    double b2 = 0.5;

    double c1 = 0.5 - sqrt3 / 6.0;
    double c2 = 0.5 + sqrt3 / 6.0;

    double curr_u = 1.0;

    std::ofstream file;
    file.open(filepath);

    for(double t = 0.0; t < t_max; t += dt){
        
        file << t << " " << curr_u << " " << N - curr_u << "\n";

        double f_val1 = f(t + c1 * dt, currU1, beta, N, gamma);
        double f_val2 = f(t + c2 * dt, currU2, beta, N, gamma);

        double m11 = 1.0 - dt * a11 * (alfa - 2.0 * beta * currU1);
        double m12 =     - dt * a12 * (alfa - 2.0 * beta * currU2);
        double m21 =     - dt * a21 * (alfa - 2.0 * beta * currU1);
        double m22 = 1.0 - dt * a22 * (alfa - 2.0 * beta * currU2);

        double F1_value = F1(currU1, currU2, curr_u, dt, a11, a12, alfa, beta);
        double F2_value = F2(currU1, currU2, curr_u, dt, a11, a12, alfa, beta);

        currU1 += DeltaU1(m11, m12, m21, m22, F1_value, F2_value);
        currU2 += DeltaU2(m11, m12, m21, m22, F1_value, F2_value);
        curr_u = Korektor(curr_u, dt, b1, b2, f_val1, f_val2); 

    }

    file.close();

}


double F1(double U1, double  U2, double un, double dt, double a11, double a12, double alfa, double beta){
    double val1 = a11 * (alfa * U1 - beta * U1 * U1);
    double val2 = a12 * (alfa * U2 - beta * U2 * U2);
    return U1 - un - dt * (val1 + val2);
}


double F2(double U1, double  U2, double un, double dt, double a21, double a22, double alfa, double beta){
    double val1 = a21 * (alfa * U1 - beta * U1 * U1);
    double val2 = a22 * (alfa * U2 - beta * U2 * U2);
    return U2 - un - dt * (val1 + val2);
}

double f(double t, double u, double beta, double N, double gamma){
    return (beta * N - gamma) * u - beta * u * u;
}


double DeltaU1(double m11, double m12, double m21, double m22, double F1, double F2){
    double licznik = F2 * m12 - F1 * m22;
    double mianownik = m11 * m22 - m12 * m21;

    return licznik / mianownik;
}

double DeltaU2(double m11, double m12, double m21, double m22, double F1, double F2){
    double licznik = F1 * m21 - F2 * m11;
    double mianownik = m11 * m22 - m12 * m21;

    return licznik / mianownik;
}



double Korektor(double un, double dt, double b1, double b2, double f_val1, double f_val2){
    return un + dt * (b1 * f_val1 + b2 * f_val2);
}

