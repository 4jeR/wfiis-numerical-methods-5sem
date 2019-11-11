#include <cmath>
#include <fstream>


// 1 ========================================= Euler =========================================

double GetY_t(double t, double lambda);
double EulerNextY(double yn, double dt, double lambda);
double EulerError(double y, double dt,double t, double lambda);
void EulerFunc(double dt, const char* filepath, const char* errorpath);
void EulerFuncAnalitics(double dt, const char* filepath);


// 2 ========================================= RK 2 =========================================

double RK2_NextY(double yn, double dt, double k1, double k2);
double RK2_GetK1(double yn, double lambda);
double RK2_GetK2(double yn, double dt, double lambda, double k1);
double RK2_Error(double y, double dt,double t, double lambda);
void RK2_Func(double dt, const char* filepath, const char* errorpath);
void RK2_FuncAnalitics(double dt, const char* filepath);

// 3 ========================================= RK 4 =========================================


double RK4_NextY(double yn, double dt, double k1, double k2);
double RK4_GetK1(double yn, double lambda);
double RK4_GetK2(double yn, double dt, double lambda, double k1);
double RK4_GetK3(double yn, double dt, double lambda, double k2);
double RK4_GetK4(double yn, double dt, double lambda, double k3);
double RK4_Error(double y, double dt,double t, double lambda, double k1, double k2, double k3, double k4);
void RK4_Func(double dt, const char* filepath, const char* errorpath);
void RK4_FuncAnalitics(double dt, const char* filepath);


// 4 ========================================= RRZ2 =========================================

void RRZ2_RLC(double ratio, const char* filepath);
double Voltage_t(double omega_v, double t);
double GetQI_n(double QIn, double dt, double k1, double k2, double k3, double k4);
double F(double t, double Q, double I);
double G(double t, double Q, double I, double omegav, double R, double L, double C);


int main (){
    // -===- metoda jawna Eulera -===-
    EulerFunc(0.01, "wynik-1a.txt", "error-1a.txt");
    EulerFunc(0.1, "wynik-1b.txt", "error-1b.txt");
    EulerFunc(1.0, "wynik-1c.txt", "error-1c.txt");
    EulerFuncAnalitics(0.01, "wynik-1-analitics.txt");


    // -===- metoda jawna RK2 -===-
    RK2_Func(0.01, "wynik-2a.txt", "error-2a.txt");
    RK2_Func(0.1, "wynik-2b.txt", "error-2b.txt");
    RK2_Func(1.0, "wynik-2c.txt", "error-2c.txt");
    RK2_FuncAnalitics(0.01, "wynik-2-analitics.txt");


    // -===- metoda jawna RK4 -===-
    RK4_Func(0.01, "wynik-3a.txt", "error-3a.txt");
    RK4_Func(0.1, "wynik-3b.txt", "error-3b.txt");
    RK4_Func(1.0, "wynik-3c.txt", "error-3c.txt");
    RK4_FuncAnalitics(0.01, "wynik-3-analitics.txt");


    // -===- RRZ 2 rzedu -===-
    RRZ2_RLC(0.5, "zad4_0.5wo.txt");
    RRZ2_RLC(0.8, "zad4_0.8wo.txt");
    RRZ2_RLC(1.0, "zad4_1.0wo.txt");
    RRZ2_RLC(1.2, "zad4_1.2wo.txt");



    return 0;
}


// ========================================= RRZ 2 rzedu =========================================


void RRZ2_RLC(double ratio, const char* filepath){
    double dt = std::pow(10.0, -4.0);
    double R = 100;
    double L = 0.1;
    double C = 0.001;
    double omega_zero = 1.0 / std::sqrt(L*C);
    double T_zero = 2.0 * M_PI / omega_zero;
    double omega_v = ratio * omega_zero;
    double curr_Q = 0.0;
    double curr_I = 0.0;

    std::ofstream file(filepath);

    for(double t = 0.0; t <= 4.0 * T_zero; t += dt){
    
        double k1_q = F(t, curr_Q, curr_I);
        double k1_i = G(t, curr_Q, curr_I, omega_v, R, L, C);

        double k2_q = F(t + dt / 2.0, curr_Q + dt * k1_q / 2.0, curr_I + dt * k1_i / 2.0);
        double k2_i = G(t + dt / 2.0, curr_Q + dt * k1_q / 2.0, curr_I + dt * k1_i / 2.0, omega_v, R, L, C);

        double k3_q = F(t + dt / 2.0, curr_Q + dt * k2_q / 2.0, curr_I + dt * k2_i / 2.0);
        double k3_i = G(t + dt / 2.0, curr_Q + dt * k2_q / 2.0, curr_I + dt * k2_i / 2.0, omega_v, R, L, C);
        
        double k4_q = F(t + dt , curr_Q + dt * k3_q , curr_I + dt * k2_i);
        double k4_i = G(t + dt , curr_Q + dt * k3_q , curr_I + dt * k2_i, omega_v, R, L, C);

        curr_Q = GetQI_n(curr_Q, dt, k1_q, k2_q, k3_q, k4_q);
        curr_I = GetQI_n(curr_I, dt, k1_i, k2_i, k3_i, k4_i);

        file << t << " " << curr_Q << " " << curr_I << "\n";
    }

    file.close();
}

double Voltage_t(double omega_v, double t){
    return 10.0 * std::sin(omega_v * t);
}


double GetQI_n(double QIn, double dt, double k1, double k2, double k3, double k4){
    return QIn + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

double F(double t, double Q, double I){
    return I;
}

double G(double t, double Q, double I, double omegav, double R, double L, double C){
    return (Voltage_t(omegav, t) / L ) - (R * I / L) - (Q / (L * C));
}

// ========================================= Euler =========================================

double EulerNextY(double yn, double dt, double lambda){
    return yn + dt * lambda * yn ;
}

double GetY_t(double t, double lambda){
    return std::exp(lambda * t);
}

double EulerError(double y, double dt,double t, double lambda){
    return EulerNextY(y, dt, lambda) - GetY_t(t, lambda);
}

void EulerFunc(double dt, const char* filepath, const char* errorpath){
    double curr_y = 1.0;
    double lambda = -1.0;
    
    std::ofstream file(filepath);
    std::ofstream file_err(errorpath);
   
    for(double t = dt; t < 5.0; t += dt){
        file << t << " " << GetY_t(t, lambda) << "\n"; 
        file_err << t << " " << EulerError(curr_y, dt, t, lambda) << "\n";
        curr_y = EulerNextY(curr_y, dt, lambda);
    }

    file_err.close();
    file.close();
}

void EulerFuncAnalitics(double dt, const char* filepath){
    double curr_y = 1.0;
    double lambda = -1.0;
    
    std::ofstream file(filepath);

    for(double t = 0.0; t < 5.0; t += dt){
        file << t << " " << EulerNextY(curr_y, dt, lambda) << "\n";
        curr_y = EulerNextY(curr_y, dt, lambda);
    }

    file.close();
}


// ========================================= RK 2 =========================================

// numeryczne
double RK2_NextY(double yn, double dt, double k1, double k2){
    return yn + dt * (k1 + k2) / 2.0;
}

double RK2_GetK1(double yn, double lambda){
    return lambda  * yn;
}

double RK2_GetK2(double yn, double dt, double lambda, double k1){
    return lambda * (yn + dt * k1);
}

double RK2_Error(double y, double dt,double t, double lambda, double k1, double k2){
    return RK2_NextY(y, dt, k1, k2) - GetY_t(t, lambda);
}

void RK2_Func(double dt, const char* filepath, const char* errorpath){
    double curr_y = 1.0;
    double lambda = -1.0;
    
    std::ofstream file(filepath);
    std::ofstream file_err(errorpath);
   
    for(double t = dt; t < 5.0; t += dt){
        double k1 = RK2_GetK1(curr_y, lambda);
        double k2 = RK2_GetK2(curr_y, dt, lambda, k1);
        file << t << " " << GetY_t(t, lambda) << "\n"; 
        file_err << t << " " << RK2_Error(curr_y, dt, t, lambda, k1, k2) << "\n";
        curr_y = RK2_NextY(curr_y, dt, k1, k2);
    }
    
    file_err.close();
    file.close();
}

void RK2_FuncAnalitics(double dt, const char* filepath){
    double curr_y = 1.0;
    double lambda = -1.0;
    
    std::ofstream file(filepath);

    for(double t = 0.0; t < 5.0; t += dt){
        double k1 = RK2_GetK1(curr_y, lambda);
        double k2 = RK2_GetK2(curr_y, dt, lambda, k1);
        file << t << " " << RK2_NextY(curr_y, dt, k1, k2) << "\n";
        curr_y = RK2_NextY(curr_y, dt, k1, k2);
    }

    file.close();
}



// ========================================= RK 4 =========================================


// numeryczne
double RK4_NextY(double yn, double dt, double k1, double k2, double k3, double k4){
    return yn + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

double RK4_GetK1(double yn, double lambda){
    return lambda  * yn;
}

double RK4_GetK2(double yn, double dt, double lambda, double k1){
    return lambda * (yn + dt * k1 / 2.0);
}

double RK4_GetK3(double yn, double dt, double lambda, double k2){
    return lambda * (yn + dt * k2 / 2.0);
}

double RK4_GetK4(double yn, double dt, double lambda, double k3){
    return lambda * (yn + dt * k3);
}

double RK4_Error(double y, double dt,double t, double lambda, double k1, double k2, double k3, double k4){
    return RK4_NextY(y, dt, k1, k2, k3 , k4) - GetY_t(t, lambda);
}

void RK4_Func(double dt, const char* filepath, const char* errorpath){
    double curr_y = 1.0;
    double lambda = -1.0;
    
    std::ofstream file(filepath);
    std::ofstream file_err(errorpath);
   
    for(double t = dt; t < 5.0; t += dt){
        double k1 = RK4_GetK1(curr_y, lambda);
        double k2 = RK4_GetK2(curr_y, dt, lambda, k1);
        double k3 = RK4_GetK3(curr_y, dt, lambda, k2);
        double k4 = RK4_GetK4(curr_y, dt, lambda, k3);
        file << t << " " << GetY_t(t, lambda) << "\n"; 
        file_err << t << " " << RK4_Error(curr_y, dt, t, lambda, k1, k2, k3, k4) << "\n";
        curr_y = RK4_NextY(curr_y, dt, k1, k2, k3, k4);
    }

    file_err.close();
    file.close();
}

void RK4_FuncAnalitics(double dt, const char* filepath){
    double curr_y = 1.0;
    double lambda = -1.0;
    
    std::ofstream file(filepath);

    for(double t = 0.0; t < 5.0; t += dt){
        double k1 = RK2_GetK1(curr_y, lambda);
        double k2 = RK2_GetK2(curr_y, dt, lambda, k1);
        file << t << " " << RK2_NextY(curr_y, dt, k1, k2) << "\n";
        curr_y = RK2_NextY(curr_y, dt, k1, k2);
    }

    file.close();
}
