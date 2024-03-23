#include <iostream>
#include <clocale>
#include <iomanip>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <locale.h>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>

#include <gtest/gtest.h>  

using namespace std;
using namespace pde_solvers;

struct pipe {
    double L;
    double d_vnesh;
    double b;
    double sher;
    double z_0;
    double z_l;
    double n;
    double density;
    double nu;
    double p_0;
    double v;
    double Q;
    double resistance;

    double get_inner_diameter() const {
        return d_vnesh - 2 * b;
    }
    double get_relative_roughness() const {
        return sher / get_inner_diameter();
    }

    double get_inner_area() const {
        double D = get_inner_diameter();
        double S = M_PI * D * D / 4;
        return v * S;
    }

    double get_v() const {
        double D = get_inner_diameter();
        return 4 * Q / (3.1415 * pow(D, 2));
    }

    double get_Re() const {
        double D = get_inner_diameter();
        return get_v() * D / nu;
    }

    double get_t_w() const {
        return resistance / 8 * density * pow(get_v(), 2);
    }

    double get_dx() const {
        return L / n;
    }

    double get_dt() const {
        return get_dx() / v;
    }

    double get_n() const {
        return (L / v) / get_dt();
    }
};

struct massiv {
    vector<double> massiv;
};

void characteristic_method(pipe myPipe, double parametr, vector<double>& current_layer, vector<double>& previous_layer) {
    current_layer[0] = parametr;
    for (size_t j = 1; j < myPipe.n; j++) {
        for (size_t i = 1; i < myPipe.n; i++) {
            current_layer[i] = previous_layer[i - 1];
        }
    }
}

void out_put(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, int i, massiv time) {
    vector<vector<double>>& current_layer = buffer.current();


    double p_0 = myPipe.p_0;

    if (i == 0) {
        ofstream outFile("block_3.csv");
        outFile << "Время,Координата,Плотность,Вязкость,Давление" << "\n";
        // Записать значения текущего слоя в файл

        for (size_t j = 1; j < current_layer[0].size(); j++) {
            double Re = myPipe.v * myPipe.get_inner_diameter() / current_layer[1][j];
            double resistance = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            current_layer[2][j] = p_0;
            double p_rachet = p_0 + myPipe.get_dx() * (resistance / myPipe.get_inner_diameter() * current_layer[0][j - 1] * pow(myPipe.v, 2) / 2 - M_G * current_layer[0][j - 1] * (myPipe.z_l - myPipe.z_0) / ((myPipe.n - 1) * myPipe.get_dx()));
            p_0 = p_rachet;

            outFile << i * myPipe.get_dt() << "," << j * myPipe.get_dx() << "," << current_layer[0][j] << "," << current_layer[1][j] << "," << current_layer[2][j] << "\n";
        }
        outFile.close();
    }
    else {
        ofstream outFile("block_3.csv", ios::app);
        // Записать значения текущего слоя в файл
        for (size_t j = 1; j < current_layer[0].size(); j++) {
            double Re = myPipe.v * myPipe.get_inner_diameter() / current_layer[1][j];
            double resistance = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            current_layer[2][j] = p_0;
            double p_rachet = p_0 + myPipe.get_dx() * (resistance / myPipe.get_inner_diameter() * current_layer[0][j - 1] * pow(myPipe.v, 2) / 2 - M_G * current_layer[0][j - 1] * (myPipe.z_l - myPipe.z_0) / ((myPipe.n - 1) * myPipe.get_dx()));
            p_0 = p_rachet;

            outFile << i * myPipe.get_dt() << "," << j * myPipe.get_dx() << "," << current_layer[0][j] << "," << current_layer[1][j] << "," << current_layer[2][j] << "\n";
        }
        outFile.close();
    }
}

int main() {
    pipe myPipe;
    myPipe.L = 100000;
    myPipe.p_0 = 6e6;
    myPipe.d_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_l = 50;
    myPipe.v = 0.5;
    myPipe.n = 10;
    myPipe.sher = 15e-6;

    massiv ro;
    ro.massiv = { 800 };

    massiv nu;
    nu.massiv = { 10e-6 };

    massiv time;
    time.massiv = { 0 };

    vector<double> ro_start(myPipe.n, 900);
    vector<double> nu_start(myPipe.n, 15e-6);
    vector<double> time_start(myPipe.n, 0);

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_start, nu_start, time_start });

    for (size_t h = 0; h < myPipe.n; h++) {

        characteristic_method(myPipe, ro.massiv[0], buffer.current()[0], buffer.previous()[0]);
        characteristic_method(myPipe, nu.massiv[0], buffer.current()[1], buffer.previous()[1]);
        out_put(myPipe, buffer, h, time);
        buffer.advance(1);

    }
}






