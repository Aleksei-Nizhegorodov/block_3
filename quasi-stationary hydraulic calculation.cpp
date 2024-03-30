#include <iomanip>
#include <iostream>
#include <clocale>
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

/// @brief Подключение пространств имен
using namespace std;
using namespace pde_solvers;

/// @brief Данные о трубопроводе
struct pipe {
    double L;           //длина, [м]
    double d_vnesh;     //внешний диаметр, [мм]
    double b;           //толщина стенки, [мм]
    double sher;        //абсюлютная шероховатость, [м]
    double z_0;         //высотная отметка в начале трубопровода, [м]
    double z_l;         //высотная отметка в конце трубопровода, [м]
    double n;           //кол-во точек расчетной сетки
    double density;     //плотность, [кг/м^3]
    double nu;          //кинематическая вязкость, [сСт]
    double p_0;         //давление в начале трубопровода, [МПа]
    double v;           //скорость течения жидкости, [м/с]
    double Q;           //объемный расход, [м^3/с]
    double resistance;  //коэфф.гидр.сопротивления
    double t_w;         //касательное напряжение трения

    //Внутренний диаметр трубопровода
    double get_inner_diameter() const {
        return d_vnesh - 2 * b;
    }

    //Относительная шероховатость
    double get_relative_roughness() const {
        return sher / get_inner_diameter();
    }

    //
    double get_inner_area() const {
        double D = get_inner_diameter();
        double S = M_PI * D * D / 4;
        return v * S;
    }

    //Скорость жидкости
    double get_v() const {
        double D = get_inner_diameter();
        return 4 * Q / (3.1415 * pow(D, 2));
    }

    //Число Рейнольдса
    double get_Re() const {
        double D = get_inner_diameter();
        return get_v() * D / nu;
    }

    //Касательное напряжение трения
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

/// @brief Массив данных
struct massiv {
    vector<double> massiv;
};

/// @brief Фун-ия расчета методом характеристик
/// @param myPipe ссылка на данные о трубопроводе
/// @param parametr 
/// @param current_layer текущий слой
/// @param previous_layer предыдущий слой
void characteristic_method(pipe myPipe, double parametr, vector<double>& current_layer, vector<double>& previous_layer) {
    current_layer[0] = parametr;
    for (size_t j = 1; j < myPipe.n; j++) {
        for (size_t i = 1; i < myPipe.n; i++) {
            current_layer[i] = previous_layer[i - 1];
        }
    }
}

/// @brief Фун-ия вывода данных расчета в excel формат
/// @param myPipe ссылка на данные о трубопроводе
/// @param buffer 
/// @param i 
/// @param time 
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
    myPipe.n = 100;
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






