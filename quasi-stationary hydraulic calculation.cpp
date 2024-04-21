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

    // кол-во сетки
    double get_dx() const {
        return L / n;
    }

    //Шаг по времени
    double get_dt() const {
        return get_dx() / v;
    }

    // шаг по координате
    double get_n() const {
        return (L / v) / get_dt();
    }
};

/// @brief Массив данных
struct massiv {
    vector<double> massiv;
};

/// @brief Ввод значений начальных условий
/// @param myPipe Ссылка на структуру начальных условий
void iniFun(pipe& myPipe) {
    myPipe.L = 100000;
    myPipe.p_0 = 6e6;
    myPipe.d_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_l = 50;
    myPipe.v = 2;
    myPipe.n = 10;
    myPipe.sher = 15e-6;
}

/// @brief инициализаци данных масивов 
/// @param myPipe ссылка на данные о трубопроводе
/// @param ro плотность
/// @param nu вязкость
/// @param time время моделирования
/// @param ro_start начальная плотность в трубе
/// @param nu_start начальная вязкость в трубе 
/// @param time_start время моделирования
void initializeVariables(pipe& myPipe, massiv& ro, massiv& nu, massiv& time, vector<double>& ro_start, vector<double>& nu_start, vector<double>& time_start) {
    ro.massiv = { 800 };
    nu.massiv = { 10e-6 };
    time.massiv = { 0 };

    ro_start = vector<double>(myPipe.n, 900);
    nu_start = vector<double>(myPipe.n, 15e-6);
    time_start = vector<double>(myPipe.n, 0);
}

/// @brief отдельная функция для метода характеристик
/// @param iniPipe объявление переменной iniPipe (Данные о параметрах трубопровода)
/// @param argument  нач.значения
/// @param current_layer вектор текущий слой 
/// @param previous_layer вектор предыущий слой
class CharacteristicMethod {
public:
    /// @brief функция расчета методом характеристик
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param current_layer текущий слой
    /// @param previous_layer предыдущий слой
    void apply(pipe& myPipe, double argument, vector<double>& current_layer, vector<double>& previous_layer) {
        for (size_t i = 1; i < current_layer.size(); i++)
            current_layer[i] = previous_layer[i - 1];
        current_layer[0] = argument;
    }
};

/// @brief Класс метод Эйлера + запись в файл после каждой итерации
class FileWriter {
private:
    std::ofstream outFile;
    std::string filename;
public:
    /// @brief функция записи результатов расчета в файл
    /// @param filename название файла для записи
    FileWriter(const std::string& filename) : filename(filename) {
        outFile.open(filename, std::ios::app);
        if (!outFile.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << std::endl;
            exit(1);
        }
    }
    ~FileWriter() {
        if (outFile.is_open()) {
            outFile.close();
        }
    }
    /// @brief функция расчета методом Эйлера
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param current_layer текущий слой
    /// @param i счетчик итераций 
    /// @param p_0 давление в начале трубопровода
    void writeData(const pipe& myPipe, std::vector<std::vector<double>>& current_layer, int i, double p_0) {
        // Создаем копию current_layer, чтобы изменять ее содержимое
        std::vector<std::vector<double>> modified_layer = current_layer;

        for (size_t j = 1; j < modified_layer[0].size(); j++) {
            double Re = myPipe.v * myPipe.get_inner_diameter() / current_layer[1][j];
            double resistance = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            current_layer[2][j] = p_0;
            modified_layer[2][j] = p_0; // Теперь изменяем копию, а не оригинал
            double p_rachet = p_0 + myPipe.get_dx() * ((-resistance) / myPipe.get_inner_diameter() * modified_layer[0][j - 1] * pow(myPipe.v, 2) / 2 - M_G * modified_layer[0][j - 1] * (myPipe.z_l - myPipe.z_0) / ((myPipe.n - 1) * myPipe.get_dx()));
            p_0 = p_rachet;

            outFile << i * myPipe.get_dt() << "," << j * myPipe.get_dx() << "," << modified_layer[0][j] << "," << modified_layer[1][j] << "," << modified_layer[2][j] << "\n";
        }
        //  Копируем измененные данные обратно в current_layer
        current_layer = modified_layer;
    }
};

/// @brief функция вывода значений результов расчета 
/// @param myPipe ссылка на данные о трубопроводе
/// @param buffer буффер данных
/// @param i счетчик итераций 
/// @param time время моделирования
void out_put(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, int i, massiv time) {
    vector<vector<double>>& current_layer = buffer.current();
    double p_0 = myPipe.p_0;

    if (i == 0) {
        ofstream outFile("block_3.csv");
        outFile << "Время,Координата,Плотность,Вязкость,Давление" << "\n";
        outFile.close();
    }

    FileWriter writer("block_3.csv");
    writer.writeData(myPipe, current_layer, i, p_0);
}

/// @brief класс расчета методом характеристик
class PipeProcessor {
private:
    CharacteristicMethod characteristicMethod;
public:
    /// @brief функция смещения предыдущего слоя и запись граничного условия
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param buffer буффер данных
    /// @param ro плоность
    /// @param nu вязкость
    /// @param time время моделирования
    void process(pipe& myPipe, ring_buffer_t<vector<vector<double>>>& buffer, massiv& ro, massiv& nu, massiv& time) {
        for (size_t h = 0; h < myPipe.n; h++) {
            for (size_t j = 0; j < buffer.current().size(); j++) {
                characteristicMethod.apply(myPipe, ro.massiv[0], buffer.current()[0], buffer.previous()[0]);
                characteristicMethod.apply(myPipe, nu.massiv[0], buffer.current()[1], buffer.previous()[1]);
                out_put(myPipe, buffer, h, time);
            }
            buffer.advance(1);
        }
    }
};

int main() {
    pipe myPipe;
    iniFun(myPipe);

    massiv ro;
    massiv nu;
    massiv time;

    vector<double> ro_start;
    vector<double> nu_start;
    vector<double> time_start;

    initializeVariables(myPipe, ro, nu, time, ro_start, nu_start, time_start);

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_start, nu_start, time_start });

    PipeProcessor pipeProcessor;
    pipeProcessor.process(myPipe, buffer, ro, nu, time);
}






