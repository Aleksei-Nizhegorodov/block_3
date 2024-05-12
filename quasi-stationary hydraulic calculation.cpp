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
    double roughness;   //абсюлютная шероховатость, [м]
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

    /// @brief Внутренний диаметр трубопровода
    /// @return 
    double get_inner_diameter() const {
        return d_vnesh - 2 * b;
    }

    /// @brief Относительная шероховатость
    /// @return 
    double get_relative_roughness() const {
        return roughness / get_inner_diameter();
    }

    /// @brief 
    /// @return 
    double get_inner_area() const {
        double D = get_inner_diameter();
        double S = M_PI * D * D / 4;
        return v * S;
    }

    /// @brief Скорость жидкости
    /// @return 
    double get_v() const {
        double D = get_inner_diameter();
        return 4 * Q / (3.1415 * pow(D, 2));
    }

    /// @brief Число Рейнольдса
    /// @return 
    double get_Re() const {
        double D = get_inner_diameter();
        return get_v() * D / nu;
    }

    /// @brief Касательное напряжение трения
    /// @return 
    double get_t_w() const {
        return resistance / 8 * density * pow(get_v(), 2);
    }

    /// @brief кол-во сетки
    /// @return 
    double get_dx() const {
        return L / n;
    }
        
    /// @brief Шаг по времени
    /// @return 
    double get_dt() const {
        return get_dx() / v;
    }

    /// @brief шаг по координате
    /// @return 
    double get_n() const {
        return (L / v) / get_dt();
    }
};

/// @brief Массив данных
struct massiv {
    vector<double> ro;
    vector<double> nu;
    vector<double> time;
};

/// @brief Ввод значений начальных условий
/// @param myPipe Ссылка на структуру начальных условий
void Pipe_1(pipe& myPipe) {
    myPipe.L = 100000;
    myPipe.p_0 = 6e6;
    myPipe.d_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_l = 50;
    myPipe.v = 2;
    myPipe.n = 10;
    myPipe.roughness = 15e-6;
}

/// @brief инициализаци данных масивов 
/// @param myPipe ссылка на данные о трубопроводе
/// @param ro плотность
/// @param nu вязкость
/// @param time время моделирования
/// @param ro_start начальная плотность в трубе
/// @param nu_start начальная вязкость в трубе 
/// @param time_start время моделирования
void initializeVariables(pipe& myPipe, massiv& ro, massiv& nu, massiv& time, vector<double>& ro_start, vector<double>& nu_start, vector<double>& time_start, double& ro_pulsing, vector<double>& nu_pulsing, vector<double>& time_pulsing) {
    
    ro.ro = vector<double>(myPipe.n, 800);
    nu.nu = vector<double>(myPipe.n, 10e-6);
    time.time = vector<double>(myPipe.n, 0); 

    ro_start = vector<double>(myPipe.n, 900);
    nu_start = vector<double>(myPipe.n, 15e-6);
    time_start = vector<double>(myPipe.n, 0);

    ro_pulsing =  990;
    nu_pulsing = 19e-6;
    time_pulsing = 0;
}

/// @brief Класс метода характеристик
class CharacteristicMethod {
public:
    /// @brief функция расчета методом характеристик
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param current_layer текущий слой
    /// @param previous_layer предыдущий слой
    void Characteristic(pipe& myPipe, double argument, vector<double>& current_layer, vector<double>& previous_layer) {
        for (size_t i = 1; i < current_layer.size(); i++)
            current_layer[i] = previous_layer[i - 1];
        current_layer[0] = argument;
    }
};

/// @brief Класс метод Эйлера
class EulerMethod {
public:
    /// @brief функция расчета методом Эйлера
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param current_layer ссылка на текущий слой
    /// @param p_0 начальное давление
    void Euler(const pipe& myPipe, std::vector<std::vector<double>>& current_layer, double p_0) {
        // Создаем копию current_layer, чтобы изменять ее содержимое
        std::vector<std::vector<double>> modified_layer = current_layer;

        for (size_t j = 1; j < modified_layer[0].size(); j++) {
            double Re = myPipe.v * myPipe.get_inner_diameter() / current_layer[1][j];
            double resistance = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            current_layer[2][j] = p_0;
            modified_layer[2][j] = p_0; // Теперь изменяем копию, а не оригинал
            double p_rachet = p_0 + myPipe.get_dx() * ((-resistance) / myPipe.get_inner_diameter() * modified_layer[0][j - 1] * pow(myPipe.v, 2) / 2 - M_G * modified_layer[0][j - 1] * (myPipe.z_l - myPipe.z_0) / ((myPipe.n - 1) * myPipe.get_dx()));
            p_0 = p_rachet;

        }
        //  Копируем измененные данные обратно в current_layer
        current_layer = modified_layer;
    }
};

/// @brief Класс записи в файл
class FileWriter {
private:
    std::ofstream outFile;
    std::string filename;
public:
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
    /// @brief функция записи в файл
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param buffer буфер о хранении данных о текущем слое
    /// @param time время моделирования
    void out_put(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, double& time) {
        
        if (time == 0) {
            ofstream outFile("block_3.csv");
            outFile << "Время,Координата,Плотность,Вязкость,Давление" << "\n";
            outFile.close();
        }
        else {

            vector<vector<double>>& current_layer = buffer.current();

            for (size_t j = 1; j < current_layer[0].size(); j++) {
                outFile << time << "," << j * myPipe.get_dx() << "," << current_layer[0][j] << "," << current_layer[1][j] << "," << current_layer[2][j] << "\n";
            }
        }
    }
};


/// @brief класс расчета методом характеристик
class PipeProcessor {
private:
    CharacteristicMethod characteristicMethod;
    FileWriter fileWriter; // Поле для хранения экземпляра FileWriter

public:
    PipeProcessor(const std::string& filename) : fileWriter(filename) {} // Конструктор, инициализирующий fileWriter
    /// @brief функция замены слоев 
    /// @param myPipe ссылка на данные о трубопроводе
    /// @param buffer буфер храннеия данных о слое
    /// @param ro плотность
    /// @param nu вязкость
    /// @param time время моделирования
    void process(pipe& myPipe, ring_buffer_t<vector<vector<double>>>& buffer, massiv& ro, massiv& nu, double& ro_start, double& nu_start, double& time_start, double& ro_pulsing, double& nu_pulsing, double& time_pulsing) {
        
         
        double time = 0;
        double total_time = 12 * 3600; 
        double pulsingStart = 5 * 3600;
        double pulsingEnd = 10 * 3600;
       
        while (time < total_time)
        {

            time += myPipe.get_dt();
            if (time >= pulsingStart && time <= pulsingEnd) {
                characteristicMethod.Characteristic(myPipe, ro_pulsing, buffer.current()[0], buffer.previous()[0]);
                characteristicMethod.Characteristic(myPipe, nu_pulsing, buffer.current()[1], buffer.previous()[1]);
            }

            else {
                characteristicMethod.Characteristic(myPipe, ro_start, buffer.current()[0], buffer.previous()[0]);
                characteristicMethod.Characteristic(myPipe, nu_start, buffer.current()[1], buffer.previous()[1]);
            }

            vector<vector<double>>& current_layer = buffer.current();
            double p_0 = myPipe.p_0;
           
            // Использование метода Эйлера
            EulerMethod eulerMethod;
            eulerMethod.Euler(myPipe, current_layer, p_0);

            fileWriter.out_put(myPipe, buffer,time);

            buffer.advance(1);

        }
                          
    }
};

int main() {
    pipe myPipe;
    Pipe_1(myPipe);

    massiv ro;
    massiv nu;
    massiv time;

    vector<double> ro_start;
    vector<double> nu_start;
    vector<double> time_start;

    vector<double> ro_pulsing;
    vector<double> nu_pulsing;
    vector<double> time_pulsing;

    initializeVariables(myPipe, ro, nu, time, ro_start, nu_start, time_start, ro_pulsing, nu_pulsing, time_pulsing);

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_start, nu_start, time_start });

    PipeProcessor pipeProcessor("block_3.csv");
    pipeProcessor.process(myPipe, buffer, ro, nu, time, ro_start, nu_start, time_start, ro_pulsing, nu_pulsing, time_pulsing);

    return 0;
}