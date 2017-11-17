#ifndef SIMPLEX_H_
#define SIMPLEX_H_
#include <vector>
#include <string>
void show_task(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>>& a);
void show_new_task(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>>& a);
void set_data(std::vector<double> &c, std::vector<double> &b, std::vector<std::vector<double>>& a, unsigned int &n, unsigned int &m);
void simplex_method(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>>& a);
bool not_reference(double ** arr, unsigned int m);
bool not_optimal(double ** arr, unsigned int m, unsigned int n);
std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n);
void set_new_data(std::vector<double> &c, std::vector<double> &b, std::vector<std::vector<double>>& a);
void calculation_of_criteria(std::vector<std::vector<double>> a);
static unsigned int ri = {0}; // Разрешающая строка
static unsigned int rj = {0}; // Разрешающий столбец
static std::string min_max; // Ищем min/max
static unsigned int chk = {0};
static unsigned int count = {0};
static std::vector<std::string> strategy; // Все стратегии с соответствующими Ui и Vj
static std::vector<double> x_strategy; // Стратегии числовые для игрока A 
static std::vector<double> y_strategy; // Стратегии числовые для игрока B
static std::vector<unsigned int> index; // Вектор индексов переменных
static double g = {0.0}; // Минимальный выигрыш игрока A
static double h = {0.0}; // Максимальный проигрыш игрока B
const static double alpha = {0.5}; // Числовой параметр для критерия Гурвица
#endif