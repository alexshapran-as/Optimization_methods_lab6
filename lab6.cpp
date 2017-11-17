#include <iostream>
#include "simplex.h"
#include <vector>
#include <string>

int main()
{
	using std::cout;
	using std::cin;
	using std::endl;
	using std::vector;
	system("chcp 65001");

	vector<double> c;
	vector<double> b;
	vector<vector<double>> a;
	unsigned int n = {0}; // количество столбцов матрицы A
	unsigned int m = {0}; // Количество строк матрицы A

	set_data(c,b,a,n,m); // Берем данные о задаче
	simplex_method(c,b,a); // Решаем ПЗ ЛП с помощью симплекс метода 
	set_new_data(c,b,a); // Меняем данные ПЗ ЛП на данные ДЗ ЛП
	simplex_method(c,b,a); // Решаем ДЗ ЛП с помощью симплекс метода 
	calculation_of_criteria(a); // Вычисляем критерии игры

	
	cout << endl << endl;
	system("pause");
	return 0;
}

