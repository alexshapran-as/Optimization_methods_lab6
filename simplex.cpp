#include <iostream>
#include <fstream>
#include "simplex.h"
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

void set_data(std::vector<double> &c, std::vector<double> &b, std::vector<std::vector<double>>& a, unsigned int &n, unsigned int &m)
{
	using std::cout;
	using std::endl;
	using std::vector;

	double temp = {0.0};
	std::ifstream fin("lab6.txt");
	if (!fin.is_open())
	{
		cout << endl << "Ошибка открытия файла.";
		exit(EXIT_FAILURE);
	}
	fin >> min_max;
	if (min_max != "min" && min_max != "max")
	{
		cout << endl << "Неправильно указано, к чему стремится ЦФ W.";
		exit(EXIT_FAILURE);
	}
	while (fin >> temp)
	{
		if (min_max == "min")
			c.push_back(temp);
		else if (min_max == "max")
			c.push_back(-temp);
		n++;
	}
	fin.clear();
	while (fin.get() != '\n');
	while (fin >> temp)
	{
		b.push_back(-temp);
		m++;
	}
	fin.clear();
	while (fin.get() != '\n');

	if (m > 2*n)
	{
		cout << endl << "Число линейно независимых уравнений больше числа переменных."
		" Такая СЛАУ несовместна.";
		exit(EXIT_FAILURE);
	}

	a.resize(m);
	for (unsigned int i = {0}; i < m; ++i)
		a[i].resize(n);
	for (unsigned int i = {0}; i < m; ++i)
	{
		for (unsigned int j = {0}; j < n; ++j)
		{
			if (!(fin >> a[i][j]))
			{
				cout << endl << "Не удалось считать данные для вектора A = [a1,..,an]!";
				exit(EXIT_FAILURE);
			}
			a[i][j] = -a[i][j];
		}
	}
	fin.close();
	show_task(c,b,a);
}

void show_task(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>>& a)
{
	using std::cout;
	using std::cin;
	using std::endl;

	cout << endl << "Постановка задачи ЛП на основе матрицы стратегий игроков. Требуется найти решение следующей задачи:" 
	<< endl << std::setw(26) << "W = 1*u -> " <<  min_max << endl 
	<< std::setw(25) << " C(i,j)u >= 1" << endl << std::setw(25) << "u >= 0" << endl
	<< "Каноническая форма задачи ЛП:" << endl;
	if (min_max == "min")
		cout << std::setw(17) << "W = ";
	else
		cout << std::setw(17) << "-W = ";
	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		if (c[i] < 0)
			cout << c[i] << "u" << i+1;
		else
			cout << "+" << c[i] << "u" << i+1;
	}
	cout << " -> ";
	if (min_max == "min")
		cout << min_max << endl << endl; 
	else
		cout << "min" << endl << endl; 
	unsigned int k = c.size() + 1;
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		cout << "                 ";
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			if (a[i][j] != 0)
			{
				if (a[i][j] != 1 && a[i][j] != -1)
				{
					if (a[i][j] < 0 || j == 0)
						cout << a[i][j] << "u" << j+1;
					else
						if (a[i][j-1] != 0)
							cout << "+" << a[i][j] << "u" << j+1;
						else 
							cout << a[i][j] << "u" << j+1;
				}
				else
				{
					if (a[i][j] < 0)
						cout << "-u" << j+1;
					else if (j != 0)
						cout << "+" << "u" << j+1;
					else 
						cout << "u" << j+1;
				}
		
			}		
		}
		cout << "+u" << k++ << " = " << b[i] << endl;
	}
	cout << endl << "             ";
	for (unsigned int i = {0}; i < c.size()+b.size(); ++i)
	{
		if (i != c.size()+b.size()-1)
			cout << "u" << i+1 << ", ";
		else
			cout << "u" << i+1;
	}
	cout << " >= 0" << endl << "Алгоритм преобразования симплекс-таблицы (жордановы исключения):";

}
 
void  simplex_method(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>>& a)
{
	using std::cout;
	using std::endl;
	using std::vector;
	using std::string;

	vector<string> x(b.size() + c.size() + 1);

	double ** arr = new double * [b.size()+1];  // Двумерный массив - симплекс матрица
	for (unsigned int i = {0}; i < b.size()+1; ++i)
		arr[i] = new double [c.size()+1];
	for (unsigned int i = {0}; i < b.size(); ++i)
		arr[i][0] = b[i]; // Заполняем первый столбец свободных членов
	arr[b.size()][0] = {0}; // Свободный элемент в функции F
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		unsigned int k = {1};
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			arr[i][k] = a[i][j]; // Заполняем все остальные элементы, начиная со второго столбца, без последней строки
			k++;
		}
	} 
	unsigned int k = {1};
	for (unsigned  int i = {0}; i < c.size(); ++i)
	{
		arr[b.size()][k] = -c[i]; // Заполняем последнюю строку со второго столбца коэфф-ми при x в ф-ии F
		k++;
	}

	x = print(arr, b.size(), c.size());

	while (not_reference(arr,b.size()))
	{
		cout << endl << "Недопустимое решение. Производим замену базиса," 
		" чтобы получить опорное решение.";
		unsigned int i = {0};
		while (arr[i][0] > 0 && i < b.size() + 1)
			i++;
		unsigned int j = {1}; // Разрешающий столбец
		unsigned int k = {0};
		for (; j < c.size() + 1; ++j)
		{
			if (arr[i][j] < 0)
			{
				k++;
				break;
			}
		}
		if (k == 0)
		{
			cout << endl << "Задача не имеет допустимых решений!";
			exit(EXIT_FAILURE);
		}
		vector<double> temp;
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			if (arr[i][0] != 0 || arr[i][j] > 0)
				temp.push_back(arr[i][0] / arr[i][j]); // Заполняем временный вектор temp отношениями свободных членов к элементам разрешающего столбца
		}
		std::sort(temp.begin(),temp.end()); // Сортируем эти отношения
		unsigned int m = {0};
		while (temp[m] < 0)
			m++; // Ищем номер первого положительного элемента в отсортированном векторе temp 
		i = {0}; // Разрешающая строка
		if (arr[ri][rj] == 0)
		{
			cout << endl << "Деление на 0 недопустимо!";
			exit(EXIT_FAILURE);
		}
		for (; temp[m] != arr[i][0] / arr[i][j]; ++i);
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n != j) // Преобразование эл-та, который не стоит ни в разрешающей строке, ни в разрешающем столбце
				{
					arr[m][n] = arr[m][n] - (arr[m][j] * arr[i][n]) / (arr[i][j]);
				}
			}
		}
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n == j) // Преобразование эл-та, который стоит в разрешающем столбце, но не в разрешающей строке
				{
					if (arr[m][n] != 0)
						arr[m][n] = -arr[m][n] / arr[i][j];
				}
				if (m == i && n != j) // Преобразование эл-тва, который стоит в разрешающей строке, но не в разрешающем столбце
				{
					arr[m][n] = arr[m][n] / arr[i][j];
				}
			}
		}
		arr[i][j] = 1.0 / arr[i][j]; // Преобразование разрешающего эл-та
		ri = i; 
		rj = j;
		x = print(arr, b.size(), c.size());
	}

	cout << endl << "Найдено опорное решение. Проводим проверку на оптимальность.";
	while (not_optimal(arr,b.size(),c.size()))
	{
		cout << endl << "Неоптимальное решение. Производим замену базиса.";
		unsigned int j = {1}; // Новый разрешающий столбец 
		for (; j < c.size() + 1; ++j)
		{
			if (arr[b.size()][j] > 0) // Ищем разрешающий столбец
				break;
		}
		vector<double> temp;
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			if (arr[i][0] != 0 || arr[i][j] > 0)
				temp.push_back(arr[i][0] / arr[i][j]); 
		}
		std::sort(temp.begin(),temp.end());
		unsigned int m = {0};
		while (temp[m] < 0) // Аналогично ищем номер первого положительного отношения
			m++;
		unsigned int i = {0}; // Новая разрешающая строка
		if (arr[ri][rj] == 0)
		{
			cout << endl << "Деление на 0 недопустимо!";
			exit(EXIT_FAILURE);
		}
		for (; temp[m] != arr[i][0] / arr[i][j]; ++i);
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n != j)
				{
					arr[m][n] = arr[m][n] - (arr[m][j] * arr[i][n]) / (arr[i][j]);
				}
			}
		}
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n == j)
				{
					if (arr[m][n] != 0)
						arr[m][n] = -arr[m][n] / arr[i][j];
				}
				if (m == i && n != j)
				{
					arr[m][n] = arr[m][n] / arr[i][j];
				}
			}
		}
		arr[i][j] = 1.0 / arr[i][j];
		ri = i;
		rj = j;
		x = print(arr, b.size(), c.size()); // Вектор строк, в котором первыми идут свободные переменные, а затем идут базисные переменные
	}
	cout << endl << "Найдено оптимальное решение!" << endl;
	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		cout << x[i] << " = 0; ";
		strategy.push_back(x[i]); // Добавляем все Ui или Vj в стратегии, чтобы выбрать затем интересующие
		strategy.push_back("0"); // Числовое значение стратегии
	}
	unsigned int l = {0};
	for (unsigned int j = c.size(); j < b.size()+c.size(); ++j)
	{
		cout << endl << x[j] << " = " << arr[l][0] << ";";
		strategy.push_back(x[j]);
		strategy.push_back(std::to_string(arr[l][0]));
		l++;
	}
	cout << endl;
	double result = {0.0};
	double sum = {0.0};
	if (min_max == "min") 
	{
		if (x[0].find("U") != -1)
		{
			cout << "W = " << arr[b.size()][0];
			result = arr[b.size()][0];
			strategy.push_back(std::to_string(result)); // Добавляем значение W
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
		else
		{
			std::swap(arr[1][0],arr[2][0]);
			cout << "-Z = " << -arr[b.size()][0];
			result = -arr[b.size()][0];
			strategy.push_back(std::to_string(result)); // Добавляем значение Z
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
	}
	else if (min_max == "max")
	{
		if (x[0].find("U") != -1)
		{
			cout << "Z = " << -arr[b.size()][0]; // Инвертируем знак W, так как ищем максимум
			result = -arr[b.size()][0];
			strategy.push_back(std::to_string(result));
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
		else
		{
			std::swap(arr[1][0],arr[2][0]);
			cout << "Z = " << arr[b.size()][0];
			result = arr[b.size()][0];
			strategy.push_back(std::to_string(result));
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
	}


	if (strategy[0].find("U") != -1)
	{
		g = 1. / (atof(strategy[strategy.size() - 1].c_str()));
		vector<string> tmp;
		for (unsigned int j = {1}; j <= c.size(); ++j)
		{
			tmp.push_back("U" + std::to_string(j));
		}
		for (unsigned int i = {0}; i < strategy.size(); ++i)
		{
			for (auto it : tmp)
			{
				if (strategy[i].find(it) != -1)
				{
					x_strategy.push_back(atof(strategy[i + 1].c_str()) * g); // Считаем интересующие стратегии в числовом виде
					strategy[i].erase(strategy[i].begin());
					index.push_back(atoi(strategy[i].c_str())); // Сохраняем индекс стратегии
				}
			}
		}
		cout << endl <<  std::setw(20) << "g = 1 / W = " << g << " - минимальный выйгрыш игрока A" << endl;
		cout << endl << "    " << "Оптимальная смешанная стратегия игрока A: " << endl;
		for (unsigned int i = {0}; i < x_strategy.size(); ++i)
		{
			cout << "       x" << index[i] << " = " << "U" << index[i] << " * g = " << x_strategy[i] << ";" << endl;
		}
		cout << endl <<  std::setw(25) << "Проверка: " << endl << "       "; // Проверка на то, что сумма вероятностей == 1
		double sum = {0.0};
		for (unsigned int i = {0}; i < x_strategy.size(); ++i)
		{
			cout << "x" << i + 1;
			if (i == x_strategy.size() - 1) cout << " = ";
			else cout << " + ";
			sum += x_strategy[i];
		}
		cout << sum;
	}
	index.clear();

	if (strategy[0].find("V") != -1)
	{
		h = 1. / (atof(strategy[strategy.size() - 1].c_str()));
		vector<string> tmp;
		for (unsigned int j = {1}; j <= c.size(); ++j)
		{
			tmp.push_back("V" + std::to_string(j));
		}
		for (unsigned int i = {0}; i < strategy.size(); ++i)
		{
			for (auto it : tmp)
			{
				if (strategy[i].find(it) != -1)
				{
					y_strategy.push_back(atof(strategy[i + 1].c_str()) * h); // Считаем интересующие стратегии в числовом виде
					strategy[i].erase(strategy[i].begin());
					index.push_back(atoi(strategy[i].c_str())); // Сохраняем индекс стратегии
				}
			}
		}
		cout << endl <<  std::setw(20) << "h = 1 / Z = " << h << " - максимальный проигрыш игрока B" << endl;
		cout << endl << "    " << "Оптимальная смешанная стратегия игрока B: " << endl;
		for (unsigned int i = {0}; i < y_strategy.size(); ++i)
		{
			cout << "       y" << index[i] << " = " << "V" << index[i] << " * h = " << y_strategy[i] << ";" << endl;
		}
		cout << endl <<  std::setw(25) << "Проверка: " << endl << "       "; // Проверка на то, что сумма вероятностей == 1
		double sum = {0.0};
		for (unsigned int i = {0}; i < y_strategy.size(); ++i)
		{
			cout << "y" << i + 1;
			if (i == y_strategy.size() - 1) cout << " = ";
			else cout << " + ";
			sum += y_strategy[i];
		}
		cout << sum;

	}

	strategy.clear();
	for (unsigned int i = {0}; i < b.size() + 1; ++i)
		delete [] arr[i];
	delete [] arr;
}

void calculation_of_criteria(std::vector<std::vector<double>> a)
{
	using std::cout;
	using std::endl;
	using std::setw;
	std::multimap<double, unsigned int> sum;
	std::vector<unsigned int> a_strategy;

	cout << endl << endl << "		Критерий Бернулли:" << endl;
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		double temp_sum = {0.0};
		for (unsigned int j = {0}; j < y_strategy.size(); ++j)
		{
			temp_sum += a[i][j];
		}
		sum.insert(std::pair<double, unsigned int>((1. / y_strategy.size()) * temp_sum, i + 1));
	}
	cout << endl;
	
	for (auto it = sum.begin(); it != sum.end(); ++it)
	{
		cout << setw(12) << "a" << it->second << ": (1 / n)*Sum(aij) = " << it->first << endl;
	}

	cout << endl << "Таким образом, необходимо руководствоваться стратегией: "
	<< endl <<  "a" << (--sum.end())->second << ", у которой соответствующее математическое"
	" ожидание выигрыша максимально и равно " << (--sum.end())->first << endl;
	a_strategy.push_back((--sum.end())->second);

	cout << endl << endl << "		Критерий Вальда:" << endl;
	std::vector<std::vector<double>> temp_array;
	std::multimap<double, unsigned int> mins;
	temp_array = a;
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		std::sort(temp_array[i].begin(), temp_array[i].end());
		mins.insert(std::pair<double, unsigned int>(temp_array[i][0], i + 1));
	}
	cout << endl << "Пессимистическая стратегия (критерий Вальда) определяет выбор a" 
	<< (--mins.end())->second << " (Нижняя цена игры равна " << (--mins.end())->first 
	<< ")." << endl;
	a_strategy.push_back((--mins.end())->second);

	cout << endl << endl << "		Критерий максимума:" << endl;
	std::multimap<double, unsigned int> maxs;
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		maxs.insert(std::pair<double, unsigned int>(temp_array[i][temp_array[i].size() - 1], i + 1));
	}
	cout << endl << "Оптимистическая стратегия соответствует выбору a" 
	<< (--maxs.end())->second << " (максимально возможный выигрыш " << (--maxs.end())->first 
	<< ")." << endl;
	a_strategy.push_back((--maxs.end())->second);

	cout << endl << endl << "		Критерий Гурвица:" << endl;
	std::multimap<double, unsigned int> temp_maxs;
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		temp_maxs.insert(std::pair<double, unsigned int>(alpha*temp_array[i][0] + (1 - alpha)*temp_array[i][temp_array[i].size() - 1], i + 1));
	}
	cout << endl << "По этому критерию наилучшая стратегия — a" << (--temp_maxs.end())->second
	<< " (ожидаемый выигрыш равен " << (--temp_maxs.end())->first << ").";
	a_strategy.push_back((--temp_maxs.end())->second);

	cout << endl << endl << "		Критерий Сэвиджа:" << endl;
	std::vector<std::vector<double>> r_array;
	std::vector<std::vector<double>> temp_r_array;
	r_array = a;
	temp_r_array.resize(y_strategy.size());
	for (unsigned int i = {0}; i < y_strategy.size(); ++i)
	{
		temp_r_array[i].resize(x_strategy.size());
	}
	for (unsigned int i = {0}; i < y_strategy.size(); ++i)
	{
		for (unsigned int j = {0}; j < x_strategy.size(); ++j)
		{
			temp_r_array[i][j] =  r_array[j][i];
		}
	}
	for (unsigned int i = {0}; i < y_strategy.size(); ++i)
	{
		std::sort(temp_r_array[i].begin(), temp_r_array[i].end());
	}
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		for (unsigned int j = {0}; j < y_strategy.size(); ++j)
		{
			r_array[i][j] = temp_r_array[j][temp_r_array[j].size() - 1] - r_array[i][j];
		}
	}
	cout << endl << "Таблица рисков имеет вид:" << endl << endl;
	cout << "Стратегии";
	for (unsigned int j = {0}; j < y_strategy.size(); ++j)
	{
		cout << setw(6) << "b" << j + 1 << "     ";
	}
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		cout << endl << setw(4) << "a" << i + 1;
		for (unsigned int j = {0}; j < y_strategy.size(); ++j)
		{
			cout << setw(12) << r_array[i][j];
		}
	}
	mins.clear();
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		std::sort(r_array[i].begin(), r_array[i].end());
		mins.insert(std::pair<double, unsigned int>(r_array[i][r_array[i].size() - 1], i + 1));
	}
	cout << endl << endl << "Таким образом, оптимальная рисковая стратегия — a" << (mins.begin())->second
	<< " (так как для этой стратегии min_i(max_j(max_i(Cij) - Cij))) = " << (mins.begin())->first 
	<< ")." << endl;
	a_strategy.push_back((mins.begin())->second);
	std::sort(a_strategy.begin(), a_strategy.end());
	std::map<unsigned int, unsigned int> a_strg_count;
	unsigned int j = {0};
	auto it = a_strategy.begin();
	for (unsigned int i = {0}; i < x_strategy.size(); ++i)
	{
		j = {0};
		for (; it != a_strategy.end(); ++it)
		{
			if (std::find(it, a_strategy.end(), i + 1) != a_strategy.end())
			{
				j++;
			}
			else
				break;
		}
		if (j != 0)
			a_strg_count.insert(std::pair<unsigned int, unsigned int>(j, a_strategy[i]));
	}
	cout << endl << "Рекомендованные стратегии:" << endl;
	for (auto it = a_strg_count.begin(); it != a_strg_count.end(); ++it)
	{
		cout << " a" << it->second << " - " << it->first << " ";
	}
	cout << endl << endl << "Окончательноо, согласно принципу большинства, следует рекомендовать"
	" выбор стратегии a" << (--a_strg_count.end())->second << " — лучшей по " 
	<< (--a_strg_count.end())->first << " из " << a_strategy.size() << " рассмоhренных критериев. ";
}

bool not_reference(double ** arr, unsigned int m)
{
	unsigned int accept = {0};
	for (unsigned int i = {0}; i < m; ++i)
	{
		if (arr[i][0] < 0) // Проверка на допустимое решение
			accept++;
	}
	if (accept > 0)
		return true;
	else
		return false;
}

bool not_optimal(double ** arr, unsigned int m, unsigned int n)
{
	unsigned int check = {0};
	for (unsigned int j = {1}; j < n + 1; ++j)
	{
		if (arr[m][j] > 0) // Проверка на оптимальное решение
			check++;
	}
	if (check > 0)
		return true;
	else 
		return false;
}

std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n)
{
	using std::cout;
	using std::endl;
	using std::vector;
	using std::string;

	cout << endl << "Симплекс-таблица:" << endl << std::setw(9) << "Si0";
	static vector<string> x(m+n+1); // Вектор Xi и F
	if (count == 0)
	{
		for (unsigned int i = {0}; i < m+n; ++i)
		{
			char temp[24]; 
			itoa(i+1,temp,10); // Переводим номер Икса в строку
			x[i] = (chk == 0) ? "U" : "V";
			x[i] += temp;
		}
		x[m+n] = (chk == 0) ? "W " : "Z ";
	}
	if (count > 0)
	{
		for (unsigned int i = {0}; i < n; ++i)
		{
			if (i == rj - 1)
			{
				for (unsigned int j = {0}; j < m+n; ++j)
				{
					if (j == ri)
					{
						string temp = x[i];
						x[i] = x[j+n];
						x[j+n] = temp; // Меняем X-ы, стоящие в разрешающей строке и столбце
					}
				}
			}
		}
	}
	count++;
	
	unsigned int k = {0};
	for (; k < n; ++k)
	{
		cout << std::setw(7) << x[k];
	}
	cout << endl;
	for (unsigned int i = {0}; i < m+1; ++i)
	{
		cout << x[k++];
		for (unsigned int j = {0}; j < n+1; ++j)
		{
			cout << std::setw(7) << std::fixed << std::setprecision(2) << arr[i][j];
		}
		cout << endl;
	}
	return x;
}

void set_new_data(std::vector<double> &c, std::vector<double> &b, std::vector<std::vector<double>>& a)
{
	using std::cout;
	using std::endl;
	using std::vector;

	chk++;
	count = 0;
	std::ofstream fout;
	fout.open("lab6.txt");
	if (!fout.is_open())
	{
		cout << endl << "Ошибка открытия файла.";
		exit(EXIT_FAILURE);
	}

	if (min_max == "max")
	{
		for (unsigned int i = {0}; i < c.size(); ++i)
			c[i] = -c[i];
		fout << "min" << endl;
	}
	else if (min_max == "min")
	{
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			for (unsigned int j = {0}; j < c.size(); ++j)
			{
				a[i][j] = -a[i][j];
			}
		}
		fout << "max" << endl;
	}
	// Транспонируем матрицу системы ограничений A
	vector<vector<double>> temp; // Задаем временный двумерный вектор, в который скопируем наш вектор a
	temp.resize(b.size());
	for (unsigned int i = {0}; i < b.size(); ++i)
		temp[i].resize(c.size());

	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			temp[i][j] = a[i][j]; // Само копирование a в temp
		}
	}

	a.clear(); // Очищаем вектор
	a.resize(c.size());
	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		a[i].resize(b.size()); // Разворачиваем матрицу, т.е. меняем ее размеры
	}

	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		for (unsigned int j = {0}; j < b.size(); ++j)
		{
			a[i][j] = temp[j][i]; // Заполняем вектор a таким образом, что его строка = столбцу вектора temp
		}
	}

	c.swap(b); // Меняем с и b, далее меняем c и b в файле lab6.txt
	// Перезаписываем значения для векторов c,b,a в файле - альтернативный метод транспонирования
	for (unsigned int i = {0}; i < c.size(); ++i)
		fout << c[i] << endl;
	fout << "q" << endl;
	for (unsigned int i = {0}; i < b.size(); ++i)
		fout << b[i] << endl;
	fout << "q" << endl;
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		for (unsigned int j = {0}; j < c.size(); ++j)
			fout << a[i][j] << endl;
	}
	fout.close();
	cout << endl << endl;
	show_new_task(c,b,a);
}

void show_new_task(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>>& a)
{
	using std::cout;
	using std::cin;
	using std::endl;

	cout << endl << "Постановка двойственной задачи (ДЗ) ЛП. Требуется найти решение следующей задачи:" 
	<< endl << std::setw(30) << "Z = b^T * v -> " <<  ((min_max == "max") ? "min" : "max") << endl 
	<< std::setw(30) << "C^T * v <= c'^T" << endl << std::setw(25) << "v >= 0" << endl
	<< "Каноническая форма ДЗ ЛП:" << endl;
	if (min_max == "max")
		cout << std::setw(17) << "Z = ";
	else
		cout << std::setw(17) << "-Z = ";
	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		if (c[i] < 0)
			cout << c[i] << "v" << i+1;
		else
			cout << "+" << c[i] << "v" << i+1;
	}
	cout << " -> ";
	if (min_max == "min")
		cout <<  min_max << endl << endl; 
	else
		cout <<  "min" << endl << endl; 
	unsigned int k = c.size() + 1;
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		cout << "                 ";
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			if (a[i][j] != 0)
			{
				if (a[i][j] != 1 && a[i][j] != -1)
				{
					if (a[i][j] < 0 || j == 0)
						cout << a[i][j] << "v" << j+1;
					else
						if (a[i][j-1] != 0)
							cout << "+" << a[i][j] << "v" << j+1;
						else 
							cout << a[i][j] << "v" << j+1;
				}
				else
				{
					if (a[i][j] < 0)
						cout << "-v" << j+1;
					else if (j != 0)
						cout << "+" << "v" << j+1;
					else 
						cout << "v" << j+1;
				}
		
			}		
		}
		cout << "+v" << k++ << " = " << b[i] << endl;
	}
	cout << endl << "             ";
	for (unsigned int i = {0}; i < c.size()+b.size(); ++i)
	{
		if (i != c.size()+b.size()-1)
			cout << "v" << i+1 << ", ";
		else
			cout << "v" << i+1;
	}
	cout << " >= 0" << endl << "Алгоритм преобразования симплекс-таблицы (жордановы исключения):";

}