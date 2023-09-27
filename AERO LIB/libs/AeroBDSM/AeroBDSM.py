import ctypes
from inspect import Attribute
import os
import struct 

_filepath = os.path.dirname(__file__)

try:
    os.add_dll_directory(_filepath)
except AttributeError:
    pass

# Определение разрядности интерпретатора
bitness = struct.calcsize("P") * 8
#  Подключение DLL
if bitness == 64:
    libstructpy = ctypes.CDLL(_filepath + '/AeroBDSM_x64.dll')
else:
    libstructpy = ctypes.CDLL(_filepath + '/AeroBDSM.dll')

# Таблица сообщений об ошибках
ErrorMessageTable = {
    1: 'число должно быть действительным',
    2: 'число не должно быть отрицательным',
    3: 'число строк/столбцов должно быть больше либо равно 1',
    4: 'ошибка выделения памяти',
    5: 'вычисление прервано пользователем',
    6: 'этот аргумент может быть равен только 0 или 1',
    7: 'угол должен быть меньше 90 градусов',
    8: 'число не должно быть больше 1',
    9: 'этот аргумент не должен быть больше следующего аргумента',
    10: 'число не должно быть меньше 1',
    11: 'этот аргумент может быть равен только 0, 1 или 2',
    12: 'количество элементов в массиве не совпадает с количеством элементов в массиве M',
    13: 'несущие поверхности не должны пересекаться',
    14: 'несущая поверхность не должна выходить за цилиндрическую часть ЛА',
    15: 'несущая поверхность не должна заходить на носовую часть ЛА',
    16: 'Модуль числа должен быть больше либо равен 1',
    17: 'Значение аргумента недопустимо большое',
    18: 'При данных значениях аргументов использовать эту функцию нельзя',
    19: 'Этот аргумент не должен быть равен 0'
}

# класс, который соответствует структуре InputComplex (AeroBDSM.h)
class InputComplex(ctypes.Structure):
    _fields_ = (
        ('Min', ctypes.c_float),
        ('Value', ctypes.c_float),
        ('Max', ctypes.c_float),
    )
    
    def __str__(self) -> str:
        return f'Min = {self.Min}; Value = {self.Value}; Max = {self.Max}'

    def __repr__(self) -> str:
        return self.__str__()

INPUT_COMPLEXES_MAXCOUNT = 9

# класс, который соответствует классу DataResult (AeroBDSM.h)
class DataResult(ctypes.Structure):
    _fields_ = (
        ('ErrorCode', ctypes.c_int),
        ('InvalidArgNumber', ctypes.c_int),
        ('Value', ctypes.c_float),
        ('ExtraValue', ctypes.c_float),
        ('InputComplexesCount', ctypes.c_int),
        ('InputComplexes', InputComplex * INPUT_COMPLEXES_MAXCOUNT),
    )

    def __float__(self):
        return self.Value
    
    def __str__(self) -> str:
        return f"Value = {self.Value}" + f"\nExtraValue = {self.ExtraValue}" + f"\nErrorCode = {self.ErrorCode}"  + f"\nInvalidArgNumber = {self.InvalidArgNumber}\n" + '\n'.join([f'x{i}: {self.InputComplexes[i]}' for i in range(self.InputComplexesCount)])

    def __repr__(self) -> str:
        return self.__str__()

    def __float__(self):        
        return self.Value

    def __lt__(self, other):        
        return self.Value < other

    def __le__(self, other):        
        return self.Value <= other

    def __eq__(self, other):        
        return self.Value == other

    def __ne__(self, other):        
        return self.Value != other

    def __gt__(self, other):        
        return self.Value > other

    def __ge__(self, other):        
        return self.Value >= other

    def __add__(self, other):        
        return self.Value + other

    def __sub__(self, other):        
        return self.Value - other

    def __mul__(self, other):        
        return self.Value * other

    def __truediv__(self, other):        
        return self.Value / other

    def __floordiv__(self, other):        
        return self.Value // other

    def __mod__(self, other):        
        return self.Value % other

    def __divmod__(self, other):        
        return divmod(self.Value,other)

    def __pow__(self, other):        
        return self.Value ** other

    def __radd__(self, other):        
        return other + self.Value 

    def __rsub__(self, other):        
        return other - self.Value

    def __rmul__(self, other):        
        return other * self.Value

    def __rtruediv__(self, other):        
        return other / self.Value

    def __rfloordiv__(self, other):        
        return other // self.Value

    def __rmod__(self, other):        
        return other % self.Value

    def __rdivmod__(self, other):        
        return divmod(other,self.Value)

    def __rpow__(self, other):        
        return other ** self.Value

 # класс, который соответствует классу GeomResult (AeroBDSM.h)
class GeomResult(ctypes.Structure):
    _fields_ = (
        ('ErrorCode', ctypes.c_int),
        ('InvalidArgNumber', ctypes.c_int),
        ('Value', ctypes.c_float),
        ('ExtraValue', ctypes.c_float),        
    )
    
    def __str__(self) -> str:
        return f"Value = {self.Value}" + f"\nExtraValue = {self.ExtraValue}" + f"\nErrorCode = {self.ErrorCode}"  + f"\nInvalidArgNumber = {self.InvalidArgNumber}\n"

    def __repr__(self) -> str:
        return self.__str__()

    def __repr__(self) -> str:
        return self.__str__()

    def __float__(self):        
        return self.Value

    def __lt__(self, other):        
        return self.Value < other

    def __le__(self, other):        
        return self.Value <= other

    def __eq__(self, other):        
        return self.Value == other

    def __ne__(self, other):        
        return self.Value != other

    def __gt__(self, other):        
        return self.Value > other

    def __ge__(self, other):        
        return self.Value >= other

    def __add__(self, other):        
        return self.Value + other

    def __sub__(self, other):        
        return self.Value - other

    def __mul__(self, other):        
        return self.Value * other

    def __truediv__(self, other):        
        return self.Value / other

    def __floordiv__(self, other):        
        return self.Value // other

    def __mod__(self, other):        
        return self.Value % other

    def __divmod__(self, other):        
        return divmod(self.Value,other)

    def __pow__(self, other):        
        return self.Value ** other

    def __radd__(self, other):        
        return other + self.Value 

    def __rsub__(self, other):        
        return other - self.Value

    def __rmul__(self, other):        
        return other * self.Value

    def __rtruediv__(self, other):        
        return other / self.Value

    def __rfloordiv__(self, other):        
        return other // self.Value

    def __rmod__(self, other):        
        return other % self.Value

    def __rdivmod__(self, other):        
        return divmod(other,self.Value)

    def __rpow__(self, other):        
        return other ** self.Value

class Vector(object):
    libstructpy.new_vector.restype = ctypes.c_void_p
    libstructpy.new_vector.argtypes = []
    libstructpy.delete_vector.restype = None
    libstructpy.delete_vector.argtypes = [ctypes.c_void_p]
    libstructpy.vector_size.restype = ctypes.c_int
    libstructpy.vector_size.argtypes = [ctypes.c_void_p]
    libstructpy.vector_get.restype = ctypes.c_float
    libstructpy.vector_get.argtypes = [ctypes.c_void_p, ctypes.c_int]
    libstructpy.vector_push_back.restype = None
    libstructpy.vector_push_back.argtypes = [ctypes.c_void_p, ctypes.c_float]

    def __init__(self, arr=None):
        # инициализация вектора
        self.vector = libstructpy.new_vector()

        # добавление элементов из arr в вектор
        if arr is not None:
            # проверяем, является ли arr iterable-объектом
            try:
                iter(arr)
            except:     # иначе, когда arr является обычным числом
                libstructpy.vector_push_back(self.vector, ctypes.c_float(arr))
            else:       # если arr - iterable-объект
                for val in arr:
                    libstructpy.vector_push_back(
                        self.vector, ctypes.c_float(val))

    # деструктор
    def __del__(self):
        libstructpy.delete_vector(self.vector)

    # длина элементов вектора
    def __len__(self):
        return libstructpy.vector_size(self.vector)

    # возвращает элемент вектора по индексу
    def __getitem__(self, ind):
        if 0 <= ind < len(self):
            return libstructpy.vector_get(self.vector, ctypes.c_int(ind))
        raise IndexError('Vector index out of range')

    # возвращает в виде строки список из элементов
    def __repr__(self):
        return '[{}]'.format(', '.join(str(self[i]) for i in range(len(self))))

    # добавляет элемент в конец вектора
    def push(self, val):
        return libstructpy.vector_push_back(self.vector, ctypes.c_float(val))

def get_Cy_alpha_ConCil(M, lambda_nos, lambda_cil):
    """
Функция для определения производной по углу атаки коэффициента подъёмной силы комбинации конус-цилиндр

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -
lambda_cil : float
    удлинение цилиндрической части фюзеляжа, -

Возврат
-------
Value : float
    производная по углу атаки коэффициента подъёмной силы изолированного фюзеляжа, 1/рад
x[0] : InputComplex
    входной комплекс sqrt(abs(sqr(M) - 1)) / lambda_nos * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_cil / lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 153, Рис.3.2.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cy_alpha_ConCil.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,        
    ]
    libstructpy.get_Cy_alpha_ConCil.restype = DataResult
    
    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cy_alpha_ConCil(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
        ctypes.c_float(lambda_cil),        
    ) 
   
    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cy_alpha_OgiCil(M, lambda_nos, lambda_cil):
    """
Функция для определения производной по углу атаки коэффициента подъёмной силы комбинации оживало-цилиндр

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -
lambda_cil : float
    удлинение цилиндрической части фюзеляжа, -

Возврат
-------
Value : float
    производная по углу атаки коэффициента подъёмной силы изолированного фюзеляжа, 1/рад
x[0] : InputComplex
    входной комплекс sqrt(abs(sqr(M) - 1)) / lambda_nos * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_cil / lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 154, Рис.3.3.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cy_alpha_OgiCil.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cy_alpha_OgiCil.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cy_alpha_OgiCil(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
        ctypes.c_float(lambda_cil),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cy_alpha_SphCil(M, lambda_nos, lambda_cil):
    """
Функция для определения производной по углу атаки коэффициента подъёмной силы цилиндра со сферическим затуплением или с плоским торцом

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа (0 - плоская, 0.5 - сферическая), -
lambda_cil : float
    удлинение цилиндрической части фюзеляжа, -

Возврат
-------
Value : float
    производная по углу атаки коэффициента подъёмной силы изолированного фюзеляжа, 1/рад
x[0] : InputComplex
    входной комплекс sqrt(abs(sqr(M) - 1)) / lambda_nos * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 154, Рис.3.4.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cy_alpha_SphCil.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cy_alpha_SphCil.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cy_alpha_SphCil(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
        ctypes.c_float(lambda_cil),  
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cy_alpha_IsP(M, lambda_c, cc, chi_05, eta_c):
    """
Функция для определения производной по углу атаки коэффициента подъёмной силы изолированной несущей поверхности

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
cc : float
    относительная толщина профиля, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
eta_c : float
    сужение консолей, -

Возврат
-------
Value : float
    производная по углу атаки коэффициента подъёмной силы изолированной несущей поверхности, 1/рад
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * pow(cc, 1/3)
x[2] : InputComplex
    входной комплекс lambda_c * tan(chi_05)

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 156, Рис.3.5.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cy_alpha_IsP.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cy_alpha_IsP.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cy_alpha_IsP(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(cc),
        ctypes.c_float(chi_05),
        ctypes.c_float(eta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_kappa_M(M):
    """
Функция для определения поправочного множителя, учитывающего влияние числа Маха при расчёте коэффициентов интерференции фюзеляжа и несущей поверхности

Параметры
---------
M : float
    число Маха, -

Возврат
-------
Value : float
    поправочный множитель, учитывающий влияние числа Маха при расчёте коэффициентов интерференции фюзеляжа и несущей поверхности, -
x[0] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 162, Рис.3.13.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_kappa_M.argtypes = [        
        ctypes.c_float,
    ]
    libstructpy.get_kappa_M.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_kappa_M(        
        ctypes.c_float(M),        
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Coordinate_zz_v(M, lambda_c, chi_05, zeta_c):
    """
Функция для определения относительной поперечной координаты вихря, сбегающего с консоли несущей поверхности

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
zeta_c : float
    обратное сужение консолей, -

Возврат
-------
Value : float
    относительная поперечная координата вихря, сбегающего с консоли несущей поверхности, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * tan(chi_05)
x[2] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 168, Рис.3.16.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Coordinate_zz_v.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Coordinate_zz_v.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Coordinate_zz_v(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_interference_v(zeta_c, D, l, y_v, z_v):
    """
Функция для определения коэффициента интерференции вихрей и задней несущей поверхности

Параметры
---------
zeta_c : float
    обратное сужение консолей, -
D : float
    диаметр фюзеляжа, м
l : float
    размах несущей поверхности, м
y_v : float
    нормальная координата вихря, м
z_v : float
    поперечная координата вихря, м

Возврат
-------
Value : float
    коэффициент интерференции вихрей и задней несущей поверхности, -
x[0] : InputComplex
    входной комплекс 2 * z_v / l_II
x[1] : InputComplex
    входной комплекс 2 * y_v / l_II
x[2] : InputComplex
    входной комплекс D_II / l_II
x[3] : InputComplex
    входной комплекс zeta_c_II

Ссылки
------
NACA Report 1307 (Chart 7)\\
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 169, Рис.3.17.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_interference_v.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_interference_v.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_interference_v(        
        ctypes.c_float(zeta_c),
        ctypes.c_float(D),
        ctypes.c_float(l),
        ctypes.c_float(y_v),
        ctypes.c_float(z_v),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_kM_Con(M, lambda_nos):
    """
Функция для определения коэффициента торможения потока, вызванного обтеканием конической носовой части

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -

Возврат
-------
Value : float
    Коэффициент торможения потока, вызванного обтеканием конической носовой части, -
x[0] : InputComplex
    входной комплекс lambda_nos
x[1] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 174, Рис.3.21.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_kM_Con.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_kM_Con.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_kM_Con(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_kM_P(M, x, b_Ac):
    """
Функция для определения коэффициента торможения потока, вызванного обтеканием передней несущей поверхности

Параметры
---------
M : float
    число Маха, -
x : float
    расстояние между несущими поверхностями, м
b_Ac : float
    средняя аэродинамическая хорда консолей передней несущей поверхности, м

Возврат
-------
Value : float
    Коэффициент торможения потока, вызванного обтеканием передней несущей поверхности, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс x / b_Ac

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 175, Рис.3.22.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_kM_P.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_kM_P.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_kM_P(        
        ctypes.c_float(M),
        ctypes.c_float(x),
        ctypes.c_float(b_Ac),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_nn_Eff1(l_r, l_c, zeta_c):
    """
Функция для определения коэффициента, характеризующего относительную эффективность концевых рулей

Параметры
---------
l_r : float
    размах рулей, м
l_c : float
    размах консолей несущей поверхности, м
zeta_c : float
    обратное сужение консолей, -

Возврат
-------
Value : float
    коэффициент, характеризующий относительную эффективность концевых рулей, -
x[0] : InputComplex
    входной комплекс l_r / l_c
x[1] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 180, Рис.3.25.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_nn_Eff1.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_nn_Eff1.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_nn_Eff1(        
        ctypes.c_float(l_r),
        ctypes.c_float(l_c),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_nn_Eff2(M, lambda_c, bb_r):
    """
Функция для определения коэффициента, характеризующего относительную эффективность рулей, расположенных вдоль задней кромки

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей несущей поверхности, -
bb_r : float
    относительная хорда руля, -

Возврат
-------
Value : float
    коэффициент, характеризующий относительную эффективность рулей, расположенных вдоль задней кромки, -
x[0] : InputComplex
    входной комплекс bb_r
x[1] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 181, Рис.3.28.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_nn_Eff2.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_nn_Eff2.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_nn_Eff2(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(bb_r),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cp1_Cp2(M, key):
    """
Функция для определения коэффициентов, определяющих давление на поверхности профиля по теории 2-ого приближения

Параметры
---------
M : float
    число Маха, -
key : float
    ключ, определяющий коэффициент, который возвращает функция (0: Cp2/Cp1, 1: Cp1, 2: Cp2)

Возврат
-------
Value : float
    Коэффициенты, определяющие давление на поверхности профиля по теории 2-ого приближения
x[0] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 183, Рис.3.29.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cp1_Cp2.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cp1_Cp2.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cp1_Cp2(        
        ctypes.c_float(M),
        ctypes.c_float(key),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cx_Cil_N(M, alpha):
    """
Функция для определения среднего по длине цилиндра коэффициента сопротивления при обтекании цилиндра по нормали к его оси

Параметры
---------
M : float
    число Маха, -
alpha : float
    угол атаки, рад

Возврат
-------
Value : float
    Средний по длине цилиндра коэффициент сопротивления при обтекании цилиндра по нормали к его оси, -
x[0] : InputComplex
    входной комплекс M * sin(alpha)

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 187, Рис.3.32.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cx_Cil_N.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cx_Cil_N.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cx_Cil_N(        
        ctypes.c_float(M),
        ctypes.c_float(alpha),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_A_IsP(M, zeta_c, Cy_alpha_IsP):
    """
Функция для определения коэффициента дополнительной нормальной силы несущей поверхности при больших углах атаки

Параметры
---------
M : float
    число Маха, -
zeta_c : float
    обратное сужение консолей, -
Cy_alpha_IsP : float
    производная по углу атаки коэффициента подъёмной силы изолированной несущей поверхности, 1/рад

Возврат
-------
Value : float
    коэффициент дополнительной нормальной силы несущей поверхности при больших углах атаки, -
x[0] : InputComplex
    входной комплекс Cy_alpha_IsP
x[1] : InputComplex
    входной комплекс zeta_c
x[2] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 188, Рис.3.35.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_A_IsP.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_A_IsP.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_A_IsP(        
        ctypes.c_float(M),
        ctypes.c_float(zeta_c),
        ctypes.c_float(Cy_alpha_IsP),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cf_M0(Re, xx_t):
    """
Функция для определения коэффициента трения плоской пластинки при М=0

Параметры
---------
Re : float
    число Рейнольдса, -
xx_t : float
    относительная координата точки перехода ламинарного пограничного слоя в турбулентный, -

Возврат
-------
Value : float
    коэффициент трения плоской пластинки при М=0, -
x[0] : InputComplex
    входной комплекс Re
x[1] : InputComplex
    входной комплекс xx_t

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 205, Рис.4.2.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cf_M0.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cf_M0.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cf_M0(        
        ctypes.c_float(Re),
        ctypes.c_float(xx_t),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_M(M, xx_t):
    """
Функция для определения поправочного множителя, учитывающего влияние числа Маха на коэффициент трения плоской пластинки

Параметры
---------
M : float
    число Маха, -
xx_t : float
    относительная координата точки перехода ламинарного пограничного слоя в турбулентный, -

Возврат
-------
Value : float
    поправочный множитель, учитывающий влияние числа Маха на коэффициент трения плоской пластинки, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс xx_t

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 205, Рис.4.3.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_M.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_M.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_M(        
        ctypes.c_float(M),
        ctypes.c_float(xx_t),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Re_t0_1(M, Re, hh):
    """
Функция для определения критического числа Рейнольдса при шероховатой поверхности тела, отсутствии теплопередачи и нулевом градиенте давления

Параметры
---------
M : float
    число Маха, -
Re : float
    число Рейнольдса, -
hh : float
    относительная высота бугорков шероховатости поверхности, -

Возврат
-------
Value : float
    критическое число Рейнольдса при шероховатой поверхности тела, отсутствии теплопередачи и нулевом градиенте давления, -
x[0] : InputComplex
    входной комплекс Re * hh
x[1] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 208, Рис.4.5.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Re_t0_1.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Re_t0_1.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Re_t0_1(        
        ctypes.c_float(M),
        ctypes.c_float(Re),
        ctypes.c_float(hh),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_TTr_TTsl(M, key):
    """
Функция для определения относительной температуры восстановления и относительной температуры поверхности гладкой стенки, при которой ламинарный пограничный слой является устойчивым

Параметры
---------
M : float
    число Маха, -
key : float
    ключ, определяющий относительную температуру, которую возвращает функция (0: Tsl/Tr, 1: Tr/T, 2: Tsl/T)

Возврат
-------
Value : float
    относительные температуры, -
x[0] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 210, Рис.4.8.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_TTr_TTsl.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_TTr_TTsl.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_TTr_TTsl(        
        ctypes.c_float(M),
        ctypes.c_float(key),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_T(M, TT_s):
    """
Функция для определения поправочного множителя, учитывающего влияние температуры поверхности тела на критическое число Рейнольдса

Параметры
---------
M : float
    число Маха, -
TT_s : float
    относительная температура поверхности тела (Ts/Tr), -
x[1] : float
    входной комплекс (TT_s - 1) / powf(M, 2)

Возврат
-------
Value : float
    поправочный множитель, учитывающий влияние температуры поверхности тела на критическое число Рейнольдса
x[0] : InputComplex
    входной комплекс (TT_s - 1) / powf(M, 2)

Ссылки
------
Краснов Н.В. Аэродинамика ракет, 1968 (с. 334, Рис.VI-1-16)\\
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 210, Рис.4.9.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_T.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_T.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_T(        
        ctypes.c_float(M),
        ctypes.c_float(TT_s),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_Con(M, lambda_nos):
    """
Функция для определения поправочного множителя, учитывающего отличие коэффициента трения конуса от коэффициента трения плоской пластинки

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -

Возврат
-------
Value : float
    поправочный множитель, учитывающий отличие коэффициента трения конуса от коэффициента трения плоской пластинки, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 211, Рис.4.10.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_Con.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_Con.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_Con(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxnos_Con(M, lambda_nos):
    """
Функция для определения коэффициента сопротивления давления конической носовой части

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -

Возврат
-------
Value : float
    коэффициент сопротивления давления конической носовой части, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 213, Рис.4.11.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxnos_Con.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxnos_Con.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxnos_Con(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxnos_Par(M, lambda_nos):
    """
Функция для определения коэффициента сопротивления давления параболической носовой части

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -

Возврат
-------
Value : float
    коэффициент сопротивления давления параболической носовой части, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 213, Рис.4.12.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxnos_Par.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxnos_Par.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxnos_Par(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxnos_Ell(M, lambda_nos):
    """
Функция для определения коэффициента сопротивления давления эллиптической носовой части

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -

Возврат
-------
Value : float
    коэффициент сопротивления давления эллиптической носовой части, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 214, Рис.4.13.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxnos_Ell.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxnos_Ell.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxnos_Ell(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Ms_Con(M, teta_Con):
    """
Функция для определения числа Маха на поверхности конуса

Параметры
---------
M : float
    число Маха, -
teta_Con : float
    полуугол при вершине конуса, рад

Возврат
-------
Value : float
    число Маха на поверхности конуса, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс teta_Con * deg

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 217, Рис.4.15.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Ms_Con.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Ms_Con.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Ms_Con(        
        ctypes.c_float(M),
        ctypes.c_float(teta_Con),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxcor_Con(M, lambda_cor, eta_cor):
    """
Функция для определения коэффициента сопротивления давления конической кормовой части

Параметры
---------
M : float
    число Маха, -
lambda_cor : float
    удлинение кормовой части фюзеляжа, -
eta_cor : float
    сужение кормовой части фюзеляжа, -

Возврат
-------
Value : float
    коэффициент сопротивления давления конической кормовой части, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс lambda_cor
x[2] : InputComplex
    входной комплекс eta_cor

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 226, Рис.4.24а.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxcor_Con.argtypes = [
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxcor_Con.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxcor_Con(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_cor),
        ctypes.c_float(eta_cor),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxcor_Par(M, lambda_cor, eta_cor):
    """
Функция для определения коэффициента сопротивления давления параболической кормовой части

Параметры
---------
M : float
    число Маха, -
lambda_cor : float
    удлинение кормовой части фюзеляжа, -
eta_cor : float
    сужение кормовой части фюзеляжа, -

Возврат
-------
Value : float
    коэффициент сопротивления давления параболической кормовой части, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс lambda_cor
x[2] : InputComplex
    входной комплекс eta_cor

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 226, Рис.4.24б.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxcor_Par.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxcor_Par.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxcor_Par(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_cor),
        ctypes.c_float(eta_cor),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxdon_Cil(M, cc, Re_f06, Re_f, xx_tf, lambda_f):
    """
Функция для определения коэффициента донного сопротивления фюзеляжа с цилиндрической кормовой частью

Параметры
---------
M : float
    число Маха, -
cc : float
    относительная толщина профиля задней несущей поверхности, -
Re_f06 : float
    число Рейнольдса для фюзеляжа при скорости, соответствующей М=0.6, -
Re_f : float
    число Рейнольдса для фюзеляжа при скорости, соответствующей М, -
xx_tf : float
    относительная координата точки перехода ламинарного пограничного слоя в турбулентный на фюзеляже, -
lambda_f : float
    удлинение фюзеляжа, -

Возврат
-------
Value : float
    коэффициент донного сопротивления фюзеляжа с цилиндрической кормовой частью, -
x[0] : InputComplex
    входной комплекс M
x[1] : InputComplex
    входной комплекс cc

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 228, Рис.4.26.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxdon_Cil.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxdon_Cil.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxdon_Cil(        
        ctypes.c_float(M),
        ctypes.c_float(cc),
        ctypes.c_float(Re_f06),
        ctypes.c_float(Re_f),
        ctypes.c_float(xx_tf),
        ctypes.c_float(lambda_f),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_eta(M, lambda_cor, eta_cor):
    """
Функция для определения поправочного множителя, учитывающего влияние сужающейся кормовой части на коэффициент донного сопротивления фюзеляжа

Параметры
---------
M : float
    число Маха, -
lambda_cor : float
    удлинение кормовой части фюзеляжа, -
eta_cor : float
    сужение кормовой части фюзеляжа, -

Возврат
-------
Value : float
    поправочный множитель, учитывающий влияние сужающейся кормовой части на коэффициент донного сопротивления фюзеляжа, -
x[0] : InputComplex
    входной комплекс (1 - eta_cor) / (2 * lambda_cor * sqr(eta_cor))
x[1] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 230, Рис.4.27.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_eta.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_eta.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_eta(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_cor),
        ctypes.c_float(eta_cor),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_c(cc, xx_t):
    """
Функция для определения поправочного множителя, учитывающего влияние толщины профиля на коэффициент профильного сопротивления несущей поверхности

Параметры
---------
cc : float
    относительная толщина профиля, -
xx_t : float
    относительная координата точки перехода ламинарного пограничного слоя в турбулентный, -

Возврат
-------
Value : float
    поправочный множитель, учитывающий влияние толщины профиля на коэффициент профильного сопротивления несущей поверхности, -
x[0] : InputComplex
    входной комплекс cc
x[1] : InputComplex
    входной комплекс xx_t

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 232, Рис.4.28.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_c.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_c.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_c(        
        ctypes.c_float(cc),
        ctypes.c_float(xx_t),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_chi(chi_0):
    """
Функция для определения поправочного множителя, учитывающего влияние угла стреловидности по передней кромке на критическое число Рейнольдса

Параметры
---------
chi_0 : float
    угол стреловидности по передней кромке несущей поверхности, рад

Возврат
-------
Value : float
    поправочный множитель, учитывающий влияние угла стреловидности по передней кромке на критическое число Рейнольдса
x[0] : InputComplex
    входной комплекс chi_0 * deg

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 232, Рис.4.29.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_chi.argtypes = [        
        ctypes.c_float,        
    ]
    libstructpy.get_Sigma_chi.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_chi(        
        ctypes.c_float(chi_0),        
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxw_Rhomb_M1(M, cc, zeta_c, chi_c, lambda_c):
    """
Функция для определения коэффициента волнового сопротивления несущей поверхности с ромбовидным профилем при М >= 1

Параметры
---------
M : float
    число Маха, -
cc : float
    относительная толщина профиля, -
zeta_c : float
    обратное сужение консолей, -
chi_c : float
    угол стреловидности по линии максимальных толщин, рад
lambda_c : float
    удлинение консолей несущей поверхности, -

Возврат
-------
Value : float
    коэффициент волнового сопротивления несущей поверхности с ромбовидным профилем при М >= 1, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(sqr(M) - 1)
x[1] : InputComplex
    входной комплекс lambda_c * pow(cc, 1/3)
x[2] : InputComplex
    входной комплекс lambda_c * tan(chi_c)
x[3] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 234, Рис.4.30.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxw_Rhomb_M1.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Cxw_Rhomb_M1.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxw_Rhomb_M1(        
        ctypes.c_float(M),
        ctypes.c_float(cc),
        ctypes.c_float(zeta_c),
        ctypes.c_float(chi_c),
        ctypes.c_float(lambda_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_pw(M, chi_c):
    """
Функция для определения поправочного множителя, учитывающего степень влияния формы профиля на коэффициент волнового сопротивления несущей поверхности

Параметры
---------
M : float
    число Маха, -
chi_c : float
    угол стреловидности по линии максимальных толщин, рад

Возврат
-------
Value : float
    поправочный множитель, учитывающий степень влияния формы профиля на коэффициент волнового сопротивления несущей поверхности, -
x[0] : InputComplex
    входной комплекс sqrt(sqr(M) - 1) - tan(chi_c)

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 240, Рис.4.32.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_pw.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_pw.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_pw(        
        ctypes.c_float(M),
        ctypes.c_float(chi_c),        
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Mcr_pr(cc, Cn, xx_c):
    """
Функция для определения критического числа Маха для симметричных дозвуковых профилей

Параметры
---------
cc : float
    относительная толщина профиля, -
Cn : float
    коэффициент нормальной силы консолей, -
xx_c : float
    относительная координата линии максимальных толщин, -

Возврат
-------
Value : float
    критическое число Маха для симметричного дозвукового профиля, -
x[0] : InputComplex
    входной комплекс Cn
x[1] : InputComplex
    входной комплекс cc
x[2] : InputComplex
    входной комплекс xx_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 241, Рис.4.34.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Mcr_pr.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Mcr_pr.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Mcr_pr(        
        ctypes.c_float(cc),
        ctypes.c_float(Cn),
        ctypes.c_float(xx_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_DELTA_Mcr0_lambda(Mcr_pr0, lambda_c):
    """
Функция для определения поправки к критическому числу Маха, учитывающей конечность размаха несущей поверхности

Параметры
---------
Mcr_pr0 : float
    критическое число Маха профиля при нулевом коэффициенте нормальной силы консолей, -
lambda_c : float
    удлинение консолей несущей поверхности, -

Возврат
-------
Value : float
    поправка к критическому числу Маха, учитывающая конечность размаха несущей поверхности, -
x[0] : InputComplex
    входной комплекс Mcr_pr0
x[1] : InputComplex
    входной комплекс lambda_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 242, Рис.4.35.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_DELTA_Mcr0_lambda.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_DELTA_Mcr0_lambda.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_DELTA_Mcr0_lambda(        
        ctypes.c_float(Mcr_pr0),
        ctypes.c_float(lambda_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_DELTA_Mcr0_chi(Mcr_pr0, chi_c):
    """
Функция для определения поправки к критическому числу Маха, учитывающей стреловидность несущей поверхности

Параметры
---------
Mcr_pr0 : float
    критическое число Маха профиля при нулевом коэффициенте нормальной силы консолей, -
chi_c : float
    угол стреловидности по линии максимальных толщин, рад

Возврат
-------
Value : float
    поправка к критическому числу Маха, учитывающая стреловидность несущей поверхности, -
x[0] : InputComplex
    входной комплекс Mcr_pr0
x[1] : InputComplex
    входной комплекс chi_c * deg

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 242, Рис.4.36.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_DELTA_Mcr0_chi.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_DELTA_Mcr0_chi.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_DELTA_Mcr0_chi(        
        ctypes.c_float(Mcr_pr0),
        ctypes.c_float(chi_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Cxdon_pr(M):
    """
Функция для определения коэффициента донного сопротивления, создаваемого затупленной задней кромкой профиля несущей поверхности

Параметры
---------
M : float
    число Маха, -

Возврат
-------
Value : float
    коэффициент донного сопротивления, создаваемого затупленной задней кромкой профиля несущей поверхности, -
x[0] : InputComplex
    входной комплекс M

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 243, Рис.4.38.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Cxdon_pr.argtypes = [        
        ctypes.c_float,
    ]
    libstructpy.get_Cxdon_pr.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Cxdon_pr(        
        ctypes.c_float(M),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_Cxinos(M, lambda_nos, key):
    """
Функция для определения множителя, учитывающего форму носовой части и число Маха при расчёте коэффициента тангенциальной силы, индуцированной перераспределением давления на носовой части при ненулевых углах атаки

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -
nos_typ : float
    тип носовой части (0: коническая, 1: оживальная), -

Возврат
-------
Value : float
    множитель, учитывающий форму носовой части и число Маха, -
x[0] : InputComplex
    входной комплекс sqrt(abs(sqr(M) - 1)) / lambda_nos * sign(M - 1)
x[1] : InputComplex
    входной комплекс nos_typ

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 245, Рис.4.40.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_Cxinos.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_Cxinos.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_Cxinos(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
        ctypes.c_float(key),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_CC_F_IsP(M, lambda_c, chi_0):
    """
Функция для определения коэффициента пропорциональности подсасывающей силы консолей несущей поверхности

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей несущей поверхности, -
chi_0 : float
    угол стреловидности по передней кромке несущей поверхности, рад

Возврат
-------
Value : float
    коэффициент пропорциональности подсасывающей силы консолей несущей поверхности, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * tan(chi_0)

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 246, Рис.4.42.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_CC_F_IsP.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_CC_F_IsP.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_CC_F_IsP(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_0),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Sigma_CF(M, chi_0, alpha):
    """
Функция для определения множителя, учитывающего неполноту реализации подсасывающей силы консолей несущей поверхности

Параметры
---------
M : float
    число Маха, -
chi_0 : float
    угол стреловидности по передней кромке несущей поверхности, рад
alpha : float
    угол атаки несущей поверхности, рад

Возврат
-------
Value : float
    множитель, учитывающий неполноту реализации подсасывающей силы консолей несущей поверхности, -
x[0] : InputComplex
    входной комплекс sqrt(abs(sqr(M) - 1)) / tan(chi_0) * sign(M - 1)
x[1] : InputComplex
    входной комплекс alpha * deg

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 246, Рис.4.43.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Sigma_CF.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Sigma_CF.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Sigma_CF(        
        ctypes.c_float(M),
        ctypes.c_float(chi_0),
        ctypes.c_float(alpha),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_xi_M(M, lambda_nos, lambda_cil):
    """
Функция для определения множителя, учитывающего влияние числа Маха на смещение фокуса комбинации носовой части и цилиндра

Параметры
---------
M : float
    число Маха, -
lambda_nos : float
    удлинение носовой части фюзеляжа, -
lambda_cil : float
    удлинение цилиндрической части фюзеляжа, -

Возврат
-------
Value : float
    множитель, учитывающий влияние числа Маха на смещение фокуса комбинации носовой части и цилиндра, -
x[0] : InputComplex
    входной комплекс sqrt(abs(sqr(M) - 1)) / lambda_nos * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_cil / lambda_nos

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 264, Рис.5.7.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_xi_M.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_xi_M.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_xi_M(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_nos),
        ctypes.c_float(lambda_cil),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Coordinate_xxF_IsP(M, lambda_c, chi_05, zeta_c):
    """
Функция для определения относительной координаты фокуса изолированной несущей поверхности

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
zeta_c : float
    обратное сужение консолей, -

Возврат
-------
Value : float
    относительная координата фокуса изолированной несущей поверхности, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс zeta_c
x[2] : InputComplex
    входной комплекс lambda_c * tan(chi_05)

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 265, Рис.5.8.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Coordinate_xxF_IsP.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Coordinate_xxF_IsP.restype = DataResult
    
    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Coordinate_xxF_IsP(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_xi_D(DD):
    """
Функция для определения множителя, учитывающего влияние относительного диаметра фюзеляжа на координату приложения дополнительной нормальной силы консолей

Параметры
---------
DD : float
    относительный диаметр фюзеляжа в районе несущей поверхности, -

Возврат
-------
Value : float
    множитель, учитывающий влияние относительного диаметра фюзеляжа на координату приложения дополнительной нормальной силы консолей, -
x[0] : InputComplex
    входной комплекс DD

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 268, Рис.5.11.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_xi_D.argtypes = [        
        ctypes.c_float,
    ]
    libstructpy.get_xi_D.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_xi_D(        
        ctypes.c_float(DD),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_mm_z_demph_isP(M, lambda_c, chi_05, zeta_c):
    """
 Функция для определения демпфирующего момента изолированной несущей поверхности относительно оси, проходящей через середину САХ
Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
zeta_c : float
    обратное сужение консолей, -
Возврат
-------
Value : float
    производная коэффициента момента по угловой скорости, отнесённая к производной коэффициента нормальной силы по углу атаки, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * tan(chi_05)
x[2] : InputComplex
    входной комплекс zeta_c
Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 278, Рис.5.15.)
"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_mm_z_demph_isP.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_mm_z_demph_isP.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_mm_z_demph_isP(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Bx_demph_IsP(M, lambda_c, chi_05):
    """
Функция для определения коэффициента, учитывающего изменение демпфирующего момента изолированной несущей поверхности с заострёнными концами относительно оси, смещённой от середины САХ
Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей несущей поверхности, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
Возврат
-------
Value : float
    коэффициент, учитывающий изменение демпфирующего момента изолированной несущей поверхности с заострёнными концами относительно оси, смещённой от середины САХ, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * tan(chi_05)
Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 278, Рис.5.16.)
"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Bx_demph_IsP.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Bx_demph_IsP.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Bx_demph_IsP(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_dCyalf_IsP_B1(M, lambda_c, chi_05):
    """
Функция для определения коэффициента B1, служащего для расчёта прирашения производной подъёмной силы ЛА при его обтекании со скольжением
Параметры
---------

M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад

Возврат
-------
Value : float
     коэффициент B1, -

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 298, Рис.6.2.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_dCyalf_IsP_B1.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_dCyalf_IsP_B1.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_dCyalf_IsP_B1(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result

        
def get_dCyalf_IsP_B2(M, lambda_c, chi_05):
    """
Функция для определения коэффициента B2, служащего для расчёта прирашения производной подъёмной силы ЛА при его обтекании со скольжением
Параметры
---------

M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад

Возврат
-------
Value : float
     коэффициент B2, -

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 298, Рис.6.3.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_dCyalf_IsP_B2.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_dCyalf_IsP_B2.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_dCyalf_IsP_B2(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Coordinate_zz_d(M, lambda_c, chi_05, zeta_c):
    """
Функция для определения относительной поперечной координаты центра давления консоли несущей поверхности

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
zeta_c : float
    обратное сужение консолей, -

Возврат
-------
Value : float
    относительная поперечная координата центра давления консоли несущей поверхности, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * tan(chi_05)
x[2] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 300, Рис.6.4.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Coordinate_zz_d.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Coordinate_zz_d.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Coordinate_zz_d(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_nu_DD(M, lambda_c, DD, zeta_c):
    """
Функция для определения коэффициента взаимного влияния правой и левой консолей,	в зависимости от относительного диаметра фюзеляжа
Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
DD : float
    относительный диаметр фюзеляжа в районе несущей поверхности, -
zeta_c : float
    обратное сужение консолей, -
Возврат
-------
Value : float
    коэффициент взаимного влияния правой и левой консолей, по мере уменьшения относительного диаметра фюзеляжа, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс DD
x[2] : InputComplex
    входной комплекс zeta_c
Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 303, Рис.6.7.)
"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_nu_DD.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_nu_DD.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_nu_DD(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(DD),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Coordinate_yy_v(xx, alpha, epsilon_m, yy_0):
    """
Функция для определения относительных координат вихрей сбегающих с поверхностей верхних и нижних консолей
Параметры
---------

xx : float
    относительная координата вихря по Ox, -
alpha : float
    угол атаки, рад
epsilon_m : float
    осреднённый угол скоса потока от горизонтальных консолей, рад
yy_0 : float
    относительная координата схода вихря по Oy, -

Возврат
-------
Value : float
     относительная координата вихря по Oy в районе задних консолей, -

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 315, Рис.6.14.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Coordinate_yy_v.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Coordinate_yy_v.restype = GeomResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Coordinate_yy_v(        
        ctypes.c_float(xx),
        ctypes.c_float(alpha),
        ctypes.c_float(epsilon_m),
        ctypes.c_float(yy_0),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_kappa_DD(DD):
    """
Функция для определения коэффициента, учитывающего взаимное влияние горизонтальных и вертикальных консолей 
Параметры
---------
DD : float
    относительный диаметр фюзеляжа в районе несущей поверхности, -
Возврат
-------
Value : float
   коэффициент, учитывающий взаимное влияние горизонтальных и вертикальных консолей, -
x[0] : InputComplex
    входной комплекс DD
Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 318, Рис.6.17.)
"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_kappa_DD.argtypes = [        
        ctypes.c_float,
    ]
    libstructpy.get_kappa_DD.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_kappa_DD(        
        ctypes.c_float(DD),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Integral_m_x(y_v, l, D, zeta_c):
    """
Функция для определения интеграла подъёмной силы по размаху, создающей момент крена и вызванной вихрем, сходящим с передней несущей поверхности

Параметры
---------
y_v : float
    координата вихря в районе задних консолей, м
l : float
    полный размах задней несущей поверхности, м
D : float
    диаметр фюзеляжа в районе задних консолей, м
zeta_c : float
    обратное сужение задних консолей, -

Возврат
-------
Value : float
    относительное значение интеграла, -
x[0] : InputComplex
    входной комплекс (y_v - D/2) / (l/2 - D/2)
x[1] : InputComplex
    входной комплекс D/l
x[2] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 320, Рис.6.18.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_Integral_m_x.argtypes = [
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_Integral_m_x.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_Integral_m_x(
        ctypes.c_float(y_v),
        ctypes.c_float(l),
        ctypes.c_float(D),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_nn_Eff4(zz_e, zeta_c):
    """
Функция для определения коэффициента, учитывающего влияние относительного размаха элеронов на их эффективность
Параметры
---------

zz_e : float
    относительная координата конца(наружного или внутреннего) элерона, -
zeta_c : float
    обратное сужение консолей, -

Возврат
-------
Value : float
     коэффициент, учитывающий влияние относительного размаха элеронов на их эффективность, -
x[0] : InputComplex
    входной комплекс zz_e
x[1] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 326, Рис.6.20.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_nn_Eff4.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_nn_Eff4.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_nn_Eff4(        
        ctypes.c_float(zz_e),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_m_x_delta_el(M, lambda_c, lambda_e, eta_c, chi_05, cc, bb_r, zz_e_nar, zz_e_vn, key):
    """
 Функция для определения производной момента крена по углу отклонения элеронов
Параметры
---------

M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
lambda_e : float
    удлинение элеронов, -
eta_c : float
    сужение консолей
chi_05 : float
    угол стреловидности по линии середин хорд, рад
cc : float
    относительная толщина профиля, -
bb_r : float
    относительная хорда руля, -
zz_e_nar : float
    относительная координата конца(наружного) элерона, -
zz_e_vn : float
    относительная координата конца(внутреннего) элерона, -
key : float
    ключ, определяющий, что возвращает функция (0: производная момента, 1: L и N)

Возврат
-------
Value : float
     производная момента крена элеронов, 1/рад


Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 327, Рис.6.21.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_m_x_delta_el.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_m_x_delta_el.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_m_x_delta_el(        
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(lambda_e),
        ctypes.c_float(eta_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(cc),
        ctypes.c_float(bb_r),
        ctypes.c_float(zz_e_nar),
        ctypes.c_float(zz_e_vn),
        ctypes.c_float(key),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_mm_x_demph_IsP(M, lambda_c, chi_05, zeta_c):
    """
Функция для определения демпфирующего момента крена изолированной несущей поверхности

Параметры
---------
M : float
    число Маха, -
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
zeta_c : float
    обратное сужение консолей, -

Возврат
-------
Value : float
    производная коэффициента момента по угловой скорости, отнесённая к производной коэффициента нормальной силы по углу атаки, -
x[0] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
x[1] : InputComplex
    входной комплекс lambda_c * tan(chi_05)
x[2] : InputComplex
    входной комплекс zeta_c

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 330, Рис.6.24.)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_mm_x_demph_IsP.argtypes = [
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_mm_x_demph_IsP.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_mm_x_demph_IsP(
        ctypes.c_float(M),
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(zeta_c),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result

    
def get_mm_x_spiral_zeta0(lambda_c, chi_05, M):
    """
Функция для определения спирального момента крена изолированной несущей поверхности с заострёнными концами
Параметры
---------
lambda_c : float
    удлинение консолей, -
chi_05 : float
    угол стреловидности по линии середин хорд, рад
M : float
    число Маха, -
Возврат
-------
Value : float
    производная коэффициента момента по угловой скорости, отнесённая к коэффициенту нормальной силы, 1/рад
x[0] : InputComplex
    входной комплекс lambda_c * tan(chi_05)
x[1] : InputComplex
    входной комплекс lambda_c
x[2] : InputComplex
    входной комплекс lambda_c * sqrt(abs(sqr(M) - 1)) * sign(M - 1)
Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 335, Рис.6.29.)
"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_mm_x_spiral_zeta0.argtypes = [
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
    ]
    libstructpy.get_mm_x_spiral_zeta0.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_mm_x_spiral_zeta0(
        ctypes.c_float(lambda_c),
        ctypes.c_float(chi_05),
        ctypes.c_float(M),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result
   

def get_k_sch(M):
    """
 Функция для определения коэффициента потерь нормальной силы рулей по сравнению с изолированным крылом за счет перетекания воздуха через щель
Параметры
---------

M : float
    число Маха, -

Возврат
-------
Value : float
     коэффициент, учитывающий потери нормальной силы рулей за счет перетекания воздуха через щель, -

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 179)

"""
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_k_sch.argtypes = [        
        ctypes.c_float,
    ]
    libstructpy.get_k_sch.restype = DataResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_k_sch(        
        ctypes.c_float(M),
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result
    

def get_S_bok(D, S_nos, lambda_cil, lambda_cor, eta_cor, M, b_b, L_hv):
    """
Функция для определения незатенённой площади боковой проекции фюзеляжа

Параметры
---------
D : float
    диаметр фюзеляжа, м
S_nos : float
    площадь боковой проекции носовой части, м2
lambda_cil : float
    удлинение цилиндрической части, -
lambda_cor : float
    удлинение кормовой части, -
eta_cor : float
    сужение кормовой части, -
M : float
    массив чисел Маха набегающего потока на каждую несущую поверхность, -
b_b : float
    массив длин бортовой хорды несущих поверхностей, м
L_hv : float
    массив расстояний от конца бортовой хорды несущей поверхности до донного среза фюзеляжа, м

Возврат
-------
Value : float
    площадь незатенённой области фюзеляжа, м2

Ссылки
------
Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 185 - 186)

"""
    M = Vector(M)
    b_b = Vector(b_b)
    L_hv = Vector(L_hv)

    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_S_bok.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    libstructpy.get_S_bok.restype = GeomResult  

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_S_bok(        
        ctypes.c_float(D),
        ctypes.c_float(S_nos),
        ctypes.c_float(lambda_cil),
        ctypes.c_float(lambda_cor),
        ctypes.c_float(eta_cor),
        M.vector,
        b_b.vector,
        L_hv.vector,
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_x_c_pl_F(D, S_nos, x_c_pl_nos, lambda_nos, lambda_cil, lambda_cor, eta_cor, M, b_b, L_hv):
    """
    Функция для определения координаты x центра тяжести незатенённой площади боковой проекции фюзеляжа

    Параметры
    ---------
    D : float
        диаметр фюзеляжа, м
    S_nos : float
        площадь боковой проекции носовой части, м2
    x_c_pl_nos : float
        координата x центра тяжести площади боковой проекции носовой части S_nos, м
    lambda_nos : float
        удлинение носовой части фюзеляжа, -
    lambda_cil : float
        удлинение цилиндрической части, -
    lambda_cor : float
        удлинение кормовой части, -
    eta_cor : float
        сужение кормовой части, -
    M : float
        массив чисел Маха набегающего потока на каждую несущую поверхность, -
    b_b : float
        массив длин бортовой хорды несущих поверхностей, м
    L_hv : float
        массив расстояний от конца бортовой хорды несущей поверхности до донного среза фюзеляжа, м

    Возврат
    -------
    Value : float
        координата x центра тяжести незатенённой площади боковой проекции фюзеляжа, м

    Ссылки
    ------
    Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 271)

    """
    M = Vector(M)
    b_b = Vector(b_b)
    L_hv = Vector(L_hv)
    
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_x_c_pl_F.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    libstructpy.get_x_c_pl_F.restype = GeomResult

    #  Вызов функций, описанных в DLL
    Result = libstructpy.get_x_c_pl_F(        
        ctypes.c_float(D),
        ctypes.c_float(S_nos),
        ctypes.c_float(x_c_pl_nos),
        ctypes.c_float(lambda_nos),
        ctypes.c_float(lambda_cil),
        ctypes.c_float(lambda_cor),
        ctypes.c_float(eta_cor),
        M.vector,
        b_b.vector,
        L_hv.vector
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_psi_eps(M_I,z_v,y_v,x_zI_II,phi,D_II,l_c_II,eta_c_II,b_b_II,chi_0_II):
    """
    Функция для определения множителя скоса потока

    Параметры
    ---------
    M_I : float
        число Маха набегающего потока в области I несущей поверхности, -
    z_v : float
        поперечная координата вихря передней консоли, м
    y_v : float
        вертикальная координата вихря передней консоли, м
    x_zI_II : float
        расстояние вдоль оси x от точки схода вихря до начала бортовой хорды II несущей поверхности, м
    phi : float
        поперечный угол задней консоли относительно передней, рад
    D_II : float
        диаметр фюзеляжа в области II несущей поверхности, м
    l_c_II : float
        длина консолей II несущей поверхности, м
    eta_c_II : float
        сужение консолей II несущей поверхности, -
    b_b_II : float
        длина бортовой хорды II несущей поверхности, м
    chi_0_II : float
        угол стреловидности по передней кромке II несущей поверхности, рад

    Возврат
    -------
    Value : float
        множитель скоса потока, -

    Ссылки
    ------
    Лебедев, Чернобровкин. Динамика полёта БПЛА, 1973 (с. 172)

    """
    #  Указание типов аргументов и возвращаемого значения
    libstructpy.get_psi_eps.argtypes = [        
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_float
    ]
    libstructpy.get_psi_eps.restype = GeomResult

     #  Вызов функций, описанных в DLL
    Result = libstructpy.get_psi_eps(        
        ctypes.c_float(M_I),
        ctypes.c_float(z_v),
        ctypes.c_float(y_v),
        ctypes.c_float(x_zI_II),
        ctypes.c_float(phi),
        ctypes.c_float(D_II),
        ctypes.c_float(l_c_II),
        ctypes.c_float(eta_c_II),
        ctypes.c_float(b_b_II),
        ctypes.c_float(chi_0_II)
    )

    if Result.ErrorCode in ErrorMessageTable:
        raise ValueError(
            f'Аргумент {Result.InvalidArgNumber}, {ErrorMessageTable[Result.ErrorCode]}.')
    else:
        return Result


def get_Info():
    """
Функция для получения информации о библиотеке
"""   
    libstructpy.get_Info.argtypes = None   
    libstructpy.get_Info.restype = ctypes.c_char_p

    #  Вызов функций, описанных в DLL
    return libstructpy.get_Info()
