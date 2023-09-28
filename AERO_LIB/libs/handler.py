from numpy import sqrt
from .Atmosphere import atmo
from scipy.special import erf

def Mach_to_v(M: float, H: float = 0) -> float:
    '''Конвертирует скорость в Махах в скорость в м/с, учитывая высоту (стандартное значение высоты = 0 м)'''
    return M * atmo.a(H)

def v_to_Mach(v: float, H: float = 0) -> float:
    '''Конвертирует скорость в м/с в скорость в Махах, учитывая высоту (стандартное значение высоты = 0 м)'''
    return v / atmo.a(H)

def Phi(z):
    '''Функция Лапласа-Гаусса (функция распределения)'''
    return 0.5 * ( 1 + erf(z / sqrt(2)) )

def q(M: float, H: float = 0) -> float:
    '''
    Расчет скоростного напора
    [Литература: А. А. Лебедев и Л.С. Чернобровкин - "Динамика полета" (стр. 33)]

    Ввод:   H: float - высота полета [м]
            M: float - число Маха [-]
    Вывод:  q: float - скоростной напор [кг/(м*с^2)]
    '''
    v = Mach_to_v(M, H)
    q = atmo.rho(H) * v * v / 2
    return q


def checkSavedValue(func: callable):
    __savedValuesKeys = __savedValues.keys()
    def wrapper(*args, **kwargs):
        nonlocal __savedValuesKeys
        funcName = func.__module__ + '.' + func.__name__
        # print(funcName)
        if not funcName in __savedValuesKeys:
            res = func(*args, **kwargs)
            __savedValues[funcName] = (args, res)
            __savedValuesKeys = __savedValues.keys()
            return res
        if args == __savedValues[funcName][0]:
            return __savedValues[funcName][1]
        res = func(*args, **kwargs)
        __savedValues[funcName] = (args, res)
        return res
    return wrapper

__savedValues = {}