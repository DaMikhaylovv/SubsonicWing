from libs import AeroBDSM
from libs.handler import *
from math import *
from typing import Dict
from .__Profile import Profile

class SubsonicProfile(Profile):
    '''
    Класс расчета АДХ дозвукового профиля
    '''
# TODO: ситуация с дозвуковыми профилями сложная. Буду думать как задавать параметры
# TODO: поместить сюда Xfoil
    def __init__(self, b: float, cc: float, xx_c: float, rr_0: float, hb: float, profile_name: str) -> None:
        '''
        Создание полей объекта класса и расчет геометрических параметров изолированного крыла бесконечного размаха c профилем, образованным дугами окружности или параболы

        Ввод:   b: float - длина хорды крыла, м
                cc: float - относительная толщина профиля крыла в долях хорды
                xx_c: float - относительная координата положения максимальной толщины профиля крыла в долях хорды
                rr_0: float - относительный радиус скругления передней кромки крыла в долях ... ?
                hb: float - относительная высота бугорков шероховатой поверхности тела в долях хорды
                profile_name: str - наименование профиля (например, NACA0012) или название файла с координатами точек профиля
        Вывод:  ...
        '''
        self.profile_name = profile_name
        super().__init__(b, cc, xx_c, rr_0, hb)

    def geometry(self) -> Dict:
        '''
        Вывод геометрических параметров профиля изолированного крыла

        Ввод:   ...
        Вывод:  
                b: float - длина хорды крыла, м
                cc: float - относительная толщина профиля крыла в долях хорды
                xx_c: float - относительная координата положения максимальной толщины профиля крыла в долях хорды
                rr_0: float - относительный радиус скругления передней кромки крыла в долях ... ?
                hb: float - относительная высота бугорков шероховатой поверхности тела в долях хорды
                profile_name: str - наименование профиля (например, NACA0012) или название файла с координатами точек профиля
        '''
        # Результат
        result = {
            'b': self.b,
            'cc': self.cc,
            'xx_c': self.xx_c,
            'rr_0': self.rr_0,
            'hb': self.hb,
            'profile_name': self.profile_name
            }

        # Результат
        return result    

    # @checkSavedValue
    def K_w(self) -> float:
        '''
        Расчет коэффициента, представляющего собой отношение коэффициента силы волнового сопротивления данного профиля к коэффициенту силы волнового сопротивления ромбовидного профиля

        Ввод:   ...
        Вывод:  K_w: float - коэффициент, представляющий собой отношение коэффициента силы волнового сопротивления данного профиля к коэффициенту силы волнового сопротивления ромбовидного профиля
        '''

        # Коэффициент, представляющий собой отношение коэффициента силы волнового сопротивления дозвукового профиля к коэффициенту силы волнового сопротивления ромбовидного профиля изолированного крыла бесконечного размаха [-]
        K_w = (2.5 + 4) / 2

        # Результат
        return K_w
