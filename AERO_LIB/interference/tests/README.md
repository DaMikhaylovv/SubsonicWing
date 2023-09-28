#### Импортирование библиотеки

Для импортирования библиотеки вставьте следующий код:
```python
import os
import sys
# если блокнот с тестом на 1 уровень ниже, чем __LA.py, то '..', если на 2 - '../..' и т.д.
AeroBDSM_OOP_dir = os.path.abspath(os.path.join('../..'))
if AeroBDSM_OOP_dir not in sys.path:
    sys.path.append(AeroBDSM_OOP_dir)
```