import numpy as np
import matplotlib.pyplot as plt

#считываем файл НЕОБХОДЖИМО ВСТАВИТЬ СВОЁ НАЗВАНИЕ и читаем его как текст
file = open("dannye.txt")
data = file.read()

#костыль для считывания файла в преподнесенном формате
#здесь я просто считываюд всё как текст, потом его записываю,
#как мссиив чисел, где каждое число по отдельности
inumber=1
listnumber=[]
for i in range(len(data)):
    if ((data[i] == ",") or (data[i] == "]")):
        listnumber.append(float(data[inumber:i]))
        inumber=i+2
        
# массив значений горизонтальной координаты        
xlist = np.linspace(1, len(listnumber), len(listnumber))
#задаем размер холста
fig, ax = plt.subplots(figsize=(15,10))
#строим x=xlist, y=listnumber 
ax.plot (xlist, listnumber, linewidth = 2, label="Вставить своё название 1")
#подпись осей
ax.set_xlabel('Вставить своё название 1', rotation=0, fontsize=25)
ax.set_ylabel('Вставить своё название 2', rotation=90, fontsize=25)
#диапазон осей
ax.set_ylim(bottom=0, top=max(listnumber)+max(listnumber)/10)
ax.set_xlim(left=0, right=len(listnumber))
#сетку grid можно убрать
ax.grid(color='black', linestyle='-', linewidth=0.2)
#отображение графика
plt.show()
#сохранение в файл рядом с данным скриптом
fig.savefig('Raspredelenie.png', dpi=170)

file.close()
