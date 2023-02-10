import json
import xlwt
from xlwt import Workbook
import os
directory = os.getcwd()
print("The current working directory of the file is : ", directory)
file = open('data_sardegna.json','r')
data = json.load(file)
diz1 = data['outputs']['hourly']
wb = Workbook()
sheet1 = wb.add_sheet('Sheet 1')
for i in range(len(diz1)):
    diz = diz1[i]
    sheet1.write(i,0,i)
    sheet1.write(i,1,diz['G(i)'])
wb.save('data_sardegna.xls')