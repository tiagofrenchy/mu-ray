import mechanicalsoup as ms
import pandas as pd

data = pd.read_csv('transcriptoma.csv')
ntot = len(data)
new_column = ['' for l in range(len(data))]
for n,row in enumerate(data["mRna - Description"]):
    print(n,'/',ntot)
    if "gene:" in row:
        genecode = row.split("gene:")[1]
        browser = ms.StatefulBrowser()
        browser.open("http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={}".format(genecode))
        page = browser.get_current_page()
        description = page.find("div", class_="rhs")
        if description:
            description2 = description.text.split(" [")[0]
            new_column[n] = description2
        else:
             new_column[n] = 'no description'
    else:
         new_column[n] = 'no genecode'
data.insert(2,"New Description",new_column)
data.to_csv("transcriptoma_descr.csv")
