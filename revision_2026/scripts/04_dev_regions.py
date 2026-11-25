# Map EDCD districts (77 + "Not specified") to Nepal's pre-2015 five development regions,
# which is the sub-national unit Google Trends uses for Nepal.
import openpyxl, pandas as pd, sys
R={'Eastern':'TAPLEJUNG PANCHTHAR ILAM JHAPA SANKHUWASABHA TERHATHUM DHANKUTA BHOJPUR MORANG SUNSARI SOLUKHUMBU OKHALDHUNGA KHOTANG UDAYAPUR SAPTARI SIRAHA',
   'Central':'DOLAKHA RAMECHHAP SINDHULI DHANUSHA DHANUSA MAHOTTARI SARLAHI KATHMANDU LALITPUR BHAKTAPUR KAVREPALANCHOK KABHREPALANCHOK SINDHUPALCHOK SINDHUPALCHOWK NUWAKOT RASUWA DHADING MAKWANPUR MAKAWANPUR CHITAWAN CHITWAN PARSA BARA RAUTAHAT',
   'Western':'GORKHA LAMJUNG TANAHU TANAHUN SYANGJA KASKI MANANG MUSTANG MYAGDI BAGLUNG PARBAT RUPANDEHI KAPILBASTU KAPILVASTU PALPA ARGHAKHANCHI GULMI NAWALPARASI',
   'Mid-Western':'PYUTHAN ROLPA RUKUM SALYAN DANG BANKE BARDIYA SURKHET DAILEKH JAJARKOT DOLPA HUMLA JUMLA KALIKOT MUGU',
   'Far-Western':'BAJURA BAJHANG ACHHAM DOTI KAILALI KANCHANPUR DADELDHURA BAITADI DARCHULA'}
def region(name):
    n=name.upper()
    for r,ds in R.items():
        if any(d in n for d in ds.split()): return r
    return None
wb=openpyxl.load_workbook(sys.argv[1],data_only=True); rows=[]
for ws in wb:
    y=int(ws.title[:4])
    for r in list(ws.iter_rows(values_only=True))[2:]:
        if r[1] in (None,'Total'): continue
        for i,v in enumerate(r[3:15]): rows.append((pd.Timestamp(y,i+1,1),r[1],region(r[1]) if r[1]!='Not specified' else 'Not specified',v or 0))
d=pd.DataFrame(rows,columns=['date','district','region','cases'])
print('unmapped:',sorted(d[d.region.isna()].district.unique()))
print('districts per region:',d.drop_duplicates('district').groupby('region').size().to_dict())
t=d.groupby(['date','region']).cases.sum().unstack(); t.to_csv('out/edcd_devregion_monthly_2019_2023.csv')
print((t.groupby(t.index.year).sum().T/ t.groupby(t.index.year).sum().sum(1)*100).round(1))
