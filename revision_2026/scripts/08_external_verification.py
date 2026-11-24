# Cross-check case files against independently published checkpoints (all cited in out/08 CSV).
# File values: EDCD district-by-month workbook (2019-2023; not public), preprint table (2024). Months treated as Gregorian.
import pandas as pd
p=pd.read_csv('out/paper_72month_table.csv',index_col=0,parse_dates=True).cases
cum=lambda y,m: int(p[str(y)][:m].sum())
rows=[
 # (year, check, file value, published value, source, expected relation)
 (2019,'Annual total',int(p['2019'].sum()),17992,'Rijal et al. 2021 Infect Dis Poverty; WHO Nepal 2022 gives 17,999','equal'),
 (2019,'Kathmandu district annual',1870,1463,'Pandey et al. 2022 Epidemiol Infect','equal'),
 (2019,'Sunsari district (to ~21 Sep)',122,3455,'Kathmandu Post 21 Sep 2019 (district health official); 3,266 by 13 Aug (Republica)','file >= published'),
 (2020,'Annual total',int(p['2020'].sum()),530,'Pandey et al. 2022 Epidemiol Infect','equal'),
 (2021,'Annual total',int(p['2021'].sum()),540,'Pandey et al. 2022 Epidemiol Infect','equal'),
 (2021,'Kathmandu 2020+2021',66,63,'Pandey et al. 2022 Epidemiol Infect','equal'),
 (2022,'Jan–Aug cumulative',cum(2022,8),28109-23514,'WHO DON412: 28,109 to 28 Sep, of which 23,514 in September','equal'),
 (2022,'Jan–28 Sep cumulative',cum(2022,8)+int(p['2022-09'].iloc[0]*28/30),28109,'WHO DON412 (real-time)','file >= published'),
 (2022,'October monthly',int(p['2022-10'].iloc[0]),17889,'EDCD via Khabarhub (2023): "17,889 cases recorded in October" 2022','equal'),
 (2022,'Jan–24 Oct cumulative',cum(2022,9)+int(p['2022-10'].iloc[0]*24/31),42504,'WHO Nepal 24 Oct 2022 (real-time)','file >= published'),
 (2022,'Jan–20 Nov cumulative',cum(2022,10)+int(p['2022-11'].iloc[0]*20/30),52557,'HRW citing government data, 20 Nov 2022','file >= published'),
 (2022,'Jan–30 Nov cumulative',cum(2022,11),53951,'Dhoubhadel et al. 2023 Front Trop Dis','file >= published'),
 (2022,'Annual total',int(p['2022'].sum()),54784,'EDCD via Viruses 2023 (Rimal et al.), Gavi','equal'),
 (2022,'Bagmati province annual',42428,42143,'MoHP 11 Dec 2022 via PMC9886556 (real-time)','file >= published'),
 (2023,'Jan–30 Jul cumulative',cum(2023,6)+int(p['2023-07'].iloc[0]*30/31),5688,'EDCD situation update 30 Jul 2023 (ReliefWeb)','file >= published'),
 (2023,'Jan–4 Sep cumulative',cum(2023,8)+int(p['2023-09'].iloc[0]*4/30),23127,'MoHP via fitfortravel 5 Sep 2023','file >= published'),
 (2023,'Jan–27 Sep cumulative',cum(2023,8)+int(p['2023-09'].iloc[0]*27/30),32897,'WHO via ReliefWeb','file >= published'),
 (2023,'Sunsari Jan–Sep',16066,16036,'WHO via ReliefWeb, to 27 Sep 2023','file >= published'),
 (2023,'August monthly',int(p['2023-08'].iloc[0]),13846,'EDCD via Khabarhub, as of 11 Oct 2023 (real-time)','file >= published'),
 (2023,'September monthly',int(p['2023-09'].iloc[0]),13472,'WHO via ReliefWeb, as of 17 Nov 2023 (11,580 on 11 Oct)','file >= published'),
 (2023,'Jan–17 Nov cumulative',cum(2023,10)+int(p['2023-11'].iloc[0]*17/30),49761,'WHO via ReliefWeb','file >= published'),
 (2023,'June–December total',int(p['2023-06':'2023-12'].sum()),51243,'IFRC DREF final report: 51,243 cases June–December 2023 (EDCD situation report)','file >= published'),
 (2024,'Annual total',int(p['2024'].sum()),41865,'Kathmandu Post 24 Feb, 29 Apr and 12 Oct 2025 (EDCD figures)','equal'),
]
t=pd.DataFrame(rows,columns=['year','check','file','published','source','expected'])
t['diff_%']=((t.file-t.published)/t.published*100).round(1)
def verdict(r):
    if r.expected=='equal': return 'MATCH' if abs(r['diff_%'])<=0.5 else ('CLOSE' if abs(r['diff_%'])<=3 else 'MISMATCH')
    return 'CONSISTENT' if r.file>=r.published*0.98 else 'INCONSISTENT'
t['verdict']=t.apply(verdict,axis=1); print(t.drop(columns='source').to_string(index=False)); t.to_csv('out/08_external_verification.csv',index=False)
