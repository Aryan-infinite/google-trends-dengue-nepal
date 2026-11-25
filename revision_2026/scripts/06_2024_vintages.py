# Compare the paper's final 2024 monthly counts (EDCD website) with cumulative totals
# reported in real time during 2024 (news citing EDCD). Real-time figures come from dated web reports.
import pandas as pd
p=pd.read_csv('out/paper_72month_table.csv',index_col=0,parse_dates=True).cases['2024']
cum=p.cumsum()
def final_by(date):  # final count up to a date, linear within month
    d=pd.Timestamp(date); m=d.to_period('M').to_timestamp(); prev=cum[cum.index<m].iloc[-1] if (cum.index<m).any() else 0
    frac=d.day/d.days_in_month; return prev+frac*p[m]
rt=[('2024-07-31',2000,'Nepalnews 1 Aug 2024 (">2,000")'),('2024-09-03',6000,'Nepalnews 3 Sep 2024 (">6,000")'),
    ('2024-09-23',12644,'Himalayan Times 26 Sep 2024'),('2024-10-15',19599,'Onlinekhabar/Republica 16 Oct 2024'),
    ('2024-10-25',27153,'Republica 28 Oct 2024'),('2024-12-31',41865,'Kathmandu Post 2025 (annual)')]
rows=[(d,r,round(final_by(d)),round(r/final_by(d)*100),s) for d,r,s in rt]
t=pd.DataFrame(rows,columns=['as_of','reported_then','final_count_by_then','pct_reported','source']); print(t.to_string(index=False))
t.to_csv('out/06_2024_reporting_lag.csv',index=False)
print('paper 2024 cumulative by month end:',cum.astype(int).tolist())
