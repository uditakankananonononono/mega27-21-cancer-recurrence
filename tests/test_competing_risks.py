import numpy as np,sys,os
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_competing_risks import code_competing

CLIN={
 'A':{'RFS_MONTHS':30,'RFS_STATUS':'1:Recurred','OS_MONTHS':50,'OS_STATUS':'1:DECEASED'},
 'B':{'RFS_MONTHS':80,'RFS_STATUS':'0:Not Recurred','OS_MONTHS':60,'OS_STATUS':'1:DECEASED'},
 'C':{'RFS_MONTHS':80,'RFS_STATUS':'0:Not Recurred','OS_MONTHS':90,'OS_STATUS':'0:LIVING'},
 'D':{'RFS_MONTHS':None,'RFS_STATUS':'0:Not Recurred','OS_MONTHS':90,'OS_STATUS':'0:LIVING'},
}

def test_recurrence_first_wins_over_later_death():
    t,ev,keep=code_competing(['A'],CLIN)
    assert ev[0]==1 and t[0]==30

def test_death_before_censoring_is_competing():
    t,ev,keep=code_competing(['B'],CLIN)
    assert ev[0]==2 and t[0]==60

def test_living_patient_censored_at_rfs_time():
    t,ev,keep=code_competing(['C'],CLIN)
    assert ev[0]==0 and t[0]==80

def test_missing_fields_dropped():
    t,ev,keep=code_competing(['D'],CLIN)
    assert keep[0]==False and len(t)==0
