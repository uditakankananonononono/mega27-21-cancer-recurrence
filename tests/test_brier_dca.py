import numpy as np,sys,os
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_brier_dca import ipcw_brier,decision_curves

def test_no_censoring_ipcW_equals_plain_mse():
    t=np.array([10.,20.,70.,80.]);e=np.array([1,1,0,0]);p=np.array([.4,.6,.3,.1])
    b,bn,p0,nk,nc,w,y=ipcw_brier(t,e,p)
    plain=float(np.mean((np.array([1,1,0,0])-p)**2))
    assert abs(b-plain)<1e-9 and nk==4 and nc==0

def test_treat_none_line_and_weights_shape():
    t=np.array([10.,70.]);e=np.array([1,0]);p=np.array([.9,.1])
    b,bn,p0,nk,nc,w,y=ipcw_brier(t,e,p)
    rows=decision_curves(t,e,p,w,y)
    assert len(rows)==12 and all(r['treat_none']==0.0 for r in rows)
    assert all('net_benefit' in r and 'treat_all' in r for r in rows)

def test_event_before_horizon_only_counts_when_observed():
    t=np.array([30.,70.]);e=np.array([0,0])  # one censored at 30, one at risk past 60
    b,bn,p0,nk,nc,w,y=ipcw_brier(t,e,np.array([.5,.5]))
    assert nc==1 and nk==1
