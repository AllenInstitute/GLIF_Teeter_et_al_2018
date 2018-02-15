import numpy as np
import allensdk.ephys.ephys_features as ft

def estimate_dv_cutoff(voltage_list, dt, start_t, end_t):
    v_set = [ v * 1e3 for v in voltage_list ]
    t_set = [ np.arange(0, len(v)) * dt for v in voltage_list ]
    
    dv_cutoff, thresh_frac = ft.estimate_adjusted_detection_parameters(v_set, t_set,
                                                                       start_t, end_t,
                                                                       filter=None)
    
    return dv_cutoff, thresh_frac