'''Check to see if any two python files have different fields or values within
a specified tolerance.
'''

import dictdiffer
import allensdk.core.json_utilities as ju
import os
relative_path=os.path.dirname(os.getcwd())

json_file1=relative_path+'/mouse_struc_data_dir/313860745_Rorb-IRES2-Cre-D/313860745_Rorb-IRES2-Cre-D_GLIF3_exp_var_ratio_10ms.json'
json_file2=relative_path+'/create_data_dir/mouse_struc_data_dir/313860745_Rorb-IRES2-Cre-D/313860745_Rorb-IRES2-Cre-D_GLIF3_exp_var_ratio_10ms.json'

# toy to see how dictdiffer works
a={'a': .01, 'b':.02, 'c':.031, 'd':.04, 'f':{'f1':.01, 'f2':.021, 'f3':.030}, 'g':['a', .010, 'c']}
b={'a': .01, 'b':.02, 'c':.030, 'e':.05, 'f':{'f1':.01, 'f2':.020, 'f3':.030}, 'g':['a', .011, 'c']}
out=list(dictdiffer.diff(a,b, tolerance=.001))
print list(out)

d1=ju.read(json_file1)
d2=ju.read(json_file2)

result=list(dictdiffer.diff(d1,d2))
print result

