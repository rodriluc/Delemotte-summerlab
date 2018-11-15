import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.datasets import make_classification

path_fv = '/data2/LucieR/Delemotte-summerlab_ERnorm_4.5A/feature_vectors/'
path_pdb = '/data2/LucieR/Delemotte-summerlab_ERnorm_4.5A/PDB_edited/'
path_gpcr = '/data2/LucieR/Delemotte-summerlab_ERnorm_4.5A/fv_gpcr/'
path_ic = '/data2/LucieR/Delemotte-summerlab_ERnorm_4.5A/fv_ionchannel/'
path_enz = '/data2/LucieR/Delemotte-summerlab_ERnorm_4.5A/fv_enzyme/'
path_kin = '/data2/LucieR/Delemotte-summerlab_ERnorm_4.5A/fv_kinase/'


#----------------------------------------------------------------------
#Necessary data
def name_base(path_pdb,file):
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			base_list.append(base)
	return base_list

def gpcr_base(path_gpcr,file):
	base_list = []
	for file in os.listdir(path_gpcr):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list

def ic_base(path_ic,file):
	base_list = []
	for file in os.listdir(path_ic):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list	

def enz_base(path_enz,file):
	base_list = []
	for file in os.listdir(path_enz):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list

def kin_base(path_kin,file):
	base_list = []
	for file in os.listdir(path_kin):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list
		
def create_vlist():
	base_list = name_base(path_pdb,file)
	temp = []
	for i in range(len(base_list)):
		with open (path_fv+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]
			x = np.array(fv)
			temp.append(x)	
	return temp
	
def create_gpcr_list():
	base_list = gpcr_base(path_gpcr,file)
	temp = []
	for i in range(len(base_list)):
		with open (path_gpcr+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]
			x = np.array(fv)
			temp.append(x)
	return temp

def create_ic_list():
	base_list = ic_base(path_ic,file)
	temp = []
	for i in range(len(base_list)):
		with open (path_ic+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]
			x = np.array(fv)
			temp.append(x)
	return temp

def create_enz_list():
	base_list = enz_base(path_enz,file)
	temp = []
	for i in range(len(base_list)):
		with open (path_enz+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]
			x = np.array(fv)
			temp.append(x)
	return temp

def create_kin_list():
	base_list = kin_base(path_kin,file)
	temp = []
	for i in range(len(base_list)):
		with open (path_kin+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]
			x = np.array(fv)
			temp.append(x)
	return temp
#----------------------------------------------------------------------



X, y = make_classification(n_samples=388, n_features=17, n_informative=2, n_redundant=0, random_state=0, shuffle=False)
#create the classifier and tune the parameters
clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
#fit the data
clf.fit(X, y) #(train, targets_train)
RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
            max_depth=2, max_features='auto', max_leaf_nodes=None,
            min_impurity_decrease=0.0, min_impurity_split=None,
            min_samples_leaf=1, min_samples_split=2,
            min_weight_fraction_leaf=0.0, n_estimators=100, n_jobs=None,
            oob_score=False, random_state=0, verbose=0, warm_start=False)

print(clf.feature_importances_)
[0.14205973 0.76664038 0.0282433  0.06305659]
print(clf.predict([[0, 0, 0, 0]]))
[1]
