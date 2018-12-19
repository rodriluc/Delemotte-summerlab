from __future__ import division #perform a close approximation of true mathematical division (by pass Zero division error)

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn import metrics
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.datasets import make_classification
from sklearn.datasets import make_blobs
from sklearn.tree import DecisionTreeClassifier
import pylab as pl
from sklearn.preprocessing import StandardScaler

import itertools

#traceback warning errors
import warnings
#warnings.simplefilter("always")
#warnings.warn("Warning...........Message")

path_fv = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/feature_vectors/'
path_pdb = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/PDB_edited/'
path_gpcr = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_gpcr/'
path_ic = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_ionchannel/'
path_enz = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_enzyme/'
path_kin = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_kinase/'


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
#Visualize Decision Tree classifier
def visualize_classifier(model, X, y, ax=None, cmap='rainbow'):
    ax = ax or plt.gca()
    
    # Plot the training points
    ax.scatter(X[:, 0], X[:, 1], c=y, s=30, cmap=cmap,
               clim=(y.min(), y.max()), zorder=3)
    ax.axis('tight')
    ax.axis('off')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # fit the estimator
    model.fit(X, y)
    xx, yy = np.meshgrid(np.linspace(*xlim, num=200),
                         np.linspace(*ylim, num=200))
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)

    # Create a color plot with the results
    n_classes = len(np.unique(y))
    contours = ax.contourf(xx, yy, Z, alpha=0.3,
                           levels=np.arange(n_classes + 1) - 0.5,
                           cmap=cmap, clim=(y.min(), y.max()),
                           zorder=1)

    ax.set(xlim=xlim, ylim=ylim)

#----------------------------------------------------------------------
#Random Forest Classifier
def RF_classifier():

	vec = create_vlist()
	vec = np.asarray(vec)
	#print vec
	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())	

	fam_list = []
	for item in vec.tolist():
		for a in gpcr.tolist():
			if item == a:
				fam_list.append(0) 
		for b in ic.tolist():
			if item == b:
				fam_list.append(1) 
		for c in kin.tolist():
			if item == c:
				fam_list.append(2) 
		for d in enz.tolist():
			if item == d:
				fam_list.append(3) 
	fam_list=np.asarray(fam_list)
	vec = StandardScaler().fit_transform(vec)
	n_samples = vec.shape[0]
	train_ind, test_ind = train_test_split(np.arange(n_samples), test_size=0.50)
	#print train_ind, test_ind
	X_train = vec[train_ind.astype(int),:]
	y_train = fam_list[train_ind.astype(int)]

	X_test = vec[test_ind.astype(int),:]
	y_test = fam_list[test_ind.astype(int)]
	#print y_train, y_test

	#create the classifier and tune the parameters
	clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
	#fit the data
	xtrain = clf.fit(X_train, y_train) 
	clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
	xtest = clf.fit(X_test, y_test)
	feat_train = xtrain.feature_importances_
	feat_test =  xtest.feature_importances_
	RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
		    max_depth=2, max_features='auto', max_leaf_nodes=None,
		    min_impurity_decrease=0.0, min_impurity_split=None,
		    min_samples_leaf=1, min_samples_split=2,
		    min_weight_fraction_leaf=0.0, n_estimators=100, n_jobs=None,
		    oob_score=False, random_state=0, verbose=0, warm_start=False)
	#Line plot of train and test feature importance
	fig,ax = plt.subplots()
	ax.plot(np.arange(feat_train.shape[0]),feat_train)
	ax.plot(np.arange(feat_test.shape[0]), feat_test)

	#Scatter plot with decision boundaries
	'''clf.fit(X_train, y_train)
        score = clf.score(X_test, y_test)

        # Plot the decision boundary. For that, we will assign a color to each
        # point in the mesh [x_min, x_max]x[y_min, y_max].
        if hasattr(clf, "decision_function"):
            Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
        else:
            Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]

        # Put the result into a color plot
        Z = Z.reshape(xx.shape)
        ax.contourf(xx, yy, Z, cmap=cm, alpha=.8)

        # Plot the training points
        ax.scatter(X_train[:, 0], X_train[:, 1], c=y_train, cmap=cm_bright,
                   edgecolors='k')
        # Plot the testing points
        ax.scatter(X_test[:, 0], X_test[:, 1], c=y_test, cmap=cm_bright,
                   edgecolors='k', alpha=0.6)

        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(yy.min(), yy.max())
        ax.set_xticks(())
        ax.set_yticks(())
        if ds_cnt == 0:
            ax.set_title(name)
        ax.text(xx.max() - .3, yy.min() + .3, ('%.2f' % score).lstrip('0'),
                size=15, horizontalalignment='right')
        i += 1'''

	plt.show()

#----------------------------------------------------------------------
#K-fold
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import KFold
def cross_validation():

	vec = create_vlist()
	vec = np.asarray(vec)
	#print vec
	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())	

	fam_list = []
	for item in vec.tolist():
		for a in gpcr.tolist():
			if item == a:
				fam_list.append(0) 
		for b in ic.tolist():
			if item == b:
				fam_list.append(1) 
		for c in kin.tolist():
			if item == c:
				fam_list.append(2) 
		for d in enz.tolist():
			if item == d:
				fam_list.append(3) 
	fam_list=np.asarray(fam_list)
	vec = StandardScaler().fit_transform(vec)
	n_samples = vec.shape[0]

	cv_list = []
	rfk = RepeatedKFold(n_splits=10, n_repeats=3)
	for train_ind, test_ind in rfk.split(vec):
		X_train = vec[train_ind.astype(int),:]
		y_train = fam_list[train_ind.astype(int)]

		X_test = vec[test_ind.astype(int),:]
		y_test = fam_list[test_ind.astype(int)]

		clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
		xtrain = clf.fit(X_train, y_train) 
		feat_train = xtrain.feature_importances_
		cv_list.append(feat_train)

	fig,ax = plt.subplots()
	df = pd.DataFrame(cv_list)
	print df
	ma = df.mean()
	mstd = df.std()
	plt.plot(ma.index, ma)
	plt.fill_between(mstd.index, ma-2*mstd, ma+2*mstd, color='b', alpha=0.2)
	plt.show()
	
#----------------------------------------------------------------------
#NOT BEING USED
def split_train_test(n_splits, samples):
#Split the data into n_splits training and test sets

	if  n_splits < 2:
		logger.info("Using all data in training and validation sets")
		all_indices = np.empty((1, len(samples)))
		for i in range(len(samples)):
			all_indices[0, i] = i
		all_indices = all_indices.astype(int)
		return all_indices, all_indices

	kf = KFold(n_splits=n_splits, shuffle=False)

	train_inds = []
	test_inds = []

	for train_ind, test_ind in kf.split(samples):
		train_inds.append(train_ind)
		test_inds.append(test_ind)
	print train_inds, test_inds

def get_train_test_set(train_ind, test_ind,samples):
# Get the train and test set given their sample/label indices

	train_set = samples[train_ind, :]
	test_set = samples[test_ind, :]

	test_labels = labels[test_ind, :]
	train_labels = labels[train_ind, :]

	return train_set, test_set, train_labels, test_labels

def avg_trainfeat():
	train_inds, test_inds = split_train_test(samples)

	# In for-loop over samples splits:

	train_set, test_set, train_labels, test_labels = get_train_test_set(train_inds[i_split],
                                                                         test_inds[i_split], samples)
	

#----------------------------------------------------------------------
#For CM analysis (below)

from sklearn import model_selection
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from sklearn.linear_model import LogisticRegression

vec = create_vlist()
vec = np.asarray(vec)
gpcr = np.asarray(create_gpcr_list())
ic = np.asarray(create_ic_list())
enz = np.asarray(create_enz_list())
kin = np.asarray(create_kin_list())	

fam_list = []
for item in vec.tolist():
	for a in gpcr.tolist():
		if item == a:
			fam_list.append(0) 
	for b in ic.tolist():
		if item == b:
			fam_list.append(1) 
	for c in kin.tolist():
		if item == c:
			fam_list.append(2) 
	for d in enz.tolist():
		if item == d:
			fam_list.append(3) 
fam_list=np.asarray(fam_list)
vec = StandardScaler().fit_transform(vec)
#print vec
n_samples = vec.shape[0]
'''train_ind, test_ind = train_test_split(np.arange(n_samples), test_size=0.30)
X_train = vec[train_ind.astype(int),:]
y_train = fam_list[train_ind.astype(int)]

X_test = vec[test_ind.astype(int),:]
y_test = fam_list[test_ind.astype(int)]'''

# Mininmized data set with rfk
rfk = RepeatedKFold(n_splits=10, n_repeats=3)
for train_ind, test_ind in rfk.split(vec):
	X_train = vec[train_ind.astype(int),:]
	y_train = fam_list[train_ind.astype(int)]

	X_test = vec[test_ind.astype(int),:]
	y_test = fam_list[test_ind.astype(int)]
#print X_train, y_train
model = LogisticRegression()
model.fit(X_train, y_train)
result = model.score(X_test, y_test)
y_pred = model.predict(X_test)
'''print("Accuracy: %.3f%%" % (result*100.0))
print("F1 Score: ", f1_score(y_test, y_pred, average="macro"))
print("Precision Score: ", precision_score(y_test, y_pred, average="macro"))
print("Recall Score: ", recall_score(y_test, y_pred, average="macro"))'''

def cm_analysis(y_true, y_pred, labels, ymap=None, figsize=(10,10)): 

	if ymap is not None:
		y_pred = [ymap[yi] for yi in y_pred]
		y_true = [ymap[yi] for yi in y_true]
		labels = [ymap[yi] for yi in labels]
	cm = confusion_matrix(y_true, y_pred, labels=labels)
	cm_sum = np.sum(cm, axis=1, keepdims=True)
	cm_perc = cm / cm_sum.astype(float) * 100
	annot = np.empty_like(cm).astype(str)
	nrows, ncols = cm.shape
	for i in range(nrows):
		for j in range(ncols):
			c = cm[i, j]
			p = cm_perc[i, j]
			if i == j:
				s = cm_sum[i]
				annot[i, j] = '%.1f%%\n%d/%d' % (p, c, s)
		    	elif c == 0:
				annot[i, j] = ''
		    	else:
				annot[i, j] = '%.1f%%\n%d' % (p, c)
	labels = 'GPCR','IC','Kinase','Enzyme'
	cm = pd.DataFrame(cm, index=labels, columns=labels)
	#print labels
	cm.index.name = 'Actual'
	cm.columns.name = 'Predicted'
	fig, ax = plt.subplots(figsize=figsize)
	sns.heatmap(cm, annot=annot, fmt='', ax=ax)
	#plt.savefig(filename)
	plt.show()

#cm_analysis(y_test, y_pred, model.classes_, ymap=None, figsize=(10,10))

#----------------------------------------------------------------------
#ROC curve

import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score, recall_score, f1_score

def roc_graph():

	vec = create_vlist()
	vec = np.asarray(vec)
	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())	

	fam_list = []
	for item in vec.tolist():
		for a in gpcr.tolist():
			if item == a:
				fam_list.append(0) 
		for b in ic.tolist():
			if item == b:
				fam_list.append(1) 
		for c in kin.tolist():
			if item == c:
				fam_list.append(2) 
		for d in enz.tolist():
			if item == d:
				fam_list.append(3) 
	fam_list=np.asarray(fam_list)
	vec = StandardScaler().fit_transform(vec)
	n_samples = vec.shape[0]
	#print fam_list

	# Binarize the output
	#y = label_binarize(fam_list, classes=[0, 1, 2, 3])
	#n_classes = y.shape[1]
	
	rfk = RepeatedKFold(n_splits=10, n_repeats=3)
	for train_ind, test_ind in rfk.split(vec):
		X_train = vec[train_ind.astype(int),:]
		y_train = fam_list[train_ind.astype(int)]

		X_test = vec[test_ind.astype(int),:]
		y_test = fam_list[test_ind.astype(int)]

		clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
		xtrain = clf.fit(X_train, y_train) 
		feat_train = xtrain.feature_importances_
		cv_list.append(feat_train)

	# Compute ROC curve and ROC area for each class
	fpr = dict()
	tpr = dict()
	roc_auc = dict()
	for i in range(n_classes):
	    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
	    roc_auc[i] = auc(fpr[i], tpr[i])

	# Compute micro-average ROC curve and ROC area
	fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
	roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

	# Plot of ROC curve for a specific class
	plt.figure()
	lw = 2
	plt.plot(fpr[2], tpr[2], color='darkorange',
		 lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right")
	plt.show()

	y_scores = classifier.decision_function(vec)
	y_pred = classifier.predict(vec)
	print y
	print y_pred
	# ROC metrics
	print 'ROC AUC score:', roc_auc_score(y, y_scores) # AUC from prediction scores
	print 'Avg. precision score:', average_precision_score(y, y_scores) #average precision from prediction scores
	print 'Precision score (micro):', precision_score(fam_list, y_pred, average='micro') #macro/unweighted option, tp/(tp+fp)
	print 'Recall score (micro):', recall_score(fam_list, y_pred, average='micro') #macro/unweighted option, tp/(tp+fn)
	print 'F1 score (weighted):', f1_score(fam_list, y_pred, average='weighted') #macro/micro option, f1=2*(prec*recall)/(prec+recall)

#----------------------------------------------------------------------
#ROC cross-validation

import numpy as np
from scipy import interp
import matplotlib.pyplot as plt

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold

def roc_cv_graph():

	vec = create_vlist()
	vec = np.asarray(vec)
	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())	

	fam_list = []
	for item in vec.tolist():
		for a in gpcr.tolist():
			if item == a:
				fam_list.append(0) 
		for b in ic.tolist():
			if item == b:
				fam_list.append(1) 
		for c in kin.tolist():
			if item == c:
				fam_list.append(2) 
		for d in enz.tolist():
			if item == d:
				fam_list.append(3) 
	fam_list=np.asarray(fam_list)
	vec = StandardScaler().fit_transform(vec)
	n_samples = vec.shape[0]

	# Binarize the output
	y = label_binarize(fam_list, classes=[0, 1, 2, 3])
	n_classes = y.shape[1]
	y = fam_list
	X = vec

	train_ind, test_ind = train_test_split(np.arange(n_samples), test_size=0.30)
	X_train = vec[train_ind.astype(int),:]
	y_train = fam_list[train_ind.astype(int)]
	X_test = vec[test_ind.astype(int),:]
	y_test = fam_list[test_ind.astype(int)]
	y_test = label_binarize(y_test, classes=[0, 1, 2, 3])
	# Learn to predict each class against the other
	classifier = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True))
	y_score = classifier.fit(X_train, y_train).decision_function(X_test)
	#print y_test, y_score

	# Compute ROC curve and ROC area for each class
	fpr = dict()
	tpr = dict()
	roc_auc = dict()
	for i in range(n_classes):
		fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
		roc_auc[i] = auc(fpr[i], tpr[i])

	# Run classifier with cross-validation and plot ROC curves
	cv = StratifiedKFold(n_splits=6)
	classifier = svm.SVC(kernel='linear', probability=True) #OneVsRestClassifier(svm.SVC(kernel='linear', probability=True))


	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)

	i = 0
	for train, test in cv.split(X, y):
		probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
		# Compute ROC curve and area the curve
		print probas_
		fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)
		plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

		i += 1
	plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	plt.plot(mean_fpr, mean_tpr, color='b',
		label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=.8)

	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')

	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")
	plt.show()

#----------------------------------------------------------------------
# MAUC implementaiton from Hand and Till 2001


vec = create_vlist()
vec = np.asarray(vec)
gpcr = np.asarray(create_gpcr_list())
ic = np.asarray(create_ic_list())
enz = np.asarray(create_enz_list())
kin = np.asarray(create_kin_list())	

fam_list = []
for item in vec.tolist():
	for a in gpcr.tolist():
		if item == a:
			fam_list.append(0) 
	for b in ic.tolist():
		if item == b:
			fam_list.append(1) 
	for c in kin.tolist():
		if item == c:
			fam_list.append(2) 
	for d in enz.tolist():
		if item == d:
			fam_list.append(3) 
fam_list=np.asarray(fam_list)
vec = StandardScaler().fit_transform(vec)
n_samples = vec.shape[0]

# Binarize the output -> nope
y = label_binarize(fam_list, classes=[0, 1, 2, 3])
n_classes = y.shape[1]
y = fam_list
X = vec

train_ind, test_ind = train_test_split(np.arange(n_samples), test_size=0.30)
X_train = vec[train_ind.astype(int),:]
y_train = fam_list[train_ind.astype(int)]
X_test = vec[test_ind.astype(int),:]
y_test = fam_list[test_ind.astype(int)]

# Random Forest classifier


rfk = RepeatedKFold(n_splits=10, n_repeats=3)
for train_ind, test_ind in rfk.split(vec):
	X_train = vec[train_ind.astype(int),:]
	y_train = fam_list[train_ind.astype(int)]
	X_test = vec[test_ind.astype(int),:]
	y_test = fam_list[test_ind.astype(int)]
	
	prob_temp = []
	y_list = []
	for item in y:
		y_list.append(item)
		clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
		probas_ = clf.fit(X_train, y_train).predict_proba(X_test)
		prob_temp.append(probas_[item].tolist())
		probabilities = zip(y_list, prob_temp)

# Run classifier with cross-validation 
'''cv = StratifiedKFold(n_splits=10)
classifier = svm.SVC(kernel='linear', probability=True) 

probabilities = []
for train, test in cv.split(X, y):
	prob_temp = []
	y_list = []
	for item in y:
		y_list.append(item)
		probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
		prob_temp.append(probas_[item].tolist())
		#probabilities.append(prob_temp)
		probabilities = zip(y_list, prob_temp)
print y_list
print probabilities'''

def a_value(probabilities, zero_label=0, one_label=1):

	# Obtain a list of the probabilities for the specified zero label class
	expanded_points = []
	for instance in probabilities:
		if instance[0] == zero_label or instance[0] == one_label:
	    		expanded_points.append((instance[0], instance[1][zero_label]))
	#print expanded_points
	sorted_ranks = sorted(expanded_points, key=lambda x: x[1])
			
	n0, n1, sum_ranks = 0, 0, 0
	# Iterate through ranks and increment counters for overall count and ranks of class 0
	for index, point in enumerate(sorted_ranks):
		if point[0] == zero_label:
			n0 += 1
			sum_ranks += index + 1  # Add 1 as ranks are one-based
		elif point[0] == one_label:
			n1 += 1
		else:
			pass  # Not interested in this class
	#try:
	return (float(sum_ranks - (n0*(n0+1)/2.0))) / float(n0 * n1)  # Eqn 3
	#except ZeroDivisionError:
		#return 0

def MAUC(data, num_classes):

	# Find all pairwise comparisons of labels
	class_pairs = [x for x in itertools.combinations(xrange(num_classes), 2)]

	# Have to take average of A value with both classes acting as label 0 as this
	# gives different outputs for more than 2 classes
	sum_avals = 0
	for pairing in class_pairs:
		sum_avals += (a_value(data, zero_label=pairing[0], one_label=pairing[1]) + a_value(data, zero_label=pairing[1], one_label=pairing[0])) / 2.0

	print sum_avals * (2 / float(num_classes * (num_classes-1))) # Eqn 7, return


#----------------------------------------------------------------------
#Precision-Recall curve

#def prec_recall():



if __name__ == '__main__':
	#visualize_classifier()
	#RF_classifier()
	#cross_validation()
	#split_train_test()
	#get_train_test_set
	#cm_analysis(y_test, y_pred, model.classes_, ymap=None, figsize=(10,10))
	#roc_graph()
	#roc_cv_graph()
	#a_value(probabilities)
	MAUC(probabilities, n_classes)
	#prec_recall()




	'''pred_train = clf.predict(X_train)
	pred_test = clf.predict(X_test)
	#print pred_test
	
	# get overall accuracy and F1 score 
	pscore = metrics.accuracy_score(y_test, pred_test)
	#score = metrics.f1_score(y_test, pred_test, pos_label=list(set(y_test)))

	# get size of the full label set
	categories = vec
	dur = len(categories)
	print "Building testing confusion matrix..."
	# initialize score matrices
	trueScores = np.zeros(shape=(dur,dur))
	predScores = np.zeros(shape=(dur,dur))
	# populate totals
	for i in xrange(len(pred_test)-1):
		trueIdx = y_test[i]
		predIdx = pred_test[i]
		trueScores[trueIdx,trueIdx] += 1
		predScores[trueIdx,predIdx] += 1
	# create %-based results
	trueSums = np.sum(trueScores,axis=0)
	#print trueSums
	conf = np.zeros(shape=predScores.shape)
	for i in xrange(len(predScores)):
		for j in xrange(dur):
			if trueSums[i] != 0:
				#print 'yes'
				conf[i,j] = predScores[i,j]/trueSums[i]
			else:
				conf[i,j] == 0
			
	# plot the confusion matrix
	hq = pl.figure(figsize=(15,15));
	aq = hq.add_subplot(1,1,1)
	aq.set_aspect(1)
	res = aq.imshow(conf,cmap=pl.get_cmap('Greens'),interpolation='nearest',vmin=-0.05,vmax=1.)
	width = len(conf)
	height = len(conf[0])
	done = []
	# label each grid cell with the misclassification rates
	for w in xrange(width):
		for h in xrange(height):
			pval = conf[w][h]
			c = 'k'
			rais = w
			if pval > 0.5: c = 'w'
			if pval > 0.001:
				if w == h:
					aq.annotate("{0:1.1f}%\n{1:1.0f}/{2:1.0f}".format(pval*100.,predScores[w][h],trueSums[w]), xy=(h, w), horizontalalignment='center', verticalalignment='center',color=c,size=10)
				else:
					aq.annotate("{0:1.1f}%\n{1:1.0f}".format(pval*100.,predScores[w][h]), xy=(h, w), horizontalalignment='center', verticalalignment='center',color=c,size=10)
	# label the axes
	pl.xticks(range(width), categories[:width],rotation=90,size=10)
	pl.yticks(range(height), categories[:height],size=10)
	# add a title with the F1 score and accuracy
	aq.set_title('Accuracy: '+'{0:2.1f}%'.format(100*pscore)+", " + str(len(y_test)) + "items)",size=10,color='k')
	aq.set_ylabel("Actual",size=10,color='k')
	aq.set_xlabel("Predicted",size=10,color='k')
	pl.grid(b=True,axis='both')

	pl.show()'''





