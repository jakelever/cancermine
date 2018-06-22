import argparse
import random

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_curve, auc

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Fit a linear model for final Cancermine results')
	parser.add_argument('--annotationData',type=str,required=True,help='TSV of data from annotation solicitation')
	args = parser.parse_args()

	X,y = [],[]

	with open(args.annotationData) as f:
		for line in f:
			yesOrNo,role,gene_id,cancer_id,sentenceCount,maxpredictprob,inference_score = line.strip().split('\t')
			sentenceCount = float(sentenceCount)
			maxpredictprob = float(maxpredictprob)
			inference_score = float(inference_score)
			if yesOrNo in ['y','n']:
				y.append(1 if yesOrNo == 'y' else 0)
				#X.append([maxpredictprob])
				#X.append([sentenceCount,maxpredictprob])#,inference_score])
				#X.append([sentenceCount,maxpredictprob,inference_score])
				X.append([sentenceCount,maxpredictprob])

	X = np.array(X)
	#y = np.array(y)

	random.seed(1)
	#trainIndices = random.sample(list(range(len(X))), int(len(X)/2))
	#testIndices
	#trainX = [ X[i] for i in trainIndices
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1)

	#clf = LogisticRegression()
	clf = LinearSVC()
	clf.fit(X_train,y_train)
	y_score = clf.decision_function(X_test)

	fpr, tpr, _ = roc_curve(y_test, y_score)
	roc_auc = auc(fpr, tpr)

	print(roc_auc)

	#print(X_test)
	#y_scores = clf.predict(X_test)
	#print(dir(clf))
	#print(clf.coef_)
	#print(y_test)
	#print(y_score)
