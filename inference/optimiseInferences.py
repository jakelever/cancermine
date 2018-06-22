import argparse
import logging
from scipy.io.matlab import loadmat
from scipy.sparse import lil_matrix
from rescal import rescal_als
from numpy import dot,zeros
import itertools
import random
import sys
from collections import defaultdict,Counter
import scipy.spatial.distance

def predict_rescal_als(T,rank,lambda_A,lambda_R):
	A, R, _, _, _ = rescal_als(
		T, rank, init='nvecs', conv=1e-4,
		lambda_A=lambda_A, lambda_R=lambda_R
	)
	n = A.shape[0]
	P = zeros((n, n, len(R)))
	for k in range(len(R)):
		P[:, :, k] = dot(A, dot(R[k], A.T))
	return P

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Make inferences on CancerMine using RESCAL')
	parser.add_argument('--trainingTensor',required=True,type=str,help='Tensor to use to train inferences')
	parser.add_argument('--evalTensor',required=True,type=str,help='Tensor to use for evaluation')
	parser.add_argument('--rank',default=10,type=int,help='Parameter for RESCAL')
	parser.add_argument('--lambda_A',default=0.0,type=float,help='Parameter for RESCAL')
	parser.add_argument('--lambda_R',default=0.0,type=float,help='Parameter for RESCAL')
	parser.add_argument('--outScores',required=True,type=str,help='Output test scores with class information for evaluation')
	args = parser.parse_args()

	with open(args.trainingTensor) as inF:
		trainDim = list(map(int,inF.readline().strip(' #').split()))
		dimx,dimy,dimz = trainDim
		training = [ lil_matrix((dimy,dimz)) for _ in range(dimx) ]
		for line in inF:
			x,y,z,val = line.strip().split()
			x,y,z = int(x),int(y),int(z)
			val = float(val)
			training[x][y,z] = val
			
	print("Running RESCAL with rank=%d, lambda_A=%f and lambda_R=%f" % (args.rank,args.lambda_A,args.lambda_R))
	preds = predict_rescal_als(training, args.rank, args.lambda_A, args.lambda_R)

	with open(args.evalTensor) as inF, open(args.outScores,'w') as outF:
		evalDim = list(map(int,inF.readline().strip(' #').split()))
		assert trainDim == evalDim, "%s != %s" % (trainDim,evalDim)
		for line in inF:
			x,y,z,c = line.strip().split()
			x,y,z = int(x),int(y),int(z)
			c = int(c)
			score = preds[y,z,x]
			outF.write("%f\t%d\n" % (score,c))

