import argparse
import sklearn.metrics

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--data',required=True,type=str,help='Tab-delimited file with first column is score and second column is 0/1 for negative/positive')
	parser.add_argument('--classbalance',default=0.5,type=float,help='Class balance for calculating precision metric')
	parser.add_argument('--prCurveData',required=False,type=str,help='Output file to give precision recall curve data')
	args = parser.parse_args()

	reweight = args.classbalance

	scoreData = []
	with open(args.data) as f:
		for line in f:
			score,posOrNeg = line.strip().split('\t')
			isPos = (posOrNeg == '1')
			scoreData.append((score,isPos))

	posCount = sum( 1 for _,isPos in scoreData if isPos )
	negCount = len(scoreData) - posCount

	scoreData = sorted(scoreData,reverse=True)
	TP,FP = 0,0
	bestFScore = -1.0

	curPoints = []

	for score,isPos in scoreData:
		if isPos:
			TP += 1
		else:
			FP += 1

		TN = negCount - FP
		FN = posCount - TP

		precision,recall,fscore = 0,0,0
		if TP+FP != 0:
			precision = reweight*TP / float(reweight*TP + (1-reweight)*FP)
		if TP+FN != 0:
			recall = TP / float(TP+FN)
		if precision+recall != 0:
			fscore = 2 * (precision*recall) / (precision+recall)

		curPoints.append((recall,precision))
		#print(score,TP,FP,TN,FN,precision,recall,fscore)
		if fscore > bestFScore:
			bestFScore = fscore
			#print(TP,FP,TN,FN,precision,recall,fscore)

	curPoints = sorted(curPoints, reverse=True)
				
	# Add the graph points in bottom left and right
	curPoints = curPoints + [(0,1)]

	# Pull out recall and precision points separately (for numpy call)
	recalls = [ r for (r,_) in curPoints ]
	precisions = [ p for (_,p) in curPoints ]

	# Calculate the area using the trapezium rule	
	areaUnderPRCurve = sklearn.metrics.auc(recalls, precisions)
			
	print(areaUnderPRCurve)

	if args.prCurveData:
		with open(args.prCurveData,'w') as f:
			f.write("recall\tprecision\n")
			for r,p in curPoints:
				f.write("%f\t%f\n" % (r,p))

