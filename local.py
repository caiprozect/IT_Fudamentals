from time import time
import numpy as np 

def local_alignment():
	sSeq1 = input("Sequence 1: ")
	sSeq2 = input("Sequence 2: ")

	sSeq1 = "N" + sSeq1.strip().upper()
	sSeq2 = "N" + sSeq2.strip().upper()

	fMatch = input("Match score: ")
	fMismatch = input("Mismatch score: ")
	fInDel = input("Indel score: ")

	fMatch = float(fMatch)
	fMismatch = float(fMismatch)
	fInDel = float(fInDel) 

	nRow = len(sSeq1)
	nCol = len(sSeq2)

	mtxScore = np.zeros((nRow, nCol))
	mtxBack = np.zeros((nRow, nCol, 3))

	for j in range(1, nCol):
		for i in range(1, nRow):
			fDiag = mtxScore[i-1, j-1]
			fDiagScore = fMatch if sSeq1[i] == sSeq2[j] else fMismatch
			fDiagScore = fDiag + fDiagScore

			fUp = mtxScore[i-1, j]
			fUpScore = fUp + fInDel

			fLeft = mtxScore[i, j-1]
			fLeftScore = fLeft + fInDel

			lScores = np.array([fDiagScore, fUpScore, fLeftScore, 0])

			mtxScore[i, j] = np.max(lScores)

			if mtxScore[i, j] != 0:
				for k in np.where(lScores == mtxScore[i, j]):
					mtxBack[i, j, k] = 1

	fMaxScore = np.max(mtxScore)

	argMaxs = np.where(mtxScore == fMaxScore)

	xs = argMaxs[0]
	ys = argMaxs[1]

	listPos = zip(xs, ys)

	listPaths = [[pos, "", ""] for pos in listPos]
	listPaths = prep_paths(listPaths, sSeq1, sSeq2)
	listPaths = backtrack(listPaths, mtxBack, sSeq1, sSeq2)
	#print(finalPaths)
	lAlignments = set([(path[1] + "," + path[2]) for path in listPaths])
	lAlignments = [path.split(',') for path in lAlignments]

	return fMaxScore, lAlignments

def prep_paths(listPaths, sSeq1, sSeq2):
	listPrepPaths = []
	for path in listPaths:
		(i, j) = path[0]
		pathSeq1 = path[1]
		pathSeq2 = path[2]
		while True:
			i = min([i+1, len(sSeq1)])
			j = min([j+1, len(sSeq2)])
			if i == len(sSeq1) and j == len(sSeq2):
				break
			if i < len(sSeq1):
				pathSeq1 = pathSeq1 + sSeq1[i]
			else:
				pathSeq1 = pathSeq1 + "-"
			if j < len(sSeq2):
				pathSeq2 = pathSeq2 + sSeq2[j]
			else:
				pathSeq2 = pathSeq2 + "-"
		listPrepPaths.append([path[0], pathSeq1, pathSeq2])
	return listPrepPaths

def backtrack(listPaths, mtxBack, sSeq1, sSeq2):
	listNewPaths = []
	listPos = [path[0] for path in listPaths]
	listVal = [mtxBack[i, j, :] for (i, j) in listPos]
	listCheck = [np.all(btVal==0) for btVal in listVal]
	if np.all(listCheck):
		return listPaths
	for path in listPaths:
		(i, j) = path[0]
		pathSeq1 = path[1]
		pathSeq2 = path[2]
		btVal = mtxBack[i, j, :]
		if np.all(btVal==0):
			listNewPaths.append(path)
		else:
			if btVal[0] == 1:
				listNewPaths.append([(i-1, j-1), sSeq1[i] + pathSeq1, sSeq2[j] + pathSeq2])
			if btVal[1] == 1:
				listNewPaths.append([(i-1, j), sSeq1[i] + pathSeq1, "-" + pathSeq2])
			if btVal[2] == 1:
				listNewPaths.append([(i, j-1), "-" + pathSeq1, sSeq2[j] + pathSeq2])
	del listPaths
	#print(listNewPaths)
	return backtrack(listNewPaths, mtxBack, sSeq1, sSeq2)

if __name__ == "__main__":
	rtime = time()
	fMaxScore, lAlignments = local_alignment()
	print()
	print("Max Score: {}".format(fMaxScore))
	print()
	for alignment in lAlignments:
		sSeq1 = alignment[0]
		sSeq2 = alignment[1]
		print(sSeq1)
		print("|" * min([len(sSeq1.strip('-')), len(sSeq2.strip('-'))]))
		print(sSeq2)
		print()
	rtime = time()-rtime
	print("{} seconds".format(rtime))