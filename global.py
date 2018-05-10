from time import time
import numpy as np 

def global_alignment():
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

	mtxScore = init_score_mtx(mtxScore, fInDel)
	mtxBack = init_bt_mtx(mtxBack)

	for j in range(1, nCol):
		for i in range(1, nRow):
			fDiag = mtxScore[i-1, j-1]
			fDiagScore = fMatch if sSeq1[i] == sSeq2[j] else fMismatch
			fDiagScore = fDiag + fDiagScore

			fUp = mtxScore[i-1, j]
			fUpScore = fUp + fInDel

			fLeft = mtxScore[i, j-1]
			fLeftScore = fLeft + fInDel

			lScores = np.array([fDiagScore, fUpScore, fLeftScore])

			mtxScore[i, j] = np.max(lScores)

			for k in np.where(lScores == mtxScore[i, j]):
				mtxBack[i, j, k] = 1

	fMaxScore = mtxScore[-1, -1]

	listPaths = [[(nRow-1, nCol-1), "", ""]]
	listPaths = backtrack(listPaths, mtxBack, sSeq1, sSeq2)
	#print(finalPaths)
	lAlignments = [(path[1], path[2]) for path in listPaths]

	return fMaxScore, lAlignments

def init_score_mtx(mtxScore, fInDel):
	(nM, nN) = mtxScore.shape

	for i in range(1, nM):
		mtxScore[i, 0] = i * fInDel

	for j in range(1, nN):
		mtxScore[0, j] = j * fInDel

	return mtxScore

def init_bt_mtx(mtxBack):
	(nM, nN, _) = mtxBack.shape

	for i in range(1, nM):
		mtxBack[i, 0, 1] = 1

	for j in range(1, nN):
		mtxBack[0, j, 2] = 1

	return mtxBack

def backtrack(listPaths, mtxBack, sSeq1, sSeq2):
	listNewPaths = []
	for path in listPaths:
		(i, j) = path[0]
		pathSeq1 = path[1]
		pathSeq2 = path[2]
		btVal = mtxBack[i, j, :]
		if np.all(btVal==0):
			#print(listPaths)
			return listPaths
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
	fMaxScore, lAlignments = global_alignment()
	print()
	print("Max Score: {}".format(fMaxScore))
	print()
	for alignment in lAlignments:
		sSeq1 = alignment[0]
		sSeq2 = alignment[1]
		print(sSeq1)
		print("|" * len(sSeq1))
		print(sSeq2)
		print()
	rtime = time()-rtime
	print("{} seconds".format(rtime))