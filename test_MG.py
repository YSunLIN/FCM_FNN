# encoding: utf-8
import fcm_fnn
import pickle


def loadDataSet(filepath):
    with open(filepath) as f:
        rawList = map(lambda line: float(line.strip()), f.readlines())
        labelSet = []
        testSet = []
        for i in range(8, 8 + 500):
            labelSet.append((rawList[i:i+4], rawList[i+1:i+5]))
        for i in range(8 + 500, 8 + 1000):
            testSet.append((rawList[i:i+4], rawList[i+1:i+5]))
        return labelSet, testSet


if __name__ == "__main__":
    labelSet, testSet = loadDataSet("./data/mg30.dat")

    with open("fnn.bin", "rb") as f:
        fnn = pickle.loads(f.read())

    # fnn = fcm_fnn.FCM_FNN(3, 3, 3, 3)

    for i in range(25):
        print "train err: %s   test err: %s" % (fnn.train(labelSet, 0.02), fnn.test(testSet))

        for oj, oconcept in enumerate(fnn.concepts):
            for omj in range(oconcept.numOfTerms):
                print "(%s, %s) C:%s W:%s Xi:%s" % (oj, omj, oconcept.C[omj], oconcept.sigma[omj], oconcept.xi[omj])
            print
    
    print "write fnn data"
    with open("fnn.bin", "wb") as f:
        f.write(pickle.dumps(fnn))

    print "write test files"
    testData = []
    for X, D in testSet:
        testData.append(D[0])
    testData += testSet[-1][1][1:4]
    with open("testData.csv", "w") as f:
        f.write("\n".join(map(lambda v:str(v), testData)))

    print "write predict files"
    predictData = []
    for X, D in testSet:
        predictData.append(fnn.predict(X)[0])
    predictData += fnn.predict(testSet[-1][0])[1:4]
    with open("predictData.csv", "w") as f:
        f.write("\n".join(map(lambda v:str(v), predictData)))
