#ifndef ZY_SingleBase_HH
#define ZY_SingleBase_HH

/*
 * record one base in the consensus
 */
class SingleBase {
public:
	char curBase;
	int curBaseNum;
	int curqualVal;

	int numA;
	int numC;
	int numG;
	int numT;
	int numN;
	int numDash;
	int qualValA;
	int qualValC;
	int qualValG;
	int qualValT;
	int qualValN;
	int qualValDash;

	SingleBase(char c, int qualVal=0);

	SingleBase(char c1, char c2, int qualVal1=0, int qualVal2=0);

	void addOneBase(char base, int qualVal=0);

	char getCurBase();

	int getTotalNumOfBase(); //return the total number of ACGTN at the current base position

	int getQualScore(); //return Phred score of current base

private:
	void changeBaseNum(char base, int qualVal=0);
	void changeCurBase();
};

#endif
