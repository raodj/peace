#ifndef ZY_Node_HH
#define ZY_Node_HH

class Node {
public:
	std::string sequence; //the bases of the est
	std::string name; //ID of the est,currently it's the starting position of the node.
	std::string comment; //comment to the est in the input est file
	std::vector<int> qualScores;
	std::string direction; // C-complemented(-1), U-uncomplemented(1). only used for ACE format.

	inline Node(std::string n, std::string c, std::string s) {
		name = n;
		comment = c;
		sequence = s;
	}

	inline Node(std::string n, std::string c, std::string s, std::vector<int> q) {
		name = n;
		comment = c;
		sequence = s;
		qualScores = q;
	}

	inline std::string getNodeStr() {
		return sequence;
	}

	/*
	 * get length of the node
	 */
	inline int getLen() {
		return sequence.length();
	}

	/*
	 * get ID of the node
	 */
	inline std::string getName() {
		return name;
	}

	/*
	 * get sequence of the node
	 */
	inline std::string getSeq() {
		return sequence;
	}

	/*
	 * get sequence of the node
	 */
	inline std::string getComment() {
		return comment;
	}

	/*
	 * get quality scores of the node
	 */
	inline std::vector<int> getQualScores() {
		return qualScores;
	}

	inline void setDirection(int i) {
		direction = (i==1? "U" : "C");
	}
};

#endif
