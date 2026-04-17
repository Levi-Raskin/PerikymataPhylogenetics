#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <sstream>
#include <string>
#include <vector>

class Node;
class RandomVariable;

typedef std::map<std::pair<Node*,Node*>,double> BranchLengths;


class Tree {

    public:
                                            Tree(void) = delete;
                                            Tree(const Tree& t); //copy constructor
                                            Tree(std::vector<std::string> taxonNames, double lambda); //random tree w/ taxon names + branches from exp dist
                                            Tree(std::string newick);
                                           ~Tree(void);
        Tree&                               operator=(const Tree& t);
        void                                checkBranchLengthsNeg(void);
        void                                dropTip(std::string tip);
        void                                dropTip(int idx);
        void                                forceBinary(void);
        std::vector<Node*>                  getAllDescendants(Node* n);
        double                              getBranchLength(Node* e1, Node* e2);
        std::vector<Node*>&                 getDownPassSequence(void) { return downPassSequence; }
        std::string                         getNewickString(void);
        int                                 getNumNodes(void) { return (int)nodes.size(); }
        int                                 getNumTaxa(void);
        Node*                               getRoot(void) { return root; }
        Node*                               getTaxonNode(std::string name);
        void                                initializeDownPassSequence(void);
        void                                keepTips(std::vector<std::string> t);
        void                                print(void);
        void                                print(std::string header);
        void                                reindexNodes(void);
        
    private:
        Node*                               addNode(void);
        void                                clone(const Tree& t);
        void                                collapseNode(Node* n);
        void                                deleteNodes(void);
        void                                initializeBranchLengthKey(std::pair<Node*,Node*>& key, Node* e1, Node* e2);
        std::vector<std::string>            parseNewickString(std::string);
        void                                passDown(Node* p, Node* from);
        void                                removeAllBranches(void);
        void                                removeBranch(Node* e1, Node* e2);
        void                                reroot(Node* r);
        double                              roundDecimal(double value, int n);
        void                                setBranch(Node* e1, Node* e2, double x);
        void                                showNode(Node* p, int indent);
        void                                writeTree(Node* p, std::stringstream& strm);
        BranchLengths                       branchLengths;
        std::vector<Node*>                  downPassSequence;
        std::vector<Node*>                  nodes;
        std::vector<Node*>                  tips;
        Node*                               freeNode;
        Node*                               root;
        double                              logLik;
        int                                 numTaxa;
};

#endif
