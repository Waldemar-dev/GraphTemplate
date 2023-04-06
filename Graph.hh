#ifndef _GRAPH_HH_
#define _GRAPH_HH_

#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <vector>

typedef void (*voidFunctionType)(void);

template <class A>
using voidClassFunctionType = void (A::*)(void);

template <class A>
class Node
{
public:
  // Node() : name_("Default"), externalID_(0) {}
  Node(std::string inName, uint inID) : name_(inName), externalID_(inID) {}
  Node(std::string inName, uint inID, std::string inType, double inValue) : name_(inName), externalID_(inID), data_(std::map<std::string, double>({inType, inValue})) {}
  Node(const Node<A> &inNode) : name_(inNode.getName()), externalID_(inNode.getID()), data_(inNode.getData()), functions_(inNode.getFunctions()), classFunctions_(inNode.getClassFunctions()) {}
  ~Node() = default;

  void addData(std::string inType, double inValue) { data_[inType] = inValue; }
  std::map<std::string, double> getData() const { return data_; }
  std::string getName() const { return name_; }
  uint getID() const { return externalID_; }
  auto getFunctions() const { return functions_; }
  auto getClassFunctions() const { return classFunctions_; }
  template <typename T>
  void insertFunction(std::string s1, T f1);
  template <typename T>
  void insertClassFunction(std::string s1, T f1);
  template <typename T, typename... Args>
  T searchAndCall(std::string s1, Args &&...args);
  template <typename T, typename... Args>
  T searchAndCall(A a, std::string s1, Args &&...args);
  void operator=(const Node<A> &);
  bool operator==(const Node<A> &) const;

private:
  std::string name_;
  uint externalID_;
  std::map<std::string, double> data_;
  std::map<std::string, std::pair<voidFunctionType, std::type_index>> functions_;
  std::map<std::string, std::pair<voidClassFunctionType<A>, std::type_index>> classFunctions_;
};

template <class A>
class Edge
{
public:
  Edge(std::shared_ptr<Node<A>> node1, std::shared_ptr<Node<A>> node2) : nodes_(std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>>(std::move(node1), std::move(node2))) {}
  Edge(const Edge<A> &inEdge) : data_(inEdge.getData()), nodes_(inEdge.getNodes()), functions_(inEdge.getFunctions()), classFunctions_(inEdge.getClassFunctions()) {}
  ~Edge() = default;

  void addData(std::string inType, double inValue) { data_[inType] = inValue; }
  std::map<std::string, double> getData() const { return data_; }
  std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>> getNodes() const { return nodes_; }
  void changeNodes(std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>);
  auto getFunctions() const { return functions_; }
  auto getClassFunctions() const { return classFunctions_; }
  template <typename T>
  void insertFunction(std::string s1, T f1);
  template <typename T>
  void insertClassFunction(std::string s1, T f1);
  template <typename T, typename... Args>
  T searchAndCall(std::string s1, Args &&...args);
  template <typename T, typename... Args>
  T searchAndCall(A a, std::string s1, Args &&...args);
  void operator=(const Edge<A> &);
  bool operator==(const Edge<A> &) const;

private:
  std::map<std::string, double> data_;
  std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>> nodes_;
  std::map<std::string, std::pair<voidFunctionType, std::type_index>> functions_;
  std::map<std::string, std::pair<voidClassFunctionType<A>, std::type_index>> classFunctions_;
};

template <class A>
class BaseGraph
{
public:
  BaseGraph() {}
  BaseGraph(const BaseGraph &in) : nodes_(in.getNodes()), edges_(in.getEdges()), endNodes_(in.getEndNodes()), adjacencyMap_(in.getAdjacencyMap()) {}
  ~BaseGraph() = default;

  std::vector<Edge<A>> getEdges() const { return edges_; }
  std::vector<Node<A>> getNodes() const { return nodes_; }
  std::vector<std::shared_ptr<Node<A>>> getEndNodes() const { return endNodes_; }
  std::multimap<uint, std::pair<uint, uint>> getAdjacencyMap() const { return adjacencyMap_; }
  void add(Node<A>, Edge<A>, Node<A>);
  void clear();
  void operator=(const BaseGraph<A> &);
  bool operator==(const BaseGraph<A> &) const;

protected:
  void findEndNodes();
  void copy(const BaseGraph *);
  std::vector<Edge<A>> edges_;
  std::vector<Node<A>> nodes_;
  std::vector<std::shared_ptr<Node<A>>> endNodes_;          // can be empty
  std::multimap<uint, std::pair<uint, uint>> adjacencyMap_; // adjacency matrix without the zeros
};

template <class A>
class Path : public BaseGraph<A>
{
public:
  Path() : BaseGraph<A>(), updateName_("Update") {}
  Path(const Path &in) : BaseGraph<A>(in), updateName_(in.getUpdateName()) {}
  ~Path() = default;

  Node<A> getRootNode() const;
  void addToEndNode(Node<A> inNode, Path<A> inGraph);
  void setUpdateName(std::string inName) { updateName_ = inName; }
  double getProbability(A);
  template <typename B>
  double getSumOfValues(std::string, A);
  std::string getUpdateName() const { return updateName_; }

  void clear();
  void operator=(const Path<A> &);
  bool operator==(const Path<A> &) const;

private:
  std::string updateName_;
};

template <class A>
class DecisionTree : public BaseGraph<A>
{
public:
  DecisionTree() : BaseGraph<A>() {}
  DecisionTree(const DecisionTree &in) : BaseGraph<A>(in), allPaths_(in.getAllPaths()) {}
  ~DecisionTree() = default;

  Node<A> getRootNode() const;
  void addToEndNode(Node<A> inNode, DecisionTree<A> inTree, bool refreshEndNodes = true);
  void connectGraphs(std::vector<Node<A>>, DecisionTree<A>);
  std::vector<Path<A>> getAllPaths() const { return allPaths_; }
  uint getNAllPaths() const { return allPaths_.size(); }
  void findAllPaths();

  void clear();
  void operator=(const DecisionTree<A> &);
  bool operator==(const DecisionTree<A> &) const;

private:
  void DFS(uint start, uint finish, Path<A> &);
  void DFS(Node<A> start, Node<A> finish, Path<A> *);
  std::vector<Path<A>> allPaths_;
};

// definitions for Node
template <class A>
template <typename T>
void Node<A>::insertFunction(std::string s1, T f1)
{
  auto tt = std::type_index(typeid(f1));
  functions_.insert(std::make_pair(s1, std::make_pair((voidFunctionType)f1, tt)));
}

template <class A>
template <typename T>
void Node<A>::insertClassFunction(std::string s1, T f1)
{
  auto tt = std::type_index(typeid(f1));
  classFunctions_.insert(std::make_pair((std::string)s1, std::make_pair((voidClassFunctionType<A>)f1, tt)));
}

template <class A>
template <typename T, typename... Args>
T Node<A>::searchAndCall(std::string s1, Args &&...args)
{
  auto mapIter = functions_.find(s1);
  auto mapVal = mapIter->second;

  auto typeCastedFun = (T(*)(Args...))(mapVal.first);

  // compare if the types are equal or not
  assert(mapVal.second == std::type_index(typeid(typeCastedFun)));
  return typeCastedFun(std::forward<Args>(args)...);
}

template <class A>
template <typename T, typename... Args>
T Node<A>::searchAndCall(A a, std::string s1, Args &&...args)
{
  auto mapIter = classFunctions_.find(s1);
  auto mapVal = mapIter->second;

  auto typeCastedFun = (T(A::*)(Args...))(mapVal.first);

  assert(mapVal.second == std::type_index(typeid(typeCastedFun)));
  return (a.*typeCastedFun)(std::forward<Args>(args)...);
}

template <class A>
void Node<A>::operator=(const Node<A> &inNode)
{
  name_ = inNode.getName();
  externalID_ = inNode.getID();
  data_ = inNode.getData();
  functions_ = inNode.getFunctions();
  classFunctions_ = inNode.getClassFunctions();
}

template <class A>
bool Node<A>::operator==(const Node<A> &inNode) const
{
  bool names = (name_ == inNode.getName());
  bool ids = (externalID_ == inNode.getID());
  bool functions = (functions_ == inNode.getFunctions());
  bool classFunctions = (classFunctions_ == inNode.getClassFunctions());
  bool data = (data_ == inNode.getData());
  if (names && ids && data && functions && classFunctions)
  {
    return true;
  }
  return false;
}

// definitions for Edge
template <class A>
void Edge<A>::changeNodes(std::shared_ptr<Node<A>> node1, std::shared_ptr<Node<A>> node2)
{
  nodes_.first = std::move(node1);
  nodes_.second = std::move(node2);
}

template <class A>
template <typename T>
void Edge<A>::insertFunction(std::string s1, T f1)
{
  auto tt = std::type_index(typeid(f1));
  functions_.insert(std::make_pair(s1, std::make_pair((voidFunctionType)f1, tt)));
}

template <class A>
template <typename T>
void Edge<A>::insertClassFunction(std::string s1, T f1)
{
  auto tt = std::type_index(typeid(f1));
  classFunctions_.insert(std::make_pair((std::string)s1, std::make_pair((voidClassFunctionType<A>)f1, tt)));
}

template <class A>
template <typename T, typename... Args>
T Edge<A>::searchAndCall(std::string s1, Args &&...args)
{
  auto mapIter = functions_.find(s1);
  auto mapVal = mapIter->second;

  auto typeCastedFun = (T(*)(Args...))(mapVal.first);

  // compare if the types are equal or not
  assert(mapVal.second == std::type_index(typeid(typeCastedFun)));
  return typeCastedFun(std::forward<Args>(args)...);
}

template <class A>
template <typename T, typename... Args>
T Edge<A>::searchAndCall(A a, std::string s1, Args &&...args)
{
  auto mapIter = classFunctions_.find(s1);
  auto mapVal = mapIter->second;

  auto typeCastedFun = (T(A::*)(Args...))(mapVal.first);

  assert(mapVal.second == std::type_index(typeid(typeCastedFun)));
  return (a.*typeCastedFun)(std::forward<Args>(args)...);
}

template <class A>
void Edge<A>::operator=(const Edge<A> &inEdge)
{
  functions_ = inEdge.getFunctions();
  classFunctions_ = inEdge.getClassFunctions();
  data_ = inEdge.getData();
  nodes_ = inEdge.getNodes();
}

template <class A>
bool Edge<A>::operator==(const Edge<A> &inEdge) const
{
  bool nodes = (nodes_ == inEdge.getNodes());
  bool data = (data_ == inEdge.getData());
  bool functions = (functions_ == inEdge.getFunctions());
  bool classFunctions = (classFunctions_ == inEdge.getClassFunctions());
  if (nodes && data && functions && classFunctions)
  {
    return true;
  }
  return false;
}

// definitions for BaseGraph
template <class A>
void BaseGraph<A>::add(Node<A> inNode1, Edge<A> inEdge, Node<A> inNode2)
{
  bool addNode1 = true;
  bool addNode2 = true;
  std::shared_ptr<Node<A>> node1;
  std::shared_ptr<Node<A>> node2;
  for (uint i = 0; i < nodes_.size(); i++)
  {
    if (nodes_[i] == inNode1)
    {
      addNode1 = false;
      node1 = std::make_shared<Node<A>>(nodes_[i]);
    }
    else if (nodes_[i] == inNode2)
    {
      addNode2 = false;
      node2 = std::make_shared<Node<A>>(nodes_[i]);
    }
    if (!addNode1 && !addNode2)
    {
      break;
    }
  }
  if (addNode1)
  {
    inNode1.addData("ID", nodes_.size());
    nodes_.push_back(inNode1);
    node1 = std::make_shared<Node<A>>(nodes_.back());
  }
  if (addNode2)
  {
    inNode2.addData("ID", nodes_.size());
    nodes_.push_back(inNode2);
    node2 = std::make_shared<Node<A>>(nodes_.back());
  }
  Edge<A> inEdgeCopy(inEdge);
  inEdgeCopy.changeNodes(node1, node2);
  std::map<std::string, double> tempData(inEdge.getData());
  for (std::map<std::string, double>::iterator it = tempData.begin();
       it != tempData.end(); it++)
  {
    inEdgeCopy.addData(it->first, it->second);
  }
  edges_.push_back(inEdgeCopy);
}

template <class A>
void BaseGraph<A>::findEndNodes()
{
  std::vector<std::shared_ptr<Node<A>>> endPoints{
      std::make_shared<Node<A>>(nodes_[0])};
  for (uint i = 0; i < endPoints.size(); i++)
  {
    if (BaseGraph<A>::adjacencyMap_.find(endPoints[i]->getData()["ID"]) != adjacencyMap_.end())
    {
      std::pair<std::multimap<uint, std::pair<uint, uint>>::iterator, std::multimap<uint, std::pair<uint, uint>>::iterator>
          range = adjacencyMap_.equal_range(endPoints[i]->getData()["ID"]);
      for (auto it = range.first; it != range.second; it++)
      {
        endPoints.push_back(std::make_shared<Node<A>>(nodes_[it->second.first]));
      }
      if (i == 0)
      {
        endPoints.erase(endPoints.begin());
      }
      else
      {
        endPoints.erase(endPoints.begin() + i);
      }
    }
  }
  endNodes_ = endPoints;
}

template <class A>
void BaseGraph<A>::copy(const BaseGraph<A> *in)
{
  edges_ = in->getEdges();
  nodes_ = in->getNodes();
  endNodes_ = in->getEndNodes();
  adjacencyMap_ = in->getAdjacencyMap();
}

template <class A>
void BaseGraph<A>::clear()
{
  nodes_.clear();
  edges_.clear();
  endNodes_.clear();
  adjacencyMap_.clear();
}

template <class A>
void BaseGraph<A>::operator=(const BaseGraph<A> &in)
{
  nodes_ = in.getNodes();
  edges_ = in.getEdges();
  endNodes_ = in.getEndNodes();
  adjacencyMap_ = in.getAdjacencyMap();
}

template <class A>
bool BaseGraph<A>::operator==(const BaseGraph<A> &in) const // To DO: sort nodes and edges
{
  bool nodes = (nodes_ == in.getNodes());
  bool edges = (edges_ == in.getEdges());
  bool endNodes = (endNodes_ == in.getEndNodes());
  bool adMap = (adjacencyMap_ == in.getAdjacencyMap());
  if (nodes && edges && endNodes && adMap)
  {
    return true;
  }
  return false;
}

// definitions for DecisionTree
template <class A>
void DecisionTree<A>::clear()
{
  BaseGraph<A>::clear();
  allPaths_.clear();
}

template <class A>
void DecisionTree<A>::DFS(uint startID, uint finishID, Path<A> &result)
{
  if (startID == finishID)
  {
    allPaths_.push_back(result);
  }
  else
  {
    if (BaseGraph<A>::adjacencyMap_.find(startID) != BaseGraph<A>::adjacencyMap_.end())
    {
      std::pair<std::multimap<uint, std::pair<uint, uint>>::iterator,
                std::multimap<uint, std::pair<uint, uint>>::iterator>
          range = BaseGraph<A>::adjacencyMap_.equal_range(startID);
      for (std::multimap<uint, std::pair<uint, uint>>::iterator it =
               range.first;
           it != range.second; it++)
      {
        Path<A> tempResult(result);
        tempResult.add(std::make_shared<Node<A>>(BaseGraph<A>::nodes_[it->first]),
                       std::make_shared<Edge<A>>(BaseGraph<A>::edges_[it->second.second]),
                       std::make_shared<Node<A>>(BaseGraph<A>::nodes_[it->second.first]));
        DFS(it->second.first, finishID, tempResult);
      }
    }
  }
}

template <class A>
void DecisionTree<A>::DFS(Node<A> startNode, Node<A> finishNode, Path<A> *result)
{
  if (startNode == finishNode)
  {
    allPaths_.push_back(*result);
  }
  else
  {
    uint startID = 0;
    bool foundStartID = false;
    for (uint i = 0; i < BaseGraph<A>::nodes_.size(); i++)
    {
      if (BaseGraph<A>::nodes_[i] == startNode)
      {
        startID = i;
        foundStartID = true;
        break;
      }
    }
    if (!foundStartID)
    {
      std::cout << "DFS cannot find start node." << std::endl;
      abort();
    }
    std::pair<std::multimap<uint, std::pair<uint, uint>>::iterator,
              std::multimap<uint, std::pair<uint, uint>>::iterator>
        range = BaseGraph<A>::adjacencyMap_.equal_range(startID);
    for (std::multimap<uint, std::pair<uint, uint>>::iterator it =
             range.first;
         it != range.second; it++)
    {
      Path<A> tempResult(*result);
      tempResult.add(BaseGraph<A>::nodes_[it->first],
                     BaseGraph<A>::edges_[it->second.second],
                     BaseGraph<A>::nodes_[it->second.first]);
      DFS(BaseGraph<A>::nodes_[it->second.first], finishNode, &tempResult);
    }
  }
}

template <class A>
void DecisionTree<A>::addToEndNode(Node<A> inNode, DecisionTree<A> inGraph, bool refreshEndNodes)
{
  bool foundNode = false;
  uint iNode = 0;
  for (uint i = 0; i < BaseGraph<A>::endNodes_.size(); i++)
  {
    if ((*BaseGraph<A>::endNodes_[i]) == inNode)
    {
      foundNode = true;
      iNode = i;
      break;
    }
  }
  if (!foundNode)
  {
    std::cout << "Cannot connect two graphs." << std::endl;
    abort();
  }
  else
  {
    uint nodesSize = BaseGraph<A>::nodes_.size();
    std::vector<Node<A>> newNodes(inGraph.getNodes());
    std::vector<Edge<A>> newEdges(inGraph.getEdges());
    Edge<A> link(BaseGraph<A>::endNodes_[iNode], std::make_shared<Node<A>>(newNodes[0])); // TO DO:get root node
    BaseGraph<A>::add(*BaseGraph<A>::endNodes_.at(iNode), link, newNodes.at(0));
    std::multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getAdjacencyMap());
    for (std::multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it)
    {
      BaseGraph<A>::add(newNodes[it->first], newEdges[it->second.second], newNodes[it->second.first]);
    }
  }
  if (refreshEndNodes)
  {
    BaseGraph<A>::endNodes_.clear();
    BaseGraph<A>::findEndNodes();
  }
}

template <class A>
void DecisionTree<A>::connectGraphs(std::vector<Node<A>> inNode, DecisionTree<A> inGraph)
{
  bool foundNode = false;
  std::vector<uint> iNodes;
  for (uint i = 0; i < BaseGraph<A>::endNodes_.size(); i++)
  {
    if (std::find(inNode.begin(), inNode.end(), (*BaseGraph<A>::endNodes_[i])) != inNode.end())
    {
      foundNode = true;
      iNodes.push_back(i);
    }
  }
  if (!foundNode)
  {
    std::cout << "Cannot connect two graphs." << std::endl;
    abort();
  }
  else
  {
    uint nodesSize = BaseGraph<A>::nodes_.size();
    std::vector<Node<A>> newNodes(inGraph.getNodes());
    std::vector<Edge<A>> newEdges(inGraph.getEdges());
    for (uint iNode : iNodes)
    {
      Edge<A> link(BaseGraph<A>::endNodes_[iNode], std::make_shared<Node<A>>(newNodes[0])); // TO DO: get root node
      BaseGraph<A>::add(*BaseGraph<A>::endNodes_[iNode], link, newNodes[0]);
    }
    std::multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getAdjacencyMap());
    for (std::multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it++)
    {
      BaseGraph<A>::add(newNodes[it->first], newEdges[it->second.second], newNodes[it->second.first]);
    }
  }
  BaseGraph<A>::endNodes_.clear();
  BaseGraph<A>::findEndNodes();
}

template <class A>
void DecisionTree<A>::operator=(const DecisionTree<A> &in)
{
  BaseGraph<A>::operator=(in);
  allPaths_ = in.getAllPaths();
}

template <class A>
bool DecisionTree<A>::operator==(const DecisionTree<A> &in) const
{
  bool base = (BaseGraph<A>::operator==(in));
  bool paths = (allPaths_ == in.getAllPaths());
  if (base && paths)
  {
    return true;
  }
  return false;
}

template <class A>
void DecisionTree<A>::findAllPaths()
{
  for (std::shared_ptr<Node<A>> endNode : BaseGraph<A>::endNodes_)
  {
    Path<A> newPath;
    DFS(getRootNode(), *endNode, &newPath);
    allPaths_.push_back(newPath);
  }
}

template <class A>
Node<A> DecisionTree<A>::getRootNode() const
{
  uint lastNodeIndex = 0;
  auto it = BaseGraph<A>::adjacencyMap_.begin();
  uint safetyCounter = 0;
  while (it != BaseGraph<A>::adjacencyMap_.end() && safetyCounter < BaseGraph<A>::adjacencyMap_.size())
  {
    if (it->second.first == lastNodeIndex)
    {
      lastNodeIndex = it->first;
      it = BaseGraph<A>::adjacencyMap_.begin();
      safetyCounter++;
    }
    else
    {
      it++;
    }
  }
  return BaseGraph<A>::nodes_[lastNodeIndex];
}

// definitions of Path
template <class A>
double Path<A>::getProbability(A inObject)
{
  std::string functionName = "Probability";
  double result = 1;
  for (auto it = BaseGraph<A>::adjacencyMap_.begin(); it != BaseGraph<A>::adjacencyMap_.end(); it++)
  {
    BaseGraph<A>::edges_[it->second.second].template searchAndCall<void>(inObject, updateName_);

    if (BaseGraph<A>::edges_[it->second.second].getClassFunctions().find(functionName) !=
        BaseGraph<A>::edges_[it->second.second].getClassFunctions().end())
    {
      result *= BaseGraph<A>::edges_[it->second.second].template searchAndCall<double>(inObject, functionName);
    }
    else if (BaseGraph<A>::edges_[it->second.second].getFunctions().find(functionName) !=
             BaseGraph<A>::edges_[it->second.second].getFunctions().end())
    {
      result *= BaseGraph<A>::edges_[it->second.second].template searchAndCall<double>(functionName);
    }
    else if (BaseGraph<A>::edges_[it->second.second].getData().find(functionName) !=
             BaseGraph<A>::edges_[it->second.second].getData().end())
    {
      result *= BaseGraph<A>::edges_[it->second.second].getData()[functionName];
    }
    else
    {
      result = 0;
      return result;
    }
  }
  return result;
}

template <class A>
template <typename B>
double Path<A>::getSumOfValues(std::string valueName, A inObject)
{
  double result = 0;
  const Node<A> rootNode(getRootNode());
  uint rootNodeIndex = 0;
  for (uint i = 0; i < BaseGraph<A>::nodes_.size(); i++)
  {
    if (BaseGraph<A>::nodes_[i] == rootNode)
    {
      rootNodeIndex = i;
      break;
    }
  }
  bool done = false;
  uint safetyCounter = 0;
  uint lastIndex = rootNodeIndex;
  while (!done && safetyCounter < BaseGraph<A>::nodes_.size())
  {
    safetyCounter++;
    if (BaseGraph<A>::nodes_[lastIndex].getClassFunctions().find(valueName) != BaseGraph<A>::nodes_[lastIndex].getClassFunctions().end())
    {
      result += BaseGraph<A>::nodes_.at(lastIndex).template searchAndCall<B>(inObject, valueName);
    }
    else if (BaseGraph<A>::nodes_[lastIndex].getFunctions().find(valueName) != BaseGraph<A>::nodes_[lastIndex].getFunctions().end())
    {
      result += BaseGraph<A>::nodes_[lastIndex].template searchAndCall<B>(valueName);
    }
    else if (BaseGraph<A>::nodes_[lastIndex].getData().find(valueName) != BaseGraph<A>::nodes_[lastIndex].getData().end())
    {
      result += BaseGraph<A>::nodes_[lastIndex].getData()[valueName];
    }
    auto it = BaseGraph<A>::adjacencyMap_.find(lastIndex);
    uint edgeIndex = it->second.second;
    BaseGraph<A>::edges_[edgeIndex].template searchAndCall<void>(inObject, updateName_);
    if (BaseGraph<A>::edges_[edgeIndex].getClassFunctions().find(valueName) != BaseGraph<A>::edges_[edgeIndex].getClassFunctions().end())
    {
      result += BaseGraph<A>::edges_[edgeIndex].template searchAndCall<B>(inObject, valueName);
    }
    else if (BaseGraph<A>::edges_[edgeIndex].getFunctions().find(valueName) != BaseGraph<A>::edges_[edgeIndex].getFunctions().end())
    {
      result += BaseGraph<A>::edges_[edgeIndex].template searchAndCall<B>(valueName);
    }
    else if (BaseGraph<A>::edges_[edgeIndex].getData().find(valueName) != BaseGraph<A>::edges_[edgeIndex].getData().end())
    {
      result += BaseGraph<A>::edges_[edgeIndex].getData()[valueName];
    }
    if (BaseGraph<A>::adjacencyMap_.find(lastIndex) == BaseGraph<A>::adjacencyMap_.end())
    {
      done = true;
    }
    else
    {
      lastIndex = BaseGraph<A>::adjacencyMap_.find(lastIndex)->second.first;
    }
  }
  return result;
}

template <class A>
void Path<A>::operator=(const Path<A> &in)
{
  BaseGraph<A>::operator=(in);
  updateName_ = in.getUpdateName();
}

template <class A>
bool Path<A>::operator==(const Path<A> &in) const
{
  bool base = (BaseGraph<A>::operator==(in));
  bool name = (updateName_ == in.getUpdateName());
  if (base && name)
  {
    return true;
  }
  return false;
}

template <class A>
Node<A> Path<A>::getRootNode() const
{
  uint lastNodeIndex = 0;
  auto it = BaseGraph<A>::adjacencyMap_.begin();
  uint safetyCounter = 0;
  while (it != BaseGraph<A>::adjacencyMap_.end() && safetyCounter < BaseGraph<A>::adjacencyMap_.size())
  {
    if (it->second.first == lastNodeIndex)
    {
      lastNodeIndex = it->first;
      it = BaseGraph<A>::adjacencyMap_.begin();
      safetyCounter++;
    }
    else
    {
      it++;
    }
  }
  return BaseGraph<A>::nodes_[lastNodeIndex];
}
#endif