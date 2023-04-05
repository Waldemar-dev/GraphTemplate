#ifndef _GRAPH_HH_
#define _GRAPH_HH_

#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <algorithm>

template <class A>
class Path;

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
  bool operator==(const Node<A> &);

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
  bool operator==(const Edge<A> &);

private:
  std::map<std::string, double> data_;
  std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>> nodes_;
  std::map<std::string, std::pair<voidFunctionType, std::type_index>> functions_;
  std::map<std::string, std::pair<voidClassFunctionType<A>, std::type_index>> classFunctions_;
};

template <class A>
class Graph
{
public:
  Graph() {}
  Graph(const Graph &in) : nodes_(in.getNodes()), edges_(in.getEdges()), endNodes_(in.getEndNodes()), allPaths_(in.getAllPaths()) {}
  ~Graph() = default;
  const Node<A> *getRootNode() const { return &nodes_[0]; } // TO DO
  void add(Node<A>, Edge<A>, Node<A>);
  void addToEndNode(Node<A> inNode, Graph<A> inGraph, bool refreshEndNodes = true);
  void connectGraphs(std::vector<Node<A>>, Graph<A>);
  std::vector<Edge<A>> getEdges() const { return edges_; }
  std::vector<Node<A>> getNodes() const { return nodes_; }
  std::vector<std::shared_ptr<Node<A>>> getEndNodes() const { return endNodes_; }
  std::vector<Path<A>> getAllPaths() const { return allPaths_; }
  uint getNAllPaths() const { return allPaths_.size(); }
  void findAllPaths();
  void clear();
  std::unordered_multimap<uint, std::pair<uint, uint>> getGraphMap() const { return graphMap_; }
  void operator=(const Graph<A> &);

private:
  void findEndNodes();
  void DFS(uint start, uint finish, Path<A> &);
  void DFS(std::shared_ptr<Node<A>> start, std::shared_ptr<Node<A>> finish, Path<A> &);
  std::vector<Node<A>> nodes_;
  std::vector<Edge<A>> edges_;
  std::vector<std::shared_ptr<Node<A>>> endNodes_;
  std::unordered_multimap<uint, std::pair<uint, uint>> graphMap_;
  std::vector<Path<A>> allPaths_;
};

template <class A>
class Path : public Graph<A>
{
public:
  Path();
  Path(const Path &);
  ~Path() = default;

  // void add(std::shared_ptr<Node<A>>, std::shared_ptr<Edge<A>>, std::shared_ptr<Node<A>>);
  void setUpdateName(std::string inName) { updateName_ = inName; }
  double getProbability(A);
  template <typename B>
  double getSumOfValues(std::string, A) const;
  std::string getUpdateName() const { return updateName_; }
  std::vector<std::shared_ptr<Node<A>>> getPathNodes() const { return pathNodes_; }
  std::vector<std::shared_ptr<Edge<A>>> getPathEdges() const { return pathEdges_; }
  std::unordered_multimap<uint, std::pair<uint, uint>> getPathMap() const { return pathMap_; }

private:
  std::vector<std::shared_ptr<Node<A>>> pathNodes_;
  std::vector<std::shared_ptr<Edge<A>>> pathEdges_;
  std::unordered_multimap<uint, std::pair<uint, uint>> pathMap_;
  std::string updateName_;
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

  // auto typeCastedFun = reinterpret_cast<T(*)(Args ...)>(mapVal.first);
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
bool Node<A>::operator==(const Node<A> &inNode)
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

  // auto typeCastedFun = reinterpret_cast<T(*)(Args ...)>(mapVal.first);
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
bool Edge<A>::operator==(const Edge<A> &inEdge)
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

// definitions for Graph
template <class A>
void Graph<A>::add(Node<A> inNode1, Edge<A> inEdge, Node<A> inNode2)
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
void Graph<A>::findEndNodes()
{
  std::vector<std::shared_ptr<Node<A>>> endPoints{
      std::make_shared<Node<A>>(nodes_[0])};
  for (uint i = 0; i < endPoints.size(); i++)
  {
    if (graphMap_.find(endPoints[i]->getData()["ID"]) != graphMap_.end())
    {
      std::pair<
          std::multimap<
              uint, std::pair<uint, std::map<std::string, double>>>::iterator,
          std::multimap<
              uint, std::pair<uint, std::map<std::string, double>>>::iterator>
          range = graphMap_.equal_range(endPoints[i]->getData()["ID"]);
      for (std::multimap<uint, std::pair<uint, std::map<std::string, double>>>::
               iterator it = range.first;
           it != range.second; it++)
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
void Graph<A>::DFS(uint startID, uint finishID, Path<A> &result)
{
  if (startID == finishID)
  {
    allPaths_.push_back(result);
  }
  else
  {
    if (graphMap_.find(startID) != graphMap_.end())
    {
      std::pair<std::unordered_multimap<uint, std::pair<uint, uint>>::iterator,
                std::unordered_multimap<uint, std::pair<uint, uint>>::iterator>
          range = graphMap_.equal_range(startID);
      for (std::unordered_multimap<uint, std::pair<uint, uint>>::iterator it =
               range.first;
           it != range.second; it++)
      {
        Path<A> tempResult(result);
        tempResult.add(std::make_shared<Node<A>>(nodes_[it->first]),
                       std::make_shared<Edge<A>>(edges_[it->second.second]),
                       std::make_shared<Node<A>>(nodes_[it->second.first]));
        DFS(it->second.first, finishID, tempResult);
      }
    }
  }
}

template <class A>
void Graph<A>::DFS(std::shared_ptr<Node<A>> startNode, std::shared_ptr<Node<A>> finishNode, Path<A> &result)
{
  if (startNode == finishNode)
  {
    allPaths_.push_back(result);
  }
  else
  {
    uint startID = 0;
    bool foundStartID = false;
    for (uint i = 0; i < nodes_.size(); i++)
    {
      if (nodes_[i] = (*startNode))
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
    std::pair<std::unordered_multimap<uint, std::pair<uint, uint>>::iterator,
              std::unordered_multimap<uint, std::pair<uint, uint>>::iterator>
        range = graphMap_.equal_range(startID);
    for (std::unordered_multimap<uint, std::pair<uint, uint>>::iterator it =
             range.first;
         it != range.second; it++)
    {
      Path<A> tempResult(result);
      tempResult.add(std::make_shared<Node<A>>(nodes_[it->first]),
                     std::make_shared<Edge<A>>(edges_[it->second.second]),
                     std::make_shared<Node<A>>(nodes_[it->second.first]));
      DFS(std::make_shared<Node<A>>(nodes_[it->second.first]), finishNode, tempResult);
    }
  }
}

template <class A>
void Graph<A>::addToEndNode(Node<A> inNode, Graph<A> inGraph, bool refreshEndNodes)
{
  bool foundNode = false;
  uint iNode = 0;
  for (uint i = 0; i < endNodes_.size(); i++)
  {
    if ((*endNodes_[i]) == inNode)
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
    uint nodesSize = nodes_.size();
    std::vector<Node<A>> newNodes(inGraph.getNodes());
    std::vector<Edge<A>> newEdges(inGraph.getEdges());
    Edge<A> link(std::make_shared<Node<A>>(endNodes_[iNode]), std::make_shared<Node<A>>(newNodes[0]));
    add(endNodes_[iNode], link, newNodes[0]);
    std::unordered_multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getGraphMap());
    for (std::unordered_multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it)
    {
      add(newNodes[it->first], newEdges[it->second.second], newNodes[it->second.first]);
    }
  }
  if (refreshEndNodes)
  {
    endNodes_.clear();
    findEndNodes();
  }
}

template <class A>
void Graph<A>::connectGraphs(std::vector<Node<A>> inNode, Graph<A> inGraph)
{
  bool foundNode = false;
  std::vector<uint> iNodes;
  for (uint i = 0; i < endNodes_.size(); i++)
  {
    if (std::find(inNode.begin(), inNode.end(), (*endNodes_[i])) != inNode.end())
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
    uint nodesSize = nodes_.size();
    std::vector<Node<A>> newNodes(inGraph.getNodes());
    std::vector<Edge<A>> newEdges(inGraph.getEdges());
    for (uint iNode : iNodes)
    {
      Edge<A> link(endNodes_[iNode], std::make_shared<Node<A>>(newNodes[0]));
      add(endNodes_[iNode], link, newNodes[0]);
    }
    std::unordered_multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getGraphMap());
    for (std::unordered_multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it++)
    {
      add(newNodes[it->first], newEdges[it->second.second], newNodes[it->second.first]);
    }
  }
  endNodes_.clear();
  findEndNodes();
}

template <class A>
void Graph<A>::operator=(const Graph<A> &in)
{
  nodes_ = in.getNodes();
  edges_ = in.getEdges();
  graphMap_ = in.getGraphMap();
  allPaths_ = in.getAllPaths();
  findEndNodes();
}

template <class A>
void Graph<A>::clear()
{
  nodes_.clear();
  edges_.clear();
  endNodes_.clear();
  graphMap_.clear();
  allPaths_.clear();
}

template <class A>
void Graph<A>::findAllPaths()
{
  for (std::shared_ptr<Node<A>> endNode : endNodes_)
  {
    Path<A> newPath;
    DFS(std::make_shared<Node<A>>(nodes_[0]), std::make_shared<Node<A>>(endNode), &newPath);
  }
}

// definitions of Path
template <class A>
Path<A>::Path() : updateName_("Update") {}

template <class A>
Path<A>::Path(const Path &inPath) : pathNodes_(inPath.getPathNodes()), pathEdges_(inPath.getPathEdges()), pathMap_(inPath.getPathMap()), updateName_(inPath.getUpdateName()) {}

// template <class A>
// void Path<A>::add(std::shared_ptr<Node<A>> inNode1, std::shared_ptr<Edge<A>> inEdge1,
//                   std::shared_ptr<Node<A>> inNode2)
// {
//   bool addNode1 = true;
//   bool addNode2 = true;
//   std::shared_ptr<Node<A>> node1 = 0;
//   std::shared_ptr<Node<A>> node2 = 0;
//   uint idNode1, idNode2, idEdge;
//   for (uint i = 0; i < pathNodes_.size(); i++)
//   {
//     if (pathNodes_[i] == inNode1)
//     {
//       addNode1 = false;
//       idNode1 = i;
//     }
//     else if (pathNodes_[i] == inNode2)
//     {
//       addNode2 = false;
//       idNode2 = i;
//     }
//     if (!addNode1 && !addNode2)
//     {
//       break;
//     }
//   }
//   if (addNode1)
//   {
//     pathNodes_.push_back(inNode1);
//   }
//   if (addNode2)
//   {
//     pathNodes_.push_back(inNode2);
//   }
//   pathEdges_.push_back(inEdge1);
//   idEdge = pathEdges_.size() - 1;
//   std::pair<uint, uint> node2EdgePair = {idNode2, idEdge};
//   pathMap_.insert(
//       std::pair<uint, std::pair<uint, uint>>(idNode1, node2EdgePair));
// }

template <class A>
double Path<A>::getProbability(A inObject)
{
  std::string functionName = "Probability";
  double result = 1;
  for (std::unordered_multimap<uint, std::pair<uint, uint>>::const_iterator it =
           pathMap_.begin();
       it != pathMap_.end(); it++)
  {
    pathEdges_[it->second.second]->template searchAndCall<void>(inObject, updateName_);

    if (pathEdges_[it->second.second]->getClassFunctions().find(functionName) !=
        pathEdges_[it->second.second]->getClassFunctions().end())
    {
      result *= pathEdges_[it->second.second]->template searchAndCall<double>(inObject, functionName);
    }
    else if (pathEdges_[it->second.second]->getFunctions().find(functionName) !=
             pathEdges_[it->second.second]->getFunctions().end())
    {
      result *= pathEdges_[it->second.second]->template searchAndCall<double>(functionName);
    }
    else if (pathEdges_[it->second.second]->getData().find(functionName) !=
             pathEdges_[it->second.second]->getData().end())
    {
      result *= pathEdges_[it->second.second]->getData()[functionName];
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
double Path<A>::getSumOfValues(std::string valueName, A inObject) const
{
  double result = 0;
  for (std::shared_ptr<Node<A>> node : pathNodes_)
  {
    if (node->getClassFunctions().find(valueName) != node->getClassFunctions().end())
    {
      result += node->template searchAndCall<B>(inObject, valueName);
    }
    else if (node->getFunctions().find(valueName) != node->getFunctions().end())
    {
      result += node->template searchAndCall<B>(valueName);
    }
    else if (node->getData().find(valueName) != node->getData().end())
    {
      result += node->getData()[valueName];
    }
  }
  return result;
}
#endif