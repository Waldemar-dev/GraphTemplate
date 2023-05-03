#ifndef _GRAPH_HH_
#define _GRAPH_HH_

#include <cassert>
#include <map>
#include <set>
#include <memory>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <vector>
#include <tuple>
#include <functional>
#include <iostream>
#include "/home/main_user/Tools/invoke.hpp-main/headers/invoke.hpp/invoke.hpp" // TO DO: version check

typedef void (*voidFunctionType)(void);

template <class A>
using voidClassFunctionType = void (A::*)(void);

template <class A>
class Node
{
public:
  Node(std::string inName, uint inID) : name_(inName), externalID_(inID) {}
  Node(const Node<A> &inNode) : name_(inNode.getName()), externalID_(inNode.getID()), data_(inNode.getData()), functions_(inNode.getFunctions()), classFunctions_(inNode.getClassFunctions()) {}
  ~Node() = default;

  void addData(std::string inType, double inValue) { data_[inType] = inValue; }
  void overwriteData(std::map<std::string, double> inData) { data_ = inData; }
  std::map<std::string, double> getData() const { return data_; }
  double getData(std::string inName) const { return data_.at(inName); }
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
  void print() const;
  void operator=(const Node<A> &);
  bool operator==(const Node<A> &) const;

private:
  std::string name_;
  uint externalID_;
  std::map<std::string, double> data_;
  std::map<std::string, std::pair<voidFunctionType, std::type_index>> functions_;
  std::map<std::string, std::pair<voidClassFunctionType<A>, std::type_index>> classFunctions_;
};

template <class A, typename... functionArgs>
class Edge
{
public:
  Edge(std::shared_ptr<Node<A>> node1, std::shared_ptr<Node<A>> node2) : nodes_(std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>>(std::move(node1), std::move(node2))) {}
  Edge(const Edge<A> &inEdge) : data_(inEdge.getData()), nodes_(inEdge.getNodes()), functions_(inEdge.getFunctions()), classFunctions_(inEdge.getClassFunctions()), linkages_(inEdge.getLinkages()) /*, argsMap_(inEdge.getArgsMap())*/ {}
  ~Edge() = default;

  void addData(std::string inType, double inValue) { data_[inType] = inValue; }
  void defineLinkage(std::string node1FunctionName, std::string edgeFunctionName, std::function<double(double, double)> linkage) { linkages_.insert(std::pair<std::string, std::pair<std::string, std::function<double(double, double)>>>(node1FunctionName, std::pair<std::string, std::function<double(double, double)>>(edgeFunctionName, linkage))); }
  std::map<std::string, std::pair<std::string, std::function<double(double, double)>>> getLinkages() const { return linkages_; }
  std::map<std::string, double> getData() const { return data_; }
  template <typename... Args>
  void updateData(A, Args &&...args);
  void updateData(A);
  template <typename... Args>
  void updateData(Args &&...args);
  std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>> getNodes() const { return nodes_; }
  void changeNodes(std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>);
  auto getFunctions() const { return functions_; }
  auto getClassFunctions() const { return classFunctions_; }
  template <typename T>
  void insertFunction(std::string s1, T f1);
  template <typename T>
  void insertClassFunctionVoid(std::string s1, T f1);
  template <typename T>
  void insertClassFunction(std::string, T, functionArgs &&...);
  template <typename T, typename... Args>
  T searchAndCall(std::string s1, Args &&...args);
  template <typename T, typename... Args>
  T searchAndCall(A a, std::string s1, Args &&...args);
  std::map<std::string, std::tuple<functionArgs...>> getArgsMap() const { return argsMap_; }
  void operator=(const Edge<A> &);
  bool operator==(const Edge<A> &) const;

private:
  std::map<std::string, double> data_;
  std::pair<std::shared_ptr<Node<A>>, std::shared_ptr<Node<A>>> nodes_;
  std::map<std::string, std::pair<voidFunctionType, std::type_index>> functions_;
  std::map<std::string, std::pair<voidClassFunctionType<A>, std::type_index>> classFunctions_;
  std::map<std::string, std::pair<std::string, std::function<double(double, double)>>> linkages_;
  std::map<std::string, std::tuple<functionArgs...>> argsMap_;
};

template <class A, typename... functionArgs>
class BaseGraph // not as general as it could be; works somewhat like a directed graph
{
public:
  BaseGraph() {}
  BaseGraph(const BaseGraph &in) : edges_(in.getEdges()), nodes_(in.getNodes()), endNodesIndices_(in.getEndNodesIndices()), adjacencyMap_(in.getAdjacencyMap()) {}
  ~BaseGraph() = default;

  std::vector<Edge<A, functionArgs...>> getEdges() const { return edges_; }
  std::vector<std::shared_ptr<Node<A>>> getNodes() const { return nodes_; }
  std::set<uint> getEndNodesIndices() const { return endNodesIndices_; }
  std::multimap<uint, std::pair<uint, uint>> getAdjacencyMap() const { return adjacencyMap_; }
  uint findInternalIDOfNode(Node<A> *) const;
  uint findInternalIDOfEndNode(Node<A> *) const;
  std::vector<uint> getIDsOfNode(std::shared_ptr<Node<A>>) const;
  void add(std::shared_ptr<Node<A>>, Edge<A, functionArgs...>, std::shared_ptr<Node<A>>);
  void addNew(std::shared_ptr<Node<A>>, Edge<A, functionArgs...>, std::shared_ptr<Node<A>>);
  void clear();
  template <typename... Args>
  void print(A, Args &&...);
  void print();
  void operator=(const BaseGraph<A> &);
  bool operator==(const BaseGraph<A> &) const;

protected:
  void findEndNodes(uint);
  void findEndNodes();
  void copy(const BaseGraph *);
  std::vector<Edge<A, functionArgs...>> edges_;
  std::vector<std::shared_ptr<Node<A>>> nodes_;
  std::set<uint> endNodesIndices_;                          // can be empty
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
  template <typename... functionArgs>
  void add(std::shared_ptr<Node<A>>, std::shared_ptr<Edge<A, functionArgs...>>, std::shared_ptr<Node<A>>);
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
  std::vector<std::shared_ptr<Node<A>>> pathNodes_;
  std::vector<std::shared_ptr<Edge<A>>> pathEdges_;
  std::map<uint, std::pair<uint, uint>> pathAdjacencyMap_;
};

template <class A, typename... functionArgs>
class DecisionTree : public BaseGraph<A, functionArgs...>
{
public:
  DecisionTree() : BaseGraph<A, functionArgs...>() {}
  DecisionTree(const DecisionTree &in) : BaseGraph<A, functionArgs...>(in), allPaths_(in.getAllPaths()) {}
  ~DecisionTree() = default;

  Node<A> getRootNode() const;
  uint getRootNodeIndex() const;
  void addToEndNode(Node<A> inNode, DecisionTree<A, functionArgs...> inTree, bool refreshEndNodes = true);
  void connectGraphs(std::vector<std::shared_ptr<Node<A>>>, std::vector<Edge<A, functionArgs...>> *, DecisionTree<A, functionArgs...>);
  void connectGraphs(uint, std::shared_ptr<Edge<A, functionArgs...>>, DecisionTree<A, functionArgs...>);
  void addAndCopyGraph(std::vector<std::shared_ptr<Node<A>>> *, std::vector<Edge<A, functionArgs...>> *, DecisionTree<A, functionArgs...>);
  void setInitialData(std::string inName, double inValue) { BaseGraph<A>::nodes_[getRootNodeIndex()]->addData(inName, inValue); }
  std::vector<Path<A>> getAllPaths() const { return allPaths_; }
  uint getNAllPaths() const { return allPaths_.size(); }
  void findAllPaths(uint = 0);
  double getExpectationValue(std::string, std::string);
  void updateNodeValues(uint, A);
  void clear();
  void operator=(const DecisionTree<A, functionArgs...> &);
  bool operator==(const DecisionTree<A, functionArgs...> &) const;

private:
  void DFS(uint start, uint finish, Path<A>);
  void DFS(Node<A> start, Node<A> finish, Path<A> *); // Depth-First-Search
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

template <class A>
void Node<A>::print() const
{
  std::cout << name_ << std::endl;
  std::cout << externalID_ << std::endl;
  for (auto it = data_.begin(); it != data_.end(); it++)
  {
    std::cout << it->first << ": " << it->second << std::endl;
  }
}

// definitions for Edge
template <class A, typename... Args>
void Edge<A, Args...>::changeNodes(std::shared_ptr<Node<A>> node1, std::shared_ptr<Node<A>> node2)
{
  nodes_.first = std::move(node1);
  nodes_.second = std::move(node2);
}

template <class A, typename... Args>
template <typename T>
void Edge<A, Args...>::insertFunction(std::string s1, T f1)
{
  auto tt = std::type_index(typeid(f1));
  functions_.insert(std::make_pair(s1, std::make_pair((voidFunctionType)f1, tt)));
}

template <class A, typename... functionArgs>
template <typename T>
void Edge<A, functionArgs...>::insertClassFunctionVoid(std::string s1, T f1)
{
  auto tt = std::type_index(typeid(f1));
  classFunctions_.insert(std::make_pair((std::string)s1, std::make_pair((voidClassFunctionType<A>)f1, tt)));
}

template <class A, typename... functionArgs>
template <typename T>
void Edge<A, functionArgs...>::insertClassFunction(std::string s1, T f1, functionArgs &&...args)
{
  auto tt = std::type_index(typeid(f1));
  classFunctions_.insert(std::make_pair((std::string)s1, std::make_pair((voidClassFunctionType<A>)f1, tt)));
  argsMap_[s1] = std::make_tuple(std::move(args)...);
}

template <class A, typename... functionArgs>
template <typename T, typename... Args>
T Edge<A, functionArgs...>::searchAndCall(std::string s1, Args &&...args)
{
  auto mapIter = functions_.find(s1);
  auto mapVal = mapIter->second;

  auto typeCastedFun = (T(*)(Args...))(mapVal.first);

  // compare if the types are equal or not
  assert(mapVal.second == std::type_index(typeid(typeCastedFun)));
  return typeCastedFun(std::forward<Args>(args)...);
}

template <class A, typename... functionArgs>
template <typename T, typename... Args>
T Edge<A, functionArgs...>::searchAndCall(A a, std::string s1, Args &&...args)
{
  auto mapIter = classFunctions_.find(s1);
  auto mapVal = mapIter->second;

  auto typeCastedFun = (T(A::*)(Args...))(mapVal.first);

  assert(mapVal.second == std::type_index(typeid(typeCastedFun)));
  return (a.*typeCastedFun)(std::forward<Args>(args)...);
}

template <class A, typename... Args>
void Edge<A, Args...>::operator=(const Edge<A> &inEdge)
{
  functions_ = inEdge.getFunctions();
  classFunctions_ = inEdge.getClassFunctions();
  data_ = inEdge.getData();
  nodes_ = inEdge.getNodes();
  linkages_ = inEdge.getLinkages();
}

template <class A, typename... Args>
bool Edge<A, Args...>::operator==(const Edge<A> &inEdge) const
{
  bool nodes = (nodes_ == inEdge.getNodes());
  bool data = (data_ == inEdge.getData());
  bool functions = (functions_ == inEdge.getFunctions());
  bool classFunctions = (classFunctions_ == inEdge.getClassFunctions());
  bool links = (linkages_ == inEdge.getLinkages());
  if (nodes && data && functions && classFunctions && links)
  {
    return true;
  }
  return false;
}

template <class A, typename... functionArgs>
template <typename... Args>
void Edge<A, functionArgs...>::updateData(A inObject, Args &&...args)
{
  for (auto it : nodes_.first->getData())
  {
    std::cout << it.first << ", " << it.second << std::endl;
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    if (linkages_.find(it.first) == linkages_.end())
    {
      continue;
    }
    std::string functionName = linkages_.at(it.first).first;
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    if (classFunctions_.find(functionName) != classFunctions_.end())
    {
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      nodes_.second->addData(it.first, (linkages_.at(it.first).second)(it.second, searchAndCall<double>(inObject, functionName, std::forward<Args>(args)...)));
    }
    else if (functions_.find(functionName) != functions_.end())
    {
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      nodes_.second->addData(it.first, (linkages_.at(it.first).second)(it.second, searchAndCall<double>(functionName, std::forward<Args>(args)...)));
    }
    else
    {
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      nodes_.second->addData(it.first, it.second);
    }
  }
  std::cout << __FILE__ << ": " << __LINE__ << std::endl;
}

template <class A, typename... functionArgs>
void Edge<A, functionArgs...>::updateData(A inObject)
{
  for (auto it : nodes_.first->getData())
  {
    std::cout << it.first << ", " << it.second << std::endl;
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    if (linkages_.find(it.first) == linkages_.end())
    {
      continue;
    }
    std::string functionName = linkages_.at(it.first).first;
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    if (classFunctions_.find(functionName) != classFunctions_.end())
    {
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      if (argsMap_.find(functionName) != argsMap_.end())
      {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        double args2 = invoke_hpp::apply([&](auto &&...args)
                                         { return searchAndCall<double, functionArgs...>(inObject, functionName, std::forward<functionArgs>(args)...); },
                                         argsMap_.at(functionName));
        nodes_.second->addData(it.first, (linkages_.at(it.first).second)(it.second, args2));
      }
      else
      {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        nodes_.second->addData(it.first, (linkages_.at(it.first).second)(it.second, searchAndCall<double>(inObject, functionName)));
      }
    }
    else if (functions_.find(functionName) != functions_.end())
    {
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      // if (argsMap_.find(functionName) != argsMap_.end())
      // {
      //   nodes_.second->addData(it.first, (linkages_.at(it.first).second)(it.second, std::apply([&](auto &&...args)
      //                                                                                          { searchAndCall<double>(functionName, std::forward<functionArgs>(args)...); },
      //                                                                                          argsMap_.at(functionName))));
      // }
      // else
      // {
      nodes_.second->addData(it.first, (linkages_.at(it.first).second)(it.second, searchAndCall<double>(functionName)));
      // }
    }
    else
    {
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      nodes_.second->addData(it.first, it.second);
    }
  }
  std::cout << __FILE__ << ": " << __LINE__ << std::endl;
}

template <class A, typename... functionArgs>
template <typename... Args>
void Edge<A, functionArgs...>::updateData(Args &&...args)
{
  for (auto it = nodes_.first->getData().begin(); it != nodes_.first.getData().end(); it++)
  {
    std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    std::cout << it->first << std::endl;
    std::string functionName = linkages_.at(it->first).first;
    std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    if (functions_.find(functionName) != functions_.end())
    {
      std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      nodes_.second->addData(it->first, (linkages_.at(it->first).second)(it->second, searchAndCall<double>(functionName, std::forward<Args>(args)...)));
      std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    }
    else
    {
      nodes_.second->addData(it->first, it->second);
    }
  }
}

// definitions for BaseGraph
template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::add(std::shared_ptr<Node<A>> inNode1, Edge<A, functionArgs...> inEdge, std::shared_ptr<Node<A>> inNode2)
{
  bool addNode1 = true;
  bool addNode2 = true;
  std::shared_ptr<Node<A>> node1 = inNode1;
  std::shared_ptr<Node<A>> node2 = inNode2;
  uint node1Index = 0;
  uint node2Index = 0;
  for (uint i = 0; i < nodes_.size(); i++)
  {
    if (nodes_[i] == inNode1)
    {
      addNode1 = false;
      node1 = nodes_[i];
      node1Index = i;
    }
    else if (nodes_[i] == inNode2)
    {
      addNode2 = false;
      node2 = nodes_[i];
      node2Index = i;
    }
    if (!addNode1 && !addNode2)
    {
      break;
    }
  }
  if (addNode1)
  {
    nodes_.push_back(node1);
    node1Index = nodes_.size() - 1;
    // node1 = nodes_.back();
  }
  if (addNode2)
  {
    nodes_.push_back(node2);
    node2Index = nodes_.size() - 1;
    // node2 = nodes_.back();
  }
  Edge<A, functionArgs...> inEdgeCopy(inEdge);
  inEdgeCopy.changeNodes(node1, node2);
  edges_.push_back(inEdgeCopy);
  adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(node1Index, std::pair<uint, uint>(node2Index, edges_.size() - 1)));
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::addNew(std::shared_ptr<Node<A>> inNode1, Edge<A, functionArgs...> inEdge, std::shared_ptr<Node<A>> inNode2)
{
  std::shared_ptr<Node<A>> node1;
  uint node1Index = 0;
  uint edgeIndex = 0;
  for (int i = nodes_.size() - 1; i >= 0; i--)
  {
    if (nodes_[i] == inNode1)
    {
      node1 = nodes_[i];
      node1Index = i;
      break;
    }
  }
  nodes_.push_back(inNode2);
  std::shared_ptr<Node<A>> node2 = nodes_.back();

  Edge<A> inEdgeCopy(inEdge);
  inEdgeCopy.changeNodes(node1, node2);
  edges_.push_back(inEdgeCopy);
  adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(node1Index, std::pair<uint, uint>(nodes_.size() - 1, edges_.size() - 1)));
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::findEndNodes(uint startID)
{
  if (startID == 0)
  {
    endNodesIndices_.clear();
  }
  if (adjacencyMap_.find(startID) != adjacencyMap_.end())
  {
    std::pair<std::multimap<uint, std::pair<uint, uint>>::iterator, std::multimap<uint, std::pair<uint, uint>>::iterator> range = adjacencyMap_.equal_range(startID);
    for (auto it = range.first; it != range.second; it++)
    {
      findEndNodes(it->second.first);
    }
  }
  else
  {
    endNodesIndices_.insert(startID);
  }
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::findEndNodes()
{
  endNodesIndices_.clear();
  for (uint i = 0; i < nodes_.size(); i++)
  {
    if (adjacencyMap_.find(i) == adjacencyMap_.end())
    {
      endNodesIndices_.insert(i);
    }
  }
}

template <class A, typename... functionArgs>
uint BaseGraph<A, functionArgs...>::findInternalIDOfNode(Node<A> *inNode) const
{
  uint nFound = 0;
  uint result = 0;
  for (Node<A> node : nodes_)
  {
    if (node.getName() != inNode->getName() || node.getID() != inNode->getID() || node.getClassFunctions() != inNode->getClassFunctions() || node.getFunctions() != inNode->getFunctions())
    {
      continue;
    }
    bool found = true;
    for (auto it = node.getData().begin; it != node.getData().end(); it++)
    {
      if (it->first == "ID")
      {
        continue;
      }
      if (it->second != inNode->getData()[it->first])
      {
        found = false;
        break;
      }
    }
    if (found)
    {
      nFound++;
      std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      result = node.getData().at("ID");
      std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    }
  }
  if (nFound > 1)
  {
    std::cout << "Multiple nodes with same ID found." << std::endl;
    abort();
  }
  return result;
}

template <class A, typename... functionArgs>
uint BaseGraph<A, functionArgs...>::findInternalIDOfEndNode(Node<A> *inNode) const
{
  uint nFound = 0;
  uint result = 0;
  for (std::shared_ptr<Node<A>> node : endNodesIndices_)
  {
    if (node->getName() != inNode->getName() || node->getID() != inNode->getID() || node->getClassFunctions() != inNode->getClassFunctions() || node->getFunctions() != inNode->getFunctions())
    {
      continue;
    }
    bool found = true;
    for (auto it = node->getData().begin(); it != node->getData().end(); it++)
    {
      if (it->first == "ID")
      {
        continue;
      }
      if (it->second != inNode->getData()[it->first])
      {
        found = false;
        break;
      }
    }
    if (found)
    {
      nFound++;
      std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      result = node->getData().at("ID");
      std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    }
  }
  if (nFound > 1)
  {
    std::cout << "Multiple end nodes with same ID found." << std::endl;
    abort();
  }
  return result;
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::copy(const BaseGraph<A, functionArgs...> *in)
{
  edges_ = in->getEdges();
  nodes_ = in->getNodes();
  endNodesIndices_ = in->getEndNodesIndices();
  adjacencyMap_ = in->getAdjacencyMap();
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::clear()
{
  nodes_.clear();
  edges_.clear();
  endNodesIndices_.clear();
  adjacencyMap_.clear();
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::operator=(const BaseGraph<A> &in)
{
  nodes_ = in.getNodes();
  edges_ = in.getEdges();
  endNodesIndices_ = in.getEndNodesIndices();
  adjacencyMap_ = in.getAdjacencyMap();
}

template <class A, typename... functionArgs>
bool BaseGraph<A, functionArgs...>::operator==(const BaseGraph<A> &in) const // TO DO: sort nodes and edges, topology
{
  bool nodes = (nodes_ == in.getNodes());
  bool edges = (edges_ == in.getEdges());
  bool endNodes = (endNodesIndices_ == in.getEndNodesIndices());
  bool adMap = (adjacencyMap_ == in.getAdjacencyMap());
  if (nodes && edges && endNodes && adMap)
  {
    return true;
  }
  return false;
}

template <class A, typename... functionArgs>
std::vector<uint> BaseGraph<A, functionArgs...>::getIDsOfNode(std::shared_ptr<Node<A>> inNode) const
{
  std::vector<uint> result;
  for (uint i = 0; i < nodes_.size(); i++)
  {
    if (*inNode == *nodes_[i])
    {
      result.push_back(i);
    }
  }
  return result;
}

template <class A, typename... functionArgs>
template <typename... Args>
void BaseGraph<A, functionArgs...>::print(A inObject, Args &&...args)
{
  for (auto it = adjacencyMap_.begin(); it != adjacencyMap_.end(); it++)
  {
    double probability = 0;
    std::string functionName1 = "Variable Probability";
    std::string functionName2 = "Constant Probability";
    std::string functionUpdateName = "Update";
    std::cout << "Edge " << it->second.second << std::endl;
    // if (edges_[it->second.second].getClassFunctions().find(functionUpdateName) != edges_[it->second.second].getClassFunctions().end())
    // {
    //   edges_[it->second.second].template searchAndCall<void>(inObject, functionUpdateName);
    // }
    // else if (edges_[it->second.second].getFunctions().find(functionUpdateName) != edges_[it->second.second].getFunctions().end())
    // {
    //   edges_[it->second.second].template searchAndCall<void>(functionUpdateName);
    // }
    if (edges_[it->second.second].getClassFunctions().count(functionName1) > 0)
    {
      probability = edges_[it->second.second].template searchAndCall<double>(inObject, functionName1, std::forward<Args>(args)...);
    }
    else if (edges_[it->second.second].getClassFunctions().count(functionName2) > 0)
    {
      probability = edges_[it->second.second].template searchAndCall<double>(inObject, functionName2);
    }
    else if (edges_[it->second.second].getFunctions().count(functionName1) > 0)
    {
      probability = edges_[it->second.second].template searchAndCall<double>(functionName1, std::forward<Args>(args)...);
    }
    else if (edges_[it->second.second].getFunctions().count(functionName2) > 0)
    {
      probability = edges_[it->second.second].template searchAndCall<double>(functionName2);
    }
    else if (edges_[it->second.second].getData().count(functionName1) > 0)
    {
      probability = edges_[it->second.second].getData()[functionName1];
    }
    else if (edges_[it->second.second].getData().count(functionName2) > 0)
    {
      probability = edges_[it->second.second].getData()[functionName2];
    }
    std::cout << nodes_[it->first].getName() << "(" << it->first << ")"
              << " --(" << probability * 100 << "%)-> " << nodes_[it->second.first].getName() << "(" << it->second.first << "), " << it->second.second << std::endl;
  }
}

template <class A, typename... functionArgs>
void BaseGraph<A, functionArgs...>::print()
{
  for (auto it = adjacencyMap_.begin(); it != adjacencyMap_.end(); it++)
  {
    std::cout << nodes_[it->first].getName() << "(" << it->first << ")"
              << " -> " << nodes_[it->second.first].getName() << "(" << it->second.first << "), " << it->second.second << std::endl;
  }
}

// definitions for DecisionTree
template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::clear()
{
  BaseGraph<A, functionArgs...>::clear();
  allPaths_.clear();
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::DFS(uint startID, uint finishID, Path<A> result)
{
  if (startID == finishID)
  {
    allPaths_.push_back(result);
  }
  else
  {
    if (BaseGraph<A, functionArgs...>::adjacencyMap_.find(startID) != BaseGraph<A, functionArgs...>::adjacencyMap_.end())
    {
      std::pair<std::multimap<uint, std::pair<uint, uint>>::iterator, std::multimap<uint, std::pair<uint, uint>>::iterator> itRange = BaseGraph<A, functionArgs...>::adjacencyMap_.equal_range(startID);
      for (auto it = itRange.first; it != itRange.second; it++)
      {
        if (it->second.first > finishID || find(BaseGraph<A>::endNodesIndices_.begin(), BaseGraph<A>::endNodesIndices_.end(), it->second.first) != BaseGraph<A>::endNodesIndices_.end()) // Tree has to be numerated from up to down, left to right
        {
          continue;
        }
        Path<A> tempResult(result);
        tempResult.add(BaseGraph<A>::nodes_[it->first],
                       std::make_shared<Edge<A, functionArgs...>>(BaseGraph<A, functionArgs...>::edges_[it->second.second]),
                       BaseGraph<A>::nodes_[it->second.first]);
        DFS(it->second.first, finishID, tempResult);
      }
    }
    else
    {
      std::cout << "Unexpected case: startID=" << startID << ", finishID=" << finishID << std::endl;
      abort();
    }
  }
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::DFS(Node<A> startNode, Node<A> finishNode, Path<A> *result) // TO DO: delete or fix this
{
  if (startNode == finishNode)
  {
    std::cout << startNode.getName() << " = " << finishNode.getName() << std::endl;
    allPaths_.push_back(*result);
    abort();
  }
  else
  {
    uint startID = 0;
    bool foundStartID = false;
    for (uint i = 0; i < BaseGraph<A, functionArgs...>::nodes_.size(); i++)
    {
      if (BaseGraph<A, functionArgs...>::nodes_[i] == startNode)
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
        range = BaseGraph<A, functionArgs...>::adjacencyMap_.equal_range(startID);
    for (std::multimap<uint, std::pair<uint, uint>>::iterator it =
             range.first;
         it != range.second; it++)
    {
      Path<A> tempResult(*result);
      tempResult.add(BaseGraph<A, functionArgs...>::nodes_[it->first],
                     BaseGraph<A, functionArgs...>::edges_[it->second.second],
                     BaseGraph<A, functionArgs...>::nodes_[it->second.first]);
      abort();
      DFS(BaseGraph<A, functionArgs...>::nodes_[it->second.first], finishNode, &tempResult);
    }
  }
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::addToEndNode(Node<A> inNode, DecisionTree<A, functionArgs...> inGraph, bool refreshEndNodes)
{
  BaseGraph<A, functionArgs...>::findEndNodes();
  uint nodeIndex = 0;
  uint foundIndices = 0;
  for (uint index : BaseGraph<A, functionArgs...>::endNodesIndices_)
  {
    if (BaseGraph<A, functionArgs...>::nodes_[index] == inNode)
    {
      foundIndices++;
      nodeIndex = index;
    }
  }
  if (foundIndices != 1)
  {
    std::cout << "Cannot connect graph to end node." << std::endl;
    for (auto it = BaseGraph<A, functionArgs...>::adjacencyMap_.begin(); it != BaseGraph<A, functionArgs...>::adjacencyMap_.end(); it++)
    {
      std::cout << it->first << " -> " << it->second.first << ", " << it->second.second << std::endl;
    }
    std::cout << "end nodes:" << std::endl;
    for (auto endNodeIndex : BaseGraph<A, functionArgs...>::endNodesIndices_)
    {
      for (uint i = 0; i < BaseGraph<A, functionArgs...>::nodes_.size(); i++)
      {
        if (BaseGraph<A, functionArgs...>::nodes_[endNodeIndex] == BaseGraph<A, functionArgs...>::nodes_[i])
        {
          std::cout << i << std::endl;
        }
      }
      std::cout << std::endl;
    }
    std::cout << "nodes size = " << BaseGraph<A, functionArgs...>::nodes_.size() << std::endl;
    std::cout << "end nodes size = " << BaseGraph<A, functionArgs...>::endNodesIndices_.size() << std::endl;
    abort();
  }
  else
  {
    uint nodesSize = BaseGraph<A, functionArgs...>::nodes_.size();
    std::vector<Node<A>> newNodes(inGraph.getNodes());
    std::vector<Edge<A, functionArgs...>> newEdges(inGraph.getEdges());
    Edge<A, functionArgs...> link(BaseGraph<A, functionArgs...>::nodes_[nodeIndex], newNodes[0]); // TO DO:get root node
    BaseGraph<A, functionArgs...>::addNew(BaseGraph<A, functionArgs...>::nodes_[nodeIndex], link, newNodes.at(0));
    BaseGraph<A, functionArgs...>::nodes_.insert(BaseGraph<A, functionArgs...>::nodes_.end(), newNodes.begin() + 1, newNodes.end());
    std::multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getAdjacencyMap());
    for (std::multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it++)
    {
      Edge<A, functionArgs...> link(newEdges[it->second.second]);
      link.changeNodes(BaseGraph<A, functionArgs...>::nodes_[it->first + nodesSize], BaseGraph<A, functionArgs...>::nodes_[it->second.first + nodesSize]);
      BaseGraph<A, functionArgs...>::edges_.push_back(link);
      BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(it->first + nodesSize, std::pair<uint, uint>(it->second.first + nodesSize, BaseGraph<A>::edges_.size() - 1)));
    }
  }
  if (refreshEndNodes)
  {
    BaseGraph<A>::endNodesIndices_.clear();
    BaseGraph<A>::findEndNodes();
  }
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::connectGraphs(std::vector<std::shared_ptr<Node<A>>> inNode, std::vector<Edge<A, functionArgs...>> *inEdges, DecisionTree<A, functionArgs...> inGraph)
{
  bool foundNode = false;
  std::vector<uint> iNodes;
  BaseGraph<A, functionArgs...>::findEndNodes();

  for (std::shared_ptr<Node<A>> node : inNode)
  {
    for (uint nodeIndex : BaseGraph<A, functionArgs...>::endNodesIndices_)
    {
      if (*node == *BaseGraph<A, functionArgs...>::nodes_[nodeIndex])
      {
        iNodes.push_back(nodeIndex);
        foundNode = true;
      }
    }
  }
  if (!foundNode)
  {
    std::cout << "Cannot connect two graphs." << std::endl;
    abort();
  }
  else
  {
    uint nodesSize = BaseGraph<A, functionArgs...>::nodes_.size();
    std::vector<Node<A>> newNodes(inGraph.getNodes());
    std::vector<Edge<A, functionArgs...>> newEdges(inGraph.getEdges());
    BaseGraph<A, functionArgs...>::nodes_.insert(BaseGraph<A>::nodes_.end(), newNodes.begin(), newNodes.end());
    uint edgeCounter = 0;
    for (uint iNode : iNodes)
    {
      Edge<A, functionArgs...> link(inEdges->at(edgeCounter));
      link.changeNodes(BaseGraph<A>::nodes_[iNode], BaseGraph<A, functionArgs...>::nodes_[nodesSize]); // TO DO: get root node
      BaseGraph<A, functionArgs...>::edges_.push_back(link);
      edgeCounter++;
      BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(iNode, std::pair<uint, uint>(nodesSize, BaseGraph<A>::edges_.size() - 1)));
    }
    std::multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getAdjacencyMap());
    for (std::multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it++)
    {
      uint node1Index = it->first + nodesSize;
      uint node2Index = it->second.first + nodesSize;
      Edge<A, functionArgs...> link(newEdges[it->second.second]);
      link.changeNodes(BaseGraph<A>::nodes_[node1Index], BaseGraph<A>::nodes_[node2Index]);
      BaseGraph<A, functionArgs...>::edges_.push_back(link);
      BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(node1Index, std::pair<uint, uint>(node2Index, BaseGraph<A>::edges_.size() - 1)));
    }
  }
  BaseGraph<A, functionArgs...>::endNodesIndices_.clear();
  BaseGraph<A, functionArgs...>::findEndNodes();
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::addAndCopyGraph(std::vector<std::shared_ptr<Node<A>>> *inNode, std::vector<Edge<A, functionArgs...>> *inEdges, DecisionTree<A, functionArgs...> inGraph)
{
  bool foundNode = false;
  std::vector<uint> iNodes;
  BaseGraph<A, functionArgs...>::findEndNodes();
  for (std::shared_ptr<Node<A>> node : *inNode)
  {
    for (uint nodeIndex : BaseGraph<A, functionArgs...>::endNodesIndices_)
    {
      if (*node == *BaseGraph<A, functionArgs...>::nodes_[nodeIndex])
      {
        iNodes.push_back(nodeIndex);
        foundNode = true;
      }
    }
  }
  if (!foundNode)
  {
    std::cout << "addAndCopyGraph: Cannot connect two graphs." << std::endl;
    abort();
  }
  else
  {
    uint edgeCounter = 0;
    std::multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getAdjacencyMap());
    for (uint iNode : iNodes)
    {
      uint nodesSize = BaseGraph<A, functionArgs...>::nodes_.size();
      std::vector<std::shared_ptr<Node<A>>> newNodes(inGraph.getNodes().size());
      for (uint i = 0; i < newNodes.size(); i++)
      {
        newNodes[i] = std::make_shared<Node<A>>(*inGraph.getNodes()[i]);
      }
      std::vector<Edge<A, functionArgs...>> newEdges(inGraph.getEdges());
      BaseGraph<A, functionArgs...>::nodes_.insert(BaseGraph<A, functionArgs...>::nodes_.end(), newNodes.begin(), newNodes.end());
      Edge<A, functionArgs...> link(inEdges->at(edgeCounter));
      link.changeNodes(BaseGraph<A, functionArgs...>::nodes_[iNode], BaseGraph<A, functionArgs...>::nodes_[nodesSize]); // TO DO: get root node
      BaseGraph<A, functionArgs...>::edges_.push_back(link);
      edgeCounter = (edgeCounter + 1 > inNode->size() - 1) ? 0 : edgeCounter + 1;
      BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(iNode, std::pair<uint, uint>(nodesSize, BaseGraph<A, functionArgs...>::edges_.size() - 1)));
      for (std::multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it++)
      {
        uint node1Index = it->first + nodesSize;
        uint node2Index = it->second.first + nodesSize;
        Edge<A, functionArgs...> link(newEdges[it->second.second]);
        link.changeNodes(BaseGraph<A, functionArgs...>::nodes_[node1Index], BaseGraph<A, functionArgs...>::nodes_[node2Index]);
        BaseGraph<A, functionArgs...>::edges_.push_back(link);
        BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(node1Index, std::pair<uint, uint>(node2Index, BaseGraph<A, functionArgs...>::edges_.size() - 1)));
      }
    }
  }
  BaseGraph<A, functionArgs...>::endNodesIndices_.clear();
  BaseGraph<A, functionArgs...>::findEndNodes();
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::connectGraphs(uint inNodeID, std::shared_ptr<Edge<A, functionArgs...>> inEdge, DecisionTree<A, functionArgs...> inGraph)
{
  uint nodesSize = BaseGraph<A, functionArgs...>::nodes_.size();
  std::vector<std::shared_ptr<Node<A>>> newNodes(inGraph.getNodes().size());
  for (uint i = 0; i < newNodes.size(); i++)
  {
    newNodes[i] = std::make_shared<Node<A>>(*inGraph.getNodes()[i]);
  }
  std::vector<Edge<A, functionArgs...>> newEdges(inGraph.getEdges());
  BaseGraph<A, functionArgs...>::nodes_.insert(BaseGraph<A, functionArgs...>::nodes_.end(), newNodes.begin(), newNodes.end());
  Edge<A, functionArgs...> link(*inEdge);
  link.changeNodes(BaseGraph<A, functionArgs...>::nodes_[inNodeID], BaseGraph<A, functionArgs...>::nodes_[nodesSize]); // TO DO: get root node
  BaseGraph<A, functionArgs...>::edges_.push_back(link);
  BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(inNodeID, std::pair<uint, uint>(nodesSize, BaseGraph<A, functionArgs...>::edges_.size() - 1)));
  std::multimap<uint, std::pair<uint, uint>> tempMap(inGraph.getAdjacencyMap());
  for (std::multimap<uint, std::pair<uint, uint>>::iterator it = tempMap.begin(); it != tempMap.end(); it++)
  {
    uint node1Index = it->first + nodesSize;
    uint node2Index = it->second.first + nodesSize;
    Edge<A, functionArgs...> newLink(newEdges[it->second.second]);
    newLink.changeNodes(BaseGraph<A, functionArgs...>::nodes_[node1Index], BaseGraph<A, functionArgs...>::nodes_[node2Index]);
    BaseGraph<A, functionArgs...>::edges_.push_back(newLink);
    BaseGraph<A, functionArgs...>::adjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(node1Index, std::pair<uint, uint>(node2Index, BaseGraph<A, functionArgs...>::edges_.size() - 1)));
  }
  BaseGraph<A, functionArgs...>::endNodesIndices_.clear();
  BaseGraph<A, functionArgs...>::findEndNodes();
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::operator=(const DecisionTree<A, functionArgs...> &in)
{
  BaseGraph<A, functionArgs...>::operator=(in);
  allPaths_ = in.getAllPaths();
}

template <class A, typename... functionArgs>
bool DecisionTree<A, functionArgs...>::operator==(const DecisionTree<A, functionArgs...> &in) const // TO DO: topology
{
  bool base = (BaseGraph<A, functionArgs...>::operator==(in));
  bool paths = (allPaths_ == in.getAllPaths());
  if (base && paths)
  {
    return true;
  }
  return false;
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::findAllPaths(uint startNodeID)
{
  BaseGraph<A, functionArgs...>::endNodesIndices_.clear();
  BaseGraph<A, functionArgs...>::findEndNodes();
  BaseGraph<A, functionArgs...>::print();
  for (uint endNodeIndex : BaseGraph<A, functionArgs...>::endNodesIndices_)
  {
    Path<A> newPath;
    DFS(startNodeID, endNodeIndex, newPath); // TO DO: get root node
  }
}

template <class A, typename... functionArgs>
Node<A> DecisionTree<A, functionArgs...>::getRootNode() const
{
  uint lastNodeIndex = 0;
  auto it = BaseGraph<A, functionArgs...>::adjacencyMap_.begin();
  uint safetyCounter = 0;
  while (it != BaseGraph<A, functionArgs...>::adjacencyMap_.end() && safetyCounter < BaseGraph<A, functionArgs...>::adjacencyMap_.size())
  {
    if (it->second.first == lastNodeIndex)
    {
      lastNodeIndex = it->first;
      it = BaseGraph<A, functionArgs...>::adjacencyMap_.begin();
      safetyCounter++;
    }
    else
    {
      it++;
    }
  }
  return BaseGraph<A, functionArgs...>::nodes_[lastNodeIndex];
}

template <class A, typename... functionArgs>
uint DecisionTree<A, functionArgs...>::getRootNodeIndex() const
{
  uint lastNodeIndex = 0;
  auto it = BaseGraph<A, functionArgs...>::adjacencyMap_.begin();
  uint safetyCounter = 0;
  while (it != BaseGraph<A, functionArgs...>::adjacencyMap_.end() && safetyCounter < BaseGraph<A, functionArgs...>::adjacencyMap_.size())
  {
    if (it->second.first == lastNodeIndex)
    {
      lastNodeIndex = it->first;
      it = BaseGraph<A, functionArgs...>::adjacencyMap_.begin();
      safetyCounter++;
    }
    else
    {
      it++;
    }
  }
  return lastNodeIndex;
}

template <class A, typename... functionArgs>
double DecisionTree<A, functionArgs...>::getExpectationValue(std::string inValueName, std::string probabilityName)
{
  BaseGraph<A, functionArgs...>::endNodesIndices_.clear();
  BaseGraph<A, functionArgs...>::findEndNodes();
  double result = 0;
  for (uint nodeIndex : BaseGraph<A, functionArgs...>::endNodesIndices_)
  {
    result += BaseGraph<A, functionArgs...>::nodes_[nodeIndex]->getData(inValueName) * BaseGraph<A, functionArgs...>::nodes_[nodeIndex]->getData(probabilityName);
  }
  return result;
}

template <class A, typename... functionArgs>
void DecisionTree<A, functionArgs...>::updateNodeValues(uint startID, A inObject)
{
  auto range = BaseGraph<A, functionArgs...>::adjacencyMap_.equal_range(startID);
  std::cout << "startID = " << startID << std::endl;
  for (auto it1 = range.first; it1 != range.second; it1++)
  {
    uint dataSize = BaseGraph<A, functionArgs...>::edges_[it1->second.second].getNodes().first->getData().size();
    assert(dataSize > 0);
    BaseGraph<A, functionArgs...>::edges_[it1->second.second].updateData(inObject);
    updateNodeValues(it1->second.first, inObject);
  }
  abort();
}

// definitions for Path
template <class A>
template <typename... functionArgs>
void Path<A>::add(std::shared_ptr<Node<A>> inNode1, std::shared_ptr<Edge<A, functionArgs...>> inEdge, std::shared_ptr<Node<A>> inNode2)
{
  uint node1Index = pathNodes_.size();
  uint node2Index = pathNodes_.size() + 1;
  if (pathNodes_.size() == 0)
  {
    pathNodes_.push_back(inNode1);
  }
  else
  {
    node1Index--;
    node2Index--;
  }
  pathNodes_.push_back(inNode2);
  pathEdges_.push_back(inEdge);
  pathAdjacencyMap_.insert(std::pair<uint, std::pair<uint, uint>>(node1Index, std::pair<uint, uint>(node2Index, pathEdges_.size() - 1)));
}

template <class A>
double Path<A>::getProbability(A inObject)
{
  std::string functionName1 = "Variable Probability";
  std::string functionName2 = "Constant Probability";
  double result = 1;
  for (auto it = BaseGraph<A>::adjacencyMap_.begin(); it != BaseGraph<A>::adjacencyMap_.end(); it++)
  {
    BaseGraph<A>::edges_[it->second.second].template searchAndCall<void>(inObject, updateName_);

    if (BaseGraph<A>::edges_[it->second.second].getClassFunctions().find(functionName2) !=
        BaseGraph<A>::edges_[it->second.second].getClassFunctions().end())
    {
      result *= BaseGraph<A>::edges_[it->second.second].template searchAndCall<double>(inObject, functionName2);
    }
    else if (BaseGraph<A>::edges_[it->second.second].getFunctions().find(functionName2) !=
             BaseGraph<A>::edges_[it->second.second].getFunctions().end())
    {
      result *= BaseGraph<A>::edges_[it->second.second].template searchAndCall<double>(functionName2);
    }
    else if (BaseGraph<A>::edges_[it->second.second].getData().find(functionName2) !=
             BaseGraph<A>::edges_[it->second.second].getData().end())
    {
      result *= BaseGraph<A>::edges_[it->second.second].getData()[functionName2];
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