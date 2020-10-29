/*
This file is modified from Thomas Willems's HipSTR: https://github.com/tfwillems/HipSTR
*/


#ifndef SRC_PEDIGREE_H__
#define SRC_PEDIGREE_H__

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common.h"

enum PHENOTYPE {
  PT_CONTROL = 1,
  PT_CASE = 2,
  PT_MISSING = 0
};

enum SEX {
  SEX_MALE = 1,
  SEX_FEMALE = 2,
  SEX_MISSING = 0
};

class NuclearFamily {
 public:
  NuclearFamily(const std::string& family_id,
		const std::string& mother, const std::string& father,
		const std::vector<std::string>& children,
		const std::vector<int>& children_status,
		const std::vector<int>& children_sex);
  virtual ~NuclearFamily();

  const int get_child_phenotype(const std::string& child_id);
  const int get_child_sex(const std::string& child_id);
  const std::vector<int>& GetChildrenStatus() const {return children_status_;}
  const std::vector<int>& GetChildrenSex() const {return children_sex_;}
  const std::string& get_family_id() const { return family_id_; }
  const std::string& get_mother()    const { return mother_; }
  const std::string& get_father()    const { return father_; }
  const int size()                   const { return 2 + children_.size(); }
  const int num_children()           const { return children_.size();      }
  const std::vector<std::string>& get_children() const { return children_; }

 private:
  std::string family_id_;
  std::string mother_, father_;
  std::vector<std::string> children_;  
  std::vector<int> children_status_;
  std::vector<int> children_sex_;
};

class PedigreeSet {
 public:
  PedigreeSet();
  virtual ~PedigreeSet();

  // Get family with this ID
  bool GetFamilyIndex(const std::string& familyid, size_t* family_index);

  // Load families from .fam file
  bool ExtractFamilies(const std::string& famfile,
		       const std::set<std::string>& samples,
		       const int& require_num_children);
  // Print summary of loaded families
  void PrintStatus();

  const std::vector<NuclearFamily> get_families() { return families; }

 private:
  std::vector<NuclearFamily> families;
};

class PedigreeNode {
 private:
  std::string name_;
  PedigreeNode* mother_;
  PedigreeNode* father_;
  std::vector<PedigreeNode*> children_;
  std::string family_id_;
  int status_;
  int sex_;

 public:
 PedigreeNode(const std::string& name, const std::string& family_id, const int& status, const int& sex)
   : name_(name), family_id_(family_id), status_(status), sex_(sex) {
    mother_ = NULL;
    father_ = NULL;
  }

  ~PedigreeNode(){
    children_.clear();
  }

  bool has_mother()                 const { return mother_ != NULL;  }
  bool has_father()                 const { return father_ != NULL;  }
  PedigreeNode* get_mother()        const { return mother_;          }
  PedigreeNode* get_father()        const { return father_;          }
  const std::string& get_name()     const { return name_;            }
  const int& get_status()           const { return status_;          }
  const int& get_sex()              const { return sex_;             }
  const std::string& get_family()   const { return family_id_;       }
  std::vector<PedigreeNode*>& get_children() { return children_;     }

  void set_mother(PedigreeNode* mother) { mother_ = mother;           }
  void set_father(PedigreeNode* father) { father_ = father;           }
  void add_child (PedigreeNode* child)  { children_.push_back(child); }
  void del_child (PedigreeNode* child)  {
    auto iter = std::find(children_.begin(), children_.end(), child);;
    if (iter == children_.end())
      PrintMessageDieOnError("Can't delete child from node as it is not among the existing children", M_ERROR);
    children_.erase(iter);
  }
  
  void print(std::ostream& out) const {
    out << "NAME:"     << name_
	<< "\tFATHER:" << (father_ == NULL ? "NONE" : father_->get_name())
	<< "\tMOTHER:" << (mother_ == NULL ? "NONE" : mother_->get_name()) << std::endl; 
  }
};

class PedigreeGraph {
 private:
  // Nodes that don't have any ancestors 
  std::vector<PedigreeNode*> no_ancestors_;
  
  // Nodes that don't have any descendants
  std::vector<PedigreeNode*> no_descendants_;

  // Nodes in graph sorted in topological order
  std::vector<PedigreeNode*> nodes_;

  bool topological_sort(std::vector<PedigreeNode*>& nodes);
  bool build(const std::string& input_file);
  void init_no_ancestors();
  void init_no_descendants();
  bool build_subgraph(std::vector<PedigreeNode*>& sorted_nodes);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  PedigreeGraph(const PedigreeGraph& other);
  PedigreeGraph& operator=(const PedigreeGraph& other);
    
 public:
  PedigreeGraph(){}

  explicit PedigreeGraph(const std::string& input_file){
    bool success = build(input_file);
    if (!success)
      PrintMessageDieOnError("Supplied pedigree file " + input_file + " contains cycles", M_ERROR);
    init_no_ancestors();
    init_no_descendants();
  }

  explicit PedigreeGraph(std::vector<PedigreeNode*>& subgraph_nodes){
    if (!build_subgraph(subgraph_nodes))
      PrintMessageDieOnError("Subgraph in pedigree contains a cycle", M_ERROR);
    init_no_ancestors();
    init_no_descendants();
  }

  ~PedigreeGraph(){
    for (int i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    nodes_.clear();
    no_ancestors_.clear();
    no_descendants_.clear();
  }

  int size() const { return nodes_.size(); }

  void print(std::ostream& out) const {
    out << "Pedigree graph contains " << nodes_.size() << " nodes" << std::endl;
  }

  void prune(const std::set<std::string>& sample_set);

  void split_into_connected_components(std::vector<PedigreeGraph*>& components);

  bool is_nuclear_family() const;

  NuclearFamily convert_to_nuclear_family() const;
};

#endif  // SRC_PEDIGREE_H__
