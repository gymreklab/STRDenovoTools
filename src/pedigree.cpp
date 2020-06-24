/*
Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>

This file is part of STRDenovoTools.

STRDenovoTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

STRDenovoTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with STRDenovoTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "src/common.h"
#include "src/pedigree.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

PedigreeSet::PedigreeSet() {}

bool PedigreeSet::GetFamilyIndex(const std::string& familyid, size_t* family_index) {
  for (size_t i = 0; i < families.size(); i++) {
    if (families[i].get_family_id() == familyid) {
      *family_index = i;
      return true;
    }
  }
  return false;
}

bool PedigreeSet::ExtractFamilies(const std::string& famfile,
				  const std::set<std::string>& samples,
				  const int& require_num_children) {
  families.clear();

  // Read original pedigree
  PedigreeGraph pedigree(famfile);

  // Remove irrelevant samples
  pedigree.prune(samples);

  // Identify simple nuclear families
  std::vector<PedigreeGraph*> pedigree_components;
  pedigree.split_into_connected_components(pedigree_components);
  int num_others = 0;
  for (unsigned int i = 0; i < pedigree_components.size(); i++){
    if (pedigree_components[i]->is_nuclear_family()) {
      NuclearFamily nf = pedigree_components[i]->convert_to_nuclear_family();
      if (nf.num_children() >= require_num_children) {
	families.push_back(nf);
      }
    }
    else {
      num_others++;
    }
    delete pedigree_components[i];
  }
  return true;
}

void PedigreeSet::PrintStatus() {
  int unaffected = 0;
  int affected = 0;
  int unknown = 0;
  stringstream ss;
  for (auto fam_iter = families.begin(); fam_iter != families.end(); fam_iter++) {
    std::vector<int> children_status = fam_iter->GetChildrenStatus();
    std::vector<std::string> children = fam_iter->get_children();
    ss << "Family=" << fam_iter->get_family_id() << " "
       << "Mother=" << fam_iter->get_mother() << " "
       << "Father=" << fam_iter->get_father() << " ";
    for (auto child_iter = children.begin(); 
	 child_iter != children.end(); child_iter++) {
      ss << "Child=" << *child_iter << " " ;
    }
    //PrintMessageDieOnError(ss.str(), M_PROGRESS);
    for (auto child_iter = children_status.begin();
	 child_iter != children_status.end();
	 child_iter++) {
      if ((*child_iter)==PT_CONTROL) {
	unaffected++;
      } else if ((*child_iter)==PT_CASE) {
	affected++;
      } else {
	unknown++;
      }
    }
  }
  ss.str("");
  ss << "PedigreeSet has " << families.size() << " nuclear families with STR data. " << endl
     << "  Unaffected children: " << unaffected << endl
     << "  Affected children: " << affected << endl
     << "  Unknown children: " << unknown;
  PrintMessageDieOnError(ss.str(), M_PROGRESS);
}

PedigreeSet::~PedigreeSet() {}

NuclearFamily::NuclearFamily(const std::string& family_id,
			     const std::string& mother, const std::string& father,
			     const std::vector<std::string>& children,
			     const std::vector<int>& children_status,
			     const std::vector<int>& children_sex) {
  family_id_ = family_id;
  mother_ = mother;
  father_ = father;
  children_ = children;
  children_status_ = children_status;
  children_sex_ = children_sex;
}

const int NuclearFamily::get_child_sex(const std::string& child_id) {
  size_t pos = find(children_.begin(), children_.end(), child_id) - children_.begin();
  return children_sex_[pos];
}

const int NuclearFamily::get_child_phenotype(const std::string& child_id) {
  size_t pos = find(children_.begin(), children_.end(), child_id) - children_.begin();
  return children_status_[pos];
}

NuclearFamily::~NuclearFamily() {}

void PedigreeGraph::init_no_ancestors() {
  no_ancestors_.clear();
  for (int i = 0; i < nodes_.size(); i++)
    if (!nodes_[i]->has_mother() && !nodes_[i]->has_father())
      no_ancestors_.push_back(nodes_[i]);
}

void PedigreeGraph::init_no_descendants() {
  no_descendants_.clear();
  for (int i = 0; i < nodes_.size(); i++)
    if (nodes_[i]->get_children().size() == 0)
      no_descendants_.push_back(nodes_[i]);
}

bool PedigreeGraph::topological_sort(std::vector<PedigreeNode*>& nodes){
  no_ancestors_.clear();
  no_descendants_.clear();
  nodes_.clear();
  
  std::map<PedigreeNode*, int> parent_counts;
  std::vector<PedigreeNode*>   sources;
  for (int i = 0; i < nodes.size(); i++){
    int count = nodes[i]->has_mother() + nodes[i]->has_father();
    if (count == 0)
      sources.push_back(nodes[i]);
    else
      parent_counts[nodes[i]] = count;
  }

  while (sources.size() != 0){
    PedigreeNode* source = sources.back();
    std::vector<PedigreeNode*>& children = source->get_children();
    nodes_.push_back(source);
    sources.pop_back();

    for (auto child_iter = children.begin(); child_iter != children.end(); child_iter++){
      auto count_iter = parent_counts.find(*child_iter);
      if (count_iter == parent_counts.end()){
	source->print(std::cerr);
	(*child_iter)->print(std::cerr);
	PrintMessageDieOnError("Logical error in topological_sort() for parent " + source->get_name() + " and child " + (*child_iter)->get_name(), M_ERROR);
      }
      else if (count_iter->second == 1){
	sources.push_back(*child_iter);
	parent_counts.erase(count_iter);
      }
      else
	count_iter->second -= 1;
    }
  }
  return parent_counts.size() == 0; // Only a DAG if no unprocessed individuals are left
}

bool PedigreeGraph::build(const std::string& filename) {
  std::ifstream input(filename.c_str());
  if (!input.is_open())
    PrintMessageDieOnError("Failed to open pedigree file " + filename, M_ERROR);
  
  std::map<std::string, PedigreeNode*> samples;
  std::vector<PedigreeNode*> nodes;
  std::string line;
  while (std::getline(input, line)){
    std::istringstream iss(line);
    std::string family, child, father, mother, sex, phenotype;
    if(! (iss >> family >> child >> father >> mother >> sex >> phenotype)) {
      PrintMessageDieOnError("Improperly formated .ped pedigree file " + filename, M_ERROR);
    }

    if (child.compare("0") == 0)
      PrintMessageDieOnError("Invalid individual id " + child, M_ERROR);

    // Get phenotype code
    int phenotype_code;
    if (phenotype == "1") {
      phenotype_code = PT_CONTROL;
    } else if (phenotype == "2") {
      phenotype_code = PT_CASE;
    } else {
      phenotype_code = PT_MISSING;
    }

    // Get sex code
    int sex_code;
    if (sex == "1") {
      sex_code = SEX_MALE;
    } else if (sex == "2") {
      sex_code = SEX_FEMALE;
    } else {
      sex_code = SEX_MISSING;
    }

    // Create new nodes for any previously unseen samples that have an identifier other than 0
    if (samples.find(child) == samples.end()){
      PedigreeNode* new_node = new PedigreeNode(child, family, phenotype_code, sex_code);
      nodes.push_back(new_node);
      samples[child] = new_node;
    }
    if (mother.compare("0") != 0 && samples.find(mother) == samples.end()){
      PedigreeNode* new_node = new PedigreeNode(mother, family, phenotype_code, sex_code);
      samples[mother] = new_node;
      nodes.push_back(new_node);
    }
    if (father.compare("0") != 0 && samples.find(father) == samples.end()){
      PedigreeNode* new_node = new PedigreeNode(father, family, phenotype_code, sex_code);
      samples[father] = new_node;
      nodes.push_back(new_node);
    }
    
    // Store relationships in node instance
    PedigreeNode* child_node  = samples.find(child)->second;
    PedigreeNode* mother_node = (mother.compare("0") == 0 ? NULL : samples.find(mother)->second);
    PedigreeNode* father_node = (father.compare("0") == 0 ? NULL : samples.find(father)->second);

    // Ensure that the family ids are consistent
    if (child_node->get_family().compare(family) != 0)
      PrintMessageDieOnError("Inconsistent family IDs detected in FAM file for sample " + child, M_ERROR);
    if (mother_node != NULL && mother_node->get_family().compare(family) != 0)
      PrintMessageDieOnError("Inconsistent family IDs detected in FAM file for sample " + mother, M_ERROR);
    if (father_node != NULL && father_node->get_family().compare(family) != 0)
      PrintMessageDieOnError("Inconsistent family IDs detected in FAM file for sample " + father, M_ERROR);

    child_node->set_mother(mother_node);
    child_node->set_father(father_node);
    if (mother_node != NULL) mother_node->add_child(child_node);
    if (father_node != NULL) father_node->add_child(child_node);
  }
  input.close();
  
  // Sort nodes in pedigree graph topologically
  return topological_sort(nodes);
}

void PedigreeGraph::prune(const std::set<std::string>& sample_set){
  // Determine if each node has an upstream requested sample
  std::map<PedigreeNode*, bool> upstream_status;
  for (int i = 0; i < nodes_.size(); i++){
    bool has_upstream = sample_set.find(nodes_[i]->get_name()) != sample_set.end();
    has_upstream     |= (nodes_[i]->has_father() && upstream_status[nodes_[i]->get_father()]);
    has_upstream     |= (nodes_[i]->has_mother() && upstream_status[nodes_[i]->get_mother()]);
    upstream_status[nodes_[i]] = has_upstream;
  }
  
  // Determine if each node has a downstream requested sample
  std::map<PedigreeNode*, bool> downstream_status;
  for (int i = nodes_.size()-1; i >= 0; i--) {
    bool has_downstream = sample_set.find(nodes_[i]->get_name()) != sample_set.end();
    for(auto iter = nodes_[i]->get_children().begin(); iter != nodes_[i]->get_children().end(); iter++)
      has_downstream |= downstream_status[(*iter)];
    downstream_status[nodes_[i]] = has_downstream;
  }

  // Determine if nodes have a requested sample both above and below
  // If not, mark them for removal
  std::map<PedigreeNode*, bool> removal_status;
  for (int i = 0; i < nodes_.size(); i++)
    removal_status[nodes_[i]] = (!upstream_status[nodes_[i]] || !downstream_status[nodes_[i]]);
      
  // Remove and modify nodes member data accordingly
  int insert_index = 0;
  for (int i = 0; i < nodes_.size(); i++){
    if (removal_status[nodes_[i]])
      delete nodes_[i];
    else {
      if (nodes_[i]->has_father() && removal_status[nodes_[i]->get_father()])
	nodes_[i]->set_father(NULL);
      if (nodes_[i]->has_mother() && removal_status[nodes_[i]->get_mother()])
	nodes_[i]->set_mother(NULL);
      int child_ins_index = 0;
      std::vector<PedigreeNode*>& children = nodes_[i]->get_children();;
      for (int j = 0; j < children.size(); j++)
	if (!removal_status[children[j]])
	  children[child_ins_index++] = children[j];
      children.resize(child_ins_index);
      nodes_[insert_index++] = nodes_[i];
    }
  }
  nodes_.resize(insert_index);

  // Rebuild internal structures
  no_ancestors_.clear();
  no_descendants_.clear();
  init_no_ancestors();
  init_no_descendants();
}

bool PedigreeGraph::build_subgraph(std::vector<PedigreeNode*>& subgraph_nodes){
  std::map<std::string, PedigreeNode*> samples;
  std::vector<PedigreeNode*> nodes;
  for (auto node_iter = subgraph_nodes.begin(); node_iter != subgraph_nodes.end(); node_iter++){
    PedigreeNode* child_node;
    PedigreeNode* mother_node = NULL;
    PedigreeNode* father_node = NULL;
    std::string child  = (*node_iter)->get_name();
    std::string family = (*node_iter)->get_family();
    int status = (*node_iter)->get_status();
    int sex = (*node_iter)->get_sex();

    // Create new nodes for any previously unseen samples
    if (samples.find(child) == samples.end()){
      child_node = new PedigreeNode(child, family, status, sex);
      samples[child] = child_node;
      nodes.push_back(child_node);
    }
    else
      child_node = samples[child];

    if ((*node_iter)->has_mother()){
      std::string mother = (*node_iter)->get_mother()->get_name();
      int mother_status = (*node_iter)->get_mother()->get_status();
      if (samples.find(mother) == samples.end()){
	mother_node = new PedigreeNode(mother, family, mother_status, SEX_FEMALE);
	nodes.push_back(mother_node);
	samples[mother] = mother_node;
      }
      else
	mother_node = samples[mother];
    }

    if ((*node_iter)->has_father()){
      std::string father = (*node_iter)->get_father()->get_name();
      int father_status = (*node_iter)->get_father()->get_status();
      if (samples.find(father) == samples.end()){
	father_node = new PedigreeNode(father, family, father_status, SEX_MALE);
	nodes.push_back(father_node);
	samples[father] = father_node;
      }
      else
	father_node = samples[father];
    }

    // Ensure that the family ids are consistent
    if (child_node->get_family().compare(family) != 0)
      PrintMessageDieOnError("Inconsistent family IDs detected in FAM file for sample " + child, M_ERROR);
    if (mother_node != NULL && mother_node->get_family().compare(family) != 0)
      PrintMessageDieOnError("Inconsistent family IDs detected in FAM file for sample " + mother_node->get_name(), M_ERROR);
    if (father_node != NULL && father_node->get_family().compare(family) != 0)
      PrintMessageDieOnError("Inconsistent family IDs detected in FAM file for sample " + father_node->get_name(), M_ERROR);

    // Store relationships in node instance
    child_node->set_mother(mother_node);
    child_node->set_father(father_node);
    if (mother_node != NULL) mother_node->add_child(child_node);
    if (father_node != NULL) father_node->add_child(child_node);
  }
  return topological_sort(nodes);
}

void PedigreeGraph::split_into_connected_components(std::vector<PedigreeGraph*>& components){
  assert(components.size() == 0);

  // Determine the component to which each node belongs
  std::set<PedigreeNode*> visited;
  std::vector< std::vector<PedigreeNode*> > component_nodes;
  for (auto node_iter = nodes_.begin(); node_iter != nodes_.end(); node_iter++){
    if (visited.find(*node_iter) != visited.end())
      continue;

    component_nodes.push_back(std::vector<PedigreeNode*>());
    std::vector<PedigreeNode*> to_process(1, *node_iter);
    while (to_process.size() != 0){
      PedigreeNode* seed = to_process.back();
      to_process.pop_back();
      if (visited.find(seed) != visited.end())
	continue;

      visited.insert(seed);
      component_nodes.back().push_back(seed);

      if (seed->has_mother() && (visited.find(seed->get_mother()) == visited.end()))
	to_process.push_back(seed->get_mother());
      if (seed->has_father() && (visited.find(seed->get_father()) == visited.end()))
	to_process.push_back(seed->get_father());

      std::vector<PedigreeNode*>& children = seed->get_children();
      for (auto child_iter = children.begin(); child_iter != children.end(); child_iter++)
	if (visited.find(*child_iter) == visited.end())
	  to_process.push_back(*child_iter);
    }
  }
  assert(visited.size() == nodes_.size());

  // Construct individual graphs for each component
  for (int i = 0; i < component_nodes.size(); i++)
    components.push_back(new PedigreeGraph(component_nodes[i]));
}

bool PedigreeGraph::is_nuclear_family() const{
  if (no_ancestors_.size() != 2)
    return false;
  if (no_descendants_.size() == 0)
    return false;
  if (no_ancestors_.size() + no_descendants_.size() != nodes_.size())
    return false;

  const std::string& p1 = no_ancestors_[0]->get_name();
  const std::string& p2 = no_ancestors_[1]->get_name();
  for (auto node_iter = no_descendants_.begin(); node_iter != no_descendants_.end(); node_iter++){
    if ((!(*node_iter)->has_mother()) || (!(*node_iter)->has_father()))
      return false;

    const std::string& mother = (*node_iter)->get_mother()->get_name();
    const std::string& father = (*node_iter)->get_father()->get_name();
    if ((mother.compare(p1) != 0) || (father.compare(p2) != 0))
      if ((mother.compare(p2) != 0) || (father.compare(p1) != 0))
	return false;
  }
  return true;
}

NuclearFamily PedigreeGraph::convert_to_nuclear_family() const {
  assert(is_nuclear_family());
  std::string mother = no_descendants_[0]->get_mother()->get_name();
  std::string father = no_descendants_[0]->get_father()->get_name();
  std::vector<std::string> children;
  std::vector<int> children_sex;
  std::vector<int> children_status;
  for (auto node_iter = no_descendants_.begin(); node_iter != no_descendants_.end(); node_iter++) {
    children.push_back((*node_iter)->get_name());
    children_status.push_back((*node_iter)->get_status());
    children_sex.push_back((*node_iter)->get_sex());
  }
  return NuclearFamily(no_descendants_[0]->get_family(), mother, father, children, children_status, children_sex);
}
