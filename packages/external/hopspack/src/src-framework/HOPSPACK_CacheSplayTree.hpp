// $Id: HOPSPACK_CacheSplayTree.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_CacheSplayTree.hpp $ 

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This file is part of HOPSPACK.
//
// HOPSPACK is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see http://www.gnu.org/licenses/.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
//                 or Todd Plantenga (tplante@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*! \file HOPSPACK_CacheSplayTree.hpp 
    \brief Template definitions for HOPSPACK::CacheSplayTreeNode and
    HOPSPACK::CacheSplayTree.
*/

#ifndef HOPSPACK_CACHESPLAYTREE_HPP
#define HOPSPACK_CACHESPLAYTREE_HPP

namespace HOPSPACK {

// Forward declaration 
template <class Comparable> class CacheSplayTree;


/*! 
  \brief Templated splay tree node.

  Used only in HOPSPACK::CacheSplayTree.
  \see CacheSplayTree for more details.

  \author H. Alton Patrick, Summer 2000. 
*/
template <class Comparable> 
class CacheSplayTreeNode
{

private:

  //! Default Constructor
  CacheSplayTreeNode() : left(NULL), right(NULL) {};
  /*! \brief Constructor with a comparable element and the
    specification of the left and right sub-trees. */
  CacheSplayTreeNode(const Comparable& e,
                     CacheSplayTreeNode* l = NULL,
                     CacheSplayTreeNode* r = NULL)
    : element(e), left(l), right(r) {};

  //! Destructor
  ~CacheSplayTreeNode() {};

  //! The comparable element
  Comparable element;
  //! Left subtree
  CacheSplayTreeNode* left;
  //! Right subtree
  CacheSplayTreeNode* right;

  //! Friend
  friend class CacheSplayTree<Comparable>;
};

/*! 
  \brief Templated binary storage structure. 

  A templated splay tree. The template takes is based on a
  <em>Comparable</em> class that must support the following five
  operations.

  <ul>
  <li> Empty Constructor

       <tt>Comparable();</tt>

       Called indirectly by splay, which creates a CacheSplayTreeNode that
       doesn't actually have an element in it. 

  <li> Copy Constructor

       <tt>Comparable(const Comparable& source);</tt>

       The copy constructor is used only when insert is called.

  <li> Copy Data

       <tt> void copyData(const Comparable& source);</tt>

       The copy data is used only when find is called and a match is
       found in the splay tree. It is used to copy the relevant data
       from the matching entry to the input to find.
      
  <li> Not Equal

       <tt> bool operator!=(const Comparable& source);</tt>

       This is used only when find is called.

  <li> Less Than 

       <tt> bool operator<(const Comparable& source);</tt>

       This is used by insert and repeatedly by splay.

  <li> Greater Than

       <tt> bool operator>(const Comparable& source);</tt>

       This is used by insert and repeatedly by splay.

  </ul>

  Used only in HOPSPACK::CacheManager.
  \author H. Alton Patrick, Summer 2000. 
  \author Tammy Kolda
*/
template <class Comparable> 
class CacheSplayTree
{
public:

  //! Construct an empty splay tree.
  CacheSplayTree();
  //! Destruct the splay tree.
  ~CacheSplayTree();
  /*! \brief Return true if x is in the splay tree. Furthermore, call
    x.copyData() with the splay tree's matching entry (which may not
    be exactly the same) */
  bool find(Comparable& x);
  /*! \brief Insert x into the splay tree if it is not already there. 
    \retval Returns true if x is inserted into the tree, false otherwise.
  */
  bool insert(const Comparable& x);

  //! Return the number of nodes in the tree.
  int  getNumNodes (void) const;

private:


  //! By design, there is no copy constructor.
  CacheSplayTree (const CacheSplayTree &);
  //! By design, there is no assignment operator.
  CacheSplayTree & operator= (const CacheSplayTree &);

  /*! 
    \brief Re-organize the splay tree.

    If x is in the tree rooted at r, moves x to the root.  If x is not
    in the tree rooted at r, then the node placed at the root is the one
    which would have been reached directly before x in a normal binary
    search.
    
    This top-down splay is based on D. Sleator's implementation, which
    can be found at
    ftp://ftp.cs.cmu.edu/user/sleator/splaying/top-down-splay.c.
  */
  void splay(const Comparable& x, CacheSplayTreeNode<Comparable>*& r);

  //! Return true if the splay tree is empty.
  bool isEmpty();

  //! Root node of the splay tree.
  CacheSplayTreeNode<Comparable>* root;

  //! Number of nodes created for the tree.
  int  numNodes;

};


//----------------------------------------------------------------------
//  Template class definitions
//  (Included in the header because compilers need the definition when
//   they use an instance.)
//----------------------------------------------------------------------


template <class Comparable> 
CacheSplayTree<Comparable>::CacheSplayTree()
{
  root = NULL;
  numNodes = 0;
}

template <class Comparable> 
CacheSplayTree<Comparable>::~CacheSplayTree()
{
  // Keep deleting the root node until everything is gone.
  CacheSplayTreeNode<Comparable>* newroot;

  while(!isEmpty()){
    if(root->left == NULL)
      newroot = root->right;
    else{
      // This trick for joining the left and right subtrees works because 
      // splaying on an element larger than any other in the subtree (as
      // splay(root->element, root->left) does) moves the largest element
      // in the subtree to the root of the subtree.
      newroot = root->left;
      splay(root->element, newroot);
      newroot->right = root->right;
    }
    delete root;
    numNodes--;
    root = newroot;
  }
}

template <class Comparable> 
bool CacheSplayTree<Comparable>::find(Comparable& x)
{
  // Find x in the tree.  If it is present, replace x with the matching
  // node in the tree and return true.  Otherwise, return false.

  if(isEmpty())
    return false;

  splay(x, root);
  // If x is in the tree, it will be at the root.
  if(x != root->element)
    return false;
  x.copyData(root->element);
  return true;
}

template <class Comparable> 
bool CacheSplayTree<Comparable>::insert(const Comparable& x)
{
  if (isEmpty())
  {
    root = new CacheSplayTreeNode<Comparable>(x);
    numNodes++;
    return true;
  }

  splay(x, root);
  if (x < root->element)
  {
    CacheSplayTreeNode<Comparable>* newNode
        = new CacheSplayTreeNode<Comparable>(x);
    numNodes++;
    newNode->left = root->left;
    newNode->right = root;
    root->left = NULL;    
    root = newNode;
    return true;
  }
  else if (x > root->element)
  {
    CacheSplayTreeNode<Comparable>* newNode
        = new CacheSplayTreeNode<Comparable>(x);
    numNodes++;
    newNode->right = root->right;
    newNode->left = root;
    root->right = NULL;
    root = newNode;
    return true;
  }

  // else x already in tree, do nothing
  return false;
}  

template <class Comparable> 
void CacheSplayTree<Comparable>::splay(const Comparable& x,
                                       CacheSplayTreeNode<Comparable>*& r)
{
  // header's left and right branches hold, respectively, the *right* and
  // *left* subtrees of the reorganized tree.
  CacheSplayTreeNode<Comparable> header;

  // CacheSplayTreeNode being inspected in the main tree.
  CacheSplayTreeNode<Comparable> *current = r;

  // CacheSplayTreeNodes being used in the left and right subtrees.
  CacheSplayTreeNode<Comparable> *leftcur, *rightcur;  
  CacheSplayTreeNode<Comparable> *temp;

  if (r == NULL)
    return;

  leftcur = rightcur = &header;

  for (;;){
    if (x < current->element){
      if (current->left == NULL)
        break;
      if (x < current->left->element){
        // Rotate current with its left child.
        temp = current->left;
        current->left = temp->right;
        temp->right = current;
        current = temp;
        if (current->left == NULL)
          break;
      }
      // Place current in the tree of elements greater than x.
      rightcur->left = current;
      rightcur = current;
      current = current->left;
    }
    else if (x > current->element){
      if (current->right ==  NULL)
        break;
      if (x > current->right->element){
        // Rotate current with its right child.
        temp = current->right;
        current->right = temp->left;
        temp->left = current;
        current = temp;
        if (current->right == NULL)
          break;
      }
      // Place current in the tree of elements less than x.
        // Rotate current with its right child.
      leftcur->right = current;
      leftcur = current;
      current = current->right;
    }
    else
      break;
  }

  // Assemble the tree.
  leftcur->right = current->left;
  rightcur->left = current->right;
  current->left = header.right;
  current->right = header.left;
  r = current;
  
}

template <class Comparable> 
bool CacheSplayTree<Comparable>::isEmpty()
{
  return (root == NULL);
}

template <class Comparable> 
int CacheSplayTree<Comparable>::getNumNodes (void) const
{
  return( numNodes );
}

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_CACHESPLAYTREE_HPP
