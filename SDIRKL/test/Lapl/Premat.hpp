#ifndef Premat_h
#define Premat_h
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>
/////////////////////////////////////////////////////////////////////////
/// Used to help build csr or csl matrices.
/////////////////////////////////////////////////////////////////////////
using namespace std;
enum MType {LineWise,ColumnWise};
#include "PrematOrder.hpp"

template<MType T> class PrePremat
{
  //! Class PrePremat
  /*!
    Class of "Pre Pre Matrices"; indices and coefficient are stored
    in a stl map object. 
    This allows dynamic modification of matrix objects
    in a logarithmic time (B-trees).
  */

 protected:
  typedef pair<int,int> pos;
  typedef map<pos,double,cppos<T> > themap; //used to store the pre-matrix
  themap M;
  //int ordre;
  int maxIndiceLig,maxIndiceCol;
 public:
  typedef  map<pos,double,cppos<T> > mat;
  PrePremat(){//! basic constructor.
  }
  ~PrePremat(){
    //! destructor
    M.clear();
  }
  inline int ncoeffs() const {return M.size();
  //! return number of coefficients.
  }
  inline double& operator()(int i,int j)
    {
      maxIndiceLig=max(maxIndiceLig,i);
      maxIndiceCol=max(maxIndiceCol,j);
      return M[make_pair(i,j)];
    }
  inline bool exist(int i,int j)
    {
      return M.find(make_pair(i,j))!=M.end();
    }
  //! printing utility.
  //! \param k verbosity.
  void print(int k=1) const
    {
      cout<<"PrePreMat: size= "<<M.size()<<endl;
      
      if(k>0)
	for(typename mat::const_iterator i=M.begin();i!=M.end();i++)
	  cout<<'('<<i->first.first<<' '<<i->first.second<<")= "
	      <<i->second<<endl;
    }
  //! return the map/
  inline mat& Mmap()  {return M;} 

};
class PrematLineWise: public PrePremat<LineWise>
{
//! class of "Pre Matrices"
/* 
  coefficients are created in a stl map: (i,j)-> a_{ij}
  Class is derived from PrePremat objects.
*/
  using PrePremat<LineWise>::M;
  
  indices<LineWise> Indice;
  typedef PrePremat<LineWise>::themap::const_iterator mapiter;
 public:
  typedef map<pos,double,cppos<LineWise> > map_type;
  PrematLineWise(): PrePremat<LineWise>()
    {maxIndiceLig=-1;maxIndiceCol=-1;}
  ~PrematLineWise(){}
  inline void clear() {
    //! empty the stl container.
    M.clear();
  }

  int nblico() const
    {
      int count=1;
      mapiter i=M.begin();
      int i0=first(i);
      for(;i!=M.end();i++)
	if(first(i)!=i0)
	  {
	    i0=first(i);
	    ++count;
	  }
      return count;
    }

  inline int max_indice_lig() const {return maxIndiceLig;}
  inline int max_indice_col() const {return maxIndiceCol;}
  inline int first(mapiter i) const {return Indice.first(i->first);}
  inline int second(mapiter i) const {return Indice.second(i->first);}

  
  
};
class PrematColumnWise: public PrePremat<ColumnWise>
{
//! class of "Pre Matrices"
/* 
  coefficients are created in a stl map: (i,j)-> a_{ij}
  Class is derived from PrePremat objects.
*/
  using PrePremat<ColumnWise>::M;

  indices<ColumnWise> Indice;
  typedef  PrePremat<ColumnWise>::themap::const_iterator mapiter;
 public:
  typedef map<pos,double,cppos<ColumnWise> > map_type;
  PrematColumnWise(): PrePremat<ColumnWise>()
    {maxIndiceLig=-1;maxIndiceCol=-1;}
  ~PrematColumnWise(){}
  inline void clear() {
    //! empty the stl container.
    M.clear();
  }
  int nblico() const
    {
      int count=1;
      mapiter i=M.begin();
      int i0=first(i);
      for(;i!=M.end();i++)
	if(first(i)!=i0)
	  {
	    i0=first(i);
	    ++count;
	  }
      return count;
    }

  inline int max_indice_lig() const {return maxIndiceLig;}
  inline int max_indice_col() const {return maxIndiceCol;}
  inline int first(mapiter i) const {return Indice.first(i->first);}
  inline int second(mapiter i) const {return Indice.second(i->first);}

  

};
#endif
